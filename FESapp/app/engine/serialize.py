"""Convert a :class:`~.runner.RunResult` into the JSON payload the UI renders.

Kept separate from the engine proper so the simulation has no opinion about
presentation, and so a stored run can be re-rendered without re-simulating.
"""
import math

from . import burden
from . import params as P
from .pbpk import STATE_NAMES
from .safety import ORGANS

ROUND = 5


def _num(value, nd=ROUND):
    """Round, mapping non-finite values to None.

    NaN and Infinity are not valid JSON: json.dumps emits them as bare NaN/Infinity,
    which SQLite's JSON_VALID check rejects and JSON.parse in the browser refuses.
    None becomes null, which Plotly renders as a gap in the line - the honest way to
    show a trajectory that left the model's valid range.
    """
    v = float(value)
    return round(v, nd) if math.isfinite(v) else None


def _series(values, nd=ROUND):
    return [_num(v, nd) for v in values]


def _tumor(result):
    return {
        "days": [int(d) for d in result.t],
        "log10N": _series(result.log10N, 4),
        "toxicity": _series(result.T, 4),
        "concentration": _series(result.C, 4),
        "dose_by_day": _series(result.dose_by_day, 4),
    }


def serialize(result):
    cfg = result.config
    limits = cfg.limits()

    decisions = result.decisions
    cycles = []
    for c in result.cycles:
        lim = c.limiting
        why = decisions.get(c.day)
        cycles.append({
            "index": c.index,
            "day": c.day,
            "planned": _num(c.planned_dose, 4),
            "delivered": _num(c.delivered_dose, 4),
            "adjusted": c.adjusted,
            "iterations": c.iterations,
            "infeasible": c.infeasible,
            "reduction_pct": _num(c.reduction_pct, 2),
            "limiting": lim.key if lim else None,
            "limiting_label": lim.label if lim else None,
            "limiting_peak": _num(lim.peak, 4) if lim else None,
            "limiting_limit": _num(lim.limit, 4) if lim else None,
            # Post-repair readings: what is actually delivered.
            "readings": [r.as_dict() for r in c.report.readings],
            # Why the controller asked for this dose. Absent for cycle 1, whose dose is
            # a fixed 40.0 in the original scheduler rather than a fuzzy decision.
            "why": None if why is None else {
                "fes1_dose": _num(why["fes1_dose"], 3),
                "fes2_increase": _num(why["fes2_increase"], 4),
                "raw": _num(why["raw"], 3),
                "capped": why["capped"],
                "soft_cap": why["soft_cap"],
                "shrink": why["shrink"],
                "log10_cells": _num(math.log10(max(why["n_cells"], 1e-12)), 3),
                "toxicity": _num(why["toxicity"], 2),
                "fes1_inputs": why["fes1_inputs"],
                "fes2_inputs": why["fes2_inputs"],
            },
        })

    # Worst-case utilisation per organ across all delivered cycles, for the safety panel.
    organs = []
    for i, organ in enumerate(ORGANS):
        peaks = [c["readings"][i]["peak"] for c in cycles]
        limit = limits[organ.key]
        worst = max(peaks) if peaks else 0.0
        organs.append({
            "key": organ.key,
            "label": organ.label,
            "limit": _num(limit, 4),
            "peak": _num(worst, 4),
            "utilisation": _num(worst / limit, 4) if limit else None,
            "binding": any(c["limiting"] == organ.key for c in cycles),
            "per_cycle": [_num(p, 4) for p in peaks],
        })

    pbpk_cycles = []
    for c in result.cycles:
        pbpk_cycles.append({
            "index": c.index,
            "dose": _num(c.delivered_dose, 4),
            "states": {name: _series(c.profile.y[i], 4)
                       for i, name in enumerate(STATE_NAMES)},
        })

    return {
        "patient": {
            "name": cfg.patient_name,
            "weight_kg": cfg.weight_kg,
            "bsa_m2": _num(result.bsa_m2, 4),
            "bsa_scaled": _num(result.bsa_scaled, 4),
            "burden": burden.describe(cfg.n0),
            "burden_preset": (burden.nearest_preset(cfg.n0) or {}).get("label"),
        },
        "regimen": {
            "interval_days": cfg.interval_days,
            "horizon_days": cfg.horizon_days,
            "log_tumor": cfg.log_tumor,
            "legacy_dose_gate": cfg.legacy_dose_gate,
            "allometric": cfg.allometric,
            "volume_factor": _num(result.physiology.volume_factor, 4),
            "flow_factor": _num(result.physiology.flow_factor, 4),
        },
        "metrics": {k: (_num(v) if isinstance(v, float) else v)
                    for k, v in result.metrics.items()},
        "limits": {k: _num(v, 4) for k, v in limits.items()},
        "final": _tumor(result.final),
        "planned": _tumor(result.planned),
        "cycles": cycles,
        "organs": organs,
        "pbpk": {
            "t_hours": _series(result.cycles[0].profile.t_hours, 3) if result.cycles else [],
            "state_names": list(STATE_NAMES),
            "cycles": pbpk_cycles,
        },
        "reference": {"toxicity_limit": float(P.T_max)},
    }
