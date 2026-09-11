"""Explore the regimen space for one patient.

Because a full simulation now costs ~0.2 s, the tool can evaluate a grid of cycle
intervals against a set of organ-limit presets and show the resulting trade-off between
tumour kill and peak toxicity. The frontier is reported; no regimen is ranked or
recommended - which one is acceptable is a clinical judgement, not a model output.
"""
from . import runner, safety

DEFAULT_INTERVALS = [7, 10, 14, 17, 21, 24, 28]
DEFAULT_PRESETS = ["legacy", "standard", "renal", "cardiac"]


def pareto_front(points):
    """Indices of regimens no other regimen beats on both kill and toxicity.

    A point is dominated when another achieves at least as much log-kill with no more
    peak toxicity, and is strictly better on one of the two.
    """
    front = []
    for i, p in enumerate(points):
        dominated = any(
            q["log_reduction"] >= p["log_reduction"]
            and q["peak_toxicity"] <= p["peak_toxicity"]
            and (q["log_reduction"] > p["log_reduction"] or q["peak_toxicity"] < p["peak_toxicity"])
            for j, q in enumerate(points) if j != i
        )
        if not dominated:
            front.append(i)
    return front


def explore(weight_kg, n0, intervals, preset_keys, horizon_days=120, progress=None):
    points = []
    total = max(len(intervals) * len(preset_keys), 1)
    done = 0

    for preset_key in preset_keys:
        preset = safety.PRESETS[preset_key]
        for interval in intervals:
            if progress:
                progress("sweep", 0.02 + 0.95 * done / total,
                         f"Simulating {interval}-day cycles, {preset['label']} "
                         f"({done + 1} of {total})…")
            result = runner.run(runner.RunConfig(
                patient_name="sweep", weight_kg=weight_kg, n0=n0,
                interval_days=interval, horizon_days=horizon_days,
                organ_limits=preset["limits"],
            ))
            m = result.metrics
            limiting = sorted({c.limiting_label for c in result.cycles if c.adjusted} - {None})
            points.append({
                "interval": interval,
                "preset": preset_key,
                "preset_label": preset["label"],
                "log_reduction": round(m["log_reduction"], 4),
                "peak_toxicity": round(m["peak_toxicity"], 3),
                "cycles_total": m["cycles_total"],
                "cycles_adjusted": m["cycles_adjusted"],
                "total_delivered": round(m["total_delivered"], 2),
                "limiting": ", ".join(limiting),
                "infeasible": m["any_infeasible"],
            })
            done += 1

    front = set(pareto_front(points))
    for i, p in enumerate(points):
        p["on_front"] = i in front

    if progress:
        progress("done", 1.0, "Complete")
    return {
        "points": points,
        "intervals": list(intervals),
        "presets": list(preset_keys),
        "toxicity_reference": 100.0,
    }
