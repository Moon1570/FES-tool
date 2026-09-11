"""Orchestrates a full simulation run.

Pipeline (unchanged from the original):

1.  Simulate the tumour model, letting FES-1/FES-2 choose a dose each cycle.
2.  Take the dose administered at the start of each cycle.
3.  Push each dose through the 35-state PBPK model.
4.  While any organ exceeds its ceiling, shrink that dose by 5% and re-solve.
5.  Re-simulate the tumour with the corrected schedule.

Everything is run-local: no module state is mutated, so concurrent runs are safe and
a second run can never inherit the first one's schedule.
"""
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

from . import params as P
from . import pbpk, safety
from .fuzzy import FuzzyController
from .tumor import SHRINK, TumorConfig, simulate_adaptive, simulate_replay

#: Guards on the repair loop. The original had neither, so an unsatisfiable limit
#: looped forever, each iteration costing a 35-ODE solve.
MAX_REPAIR_ITERATIONS = 80
MIN_DOSE = 0.5


def calc_bsa(weight_kg):
    """Costeff BSA, scaled by 100.

    ``(4W + 7) / (W + 90)`` gives BSA in m². The original multiplies by 100 so the
    value lands on FES-2's ``calculated_dose`` universe, whose 'normal dose' term is
    centred at 180. Both forms are returned so the UI can show m² and the controller
    can keep its scaled input.
    """
    scaled = (((weight_kg * 4) + 7) / (90 + weight_kg)) * 100
    return scaled / 100.0, scaled


@dataclass
class RunConfig:
    patient_name: str = "John Doe"
    weight_kg: float = 70.0
    interval_days: int = 14
    horizon_days: int = 120
    organ_limits: Dict[str, float] = field(default_factory=dict)
    n0: float = float(P.N_0)         # initial tumour burden, in cells
    log_tumor: bool = True
    allometric: bool = True          # scale PBPK physiology with body weight
    legacy_dose_gate: bool = False   # only for reproducing the original numbers

    def limits(self):
        return {o.key: float(self.organ_limits.get(o.key, safety.LEGACY_DEFAULT_LIMIT))
                for o in safety.ORGANS}


@dataclass
class CycleOutcome:
    index: int                     # 1-based cycle number
    day: int
    planned_dose: float            # what the fuzzy controllers asked for
    delivered_dose: float          # what survived the safety check
    adjusted: bool
    iterations: int
    infeasible: bool
    report: safety.SafetyReport          # after repair: what is actually delivered
    initial_report: safety.SafetyReport  # before repair: what forced the reduction
    profile: pbpk.PBPKResult

    @property
    def reduction_pct(self):
        if self.planned_dose <= 0:
            return 0.0
        return 100.0 * (1 - self.delivered_dose / self.planned_dose)

    @property
    def limiting(self):
        """The organ that forced this cycle's reduction (from the pre-repair check)."""
        return self.initial_report.limiting

    @property
    def limiting_label(self):
        lim = self.limiting
        return lim.label if lim else None


@dataclass
class RunResult:
    config: RunConfig
    bsa_m2: float
    bsa_scaled: float
    physiology: object
    planned: object                 # TumorResult from the adaptive pass
    final: object                   # TumorResult from the corrected pass
    cycles: List[CycleOutcome]
    decisions: Dict[int, dict]

    @property
    def metrics(self):
        n = self.final.N
        log0, logend = self.final.log10N[0], self.final.log10N[-1]
        adjusted = [c for c in self.cycles if c.adjusted]
        return {
            "n0": float(self.config.n0),
            "log10_n_start": float(log0),
            "log10_n_end": float(logend),
            "log_reduction": float(log0 - logend),
            "n_end": float(n[-1]),
            "peak_toxicity": float(self.final.T.max()),
            "toxicity_limit": float(P.T_max),
            "cycles_total": len(self.cycles),
            "cycles_adjusted": len(adjusted),
            "total_delivered": float(self.final.dose_by_day.sum()),
            "total_planned": float(self.planned.dose_by_day.sum()),
            "any_infeasible": any(c.infeasible for c in self.cycles),
        }


def _noop(stage, frac, message):
    pass


def run(config: RunConfig, progress: Optional[Callable] = None) -> RunResult:
    progress = progress or _noop
    bsa_m2, bsa_scaled = calc_bsa(config.weight_kg)
    limits = config.limits()

    tcfg = TumorConfig(
        interval_days=config.interval_days, horizon_days=config.horizon_days,
        n0=config.n0, bsa_scaled=bsa_scaled, log_tumor=config.log_tumor,
    )
    phys = pbpk.Physiology(config.weight_kg, enabled=config.allometric)

    progress("tumour", 0.05, "Simulating tumour response…")
    controller = FuzzyController(log_tumor=config.log_tumor)
    planned, schedule = simulate_adaptive(controller, tcfg)

    planned_cycles = list(planned.cycle_doses)
    total = max(len(planned_cycles), 1)
    cycles: List[CycleOutcome] = []

    for i, planned_dose in enumerate(planned_cycles):
        day = i * config.interval_days
        progress("pbpk", 0.10 + 0.75 * (i / total),
                 f"Checking organ safety for cycle {i + 1} of {total}…")

        amount = planned_dose
        profile = pbpk.simulate_pbpk(amount, config.legacy_dose_gate, phys)
        report = safety.evaluate(profile, limits)
        initial_report = report   # the breach that justifies any reduction
        iterations = 0
        infeasible = False

        while report.breached:
            if iterations >= MAX_REPAIR_ITERATIONS or amount * SHRINK < MIN_DOSE:
                infeasible = True
                break
            amount *= SHRINK
            iterations += 1
            profile = pbpk.simulate_pbpk(amount, config.legacy_dose_gate, phys)
            report = safety.evaluate(profile, limits)

        cycles.append(CycleOutcome(
            index=i + 1, day=day, planned_dose=planned_dose, delivered_dose=amount,
            adjusted=iterations > 0, iterations=iterations, infeasible=infeasible,
            report=report, initial_report=initial_report, profile=profile,
        ))

    progress("replay", 0.90, "Re-simulating with the corrected schedule…")
    final = simulate_replay([c.delivered_dose for c in cycles], tcfg)

    progress("done", 1.0, "Complete")
    return RunResult(
        config=config, bsa_m2=bsa_m2, bsa_scaled=bsa_scaled, physiology=phys,
        planned=planned, final=final, cycles=cycles, decisions=schedule.decisions,
    )
