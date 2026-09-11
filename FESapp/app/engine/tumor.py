"""Martin's 3-state tumour / pharmacokinetic / toxicity model.

States: ``C`` plasma drug concentration, ``N`` tumour cell count, ``T`` cumulative
toxicity, integrated over the treatment horizon in days.

The original had two near-identical right-hand sides (``dSdt`` and ``dS3dt``) that
differed only in where the dose came from, and both mutated a module-level ``dose``
list that was never reset between requests. Here the schedule is an object owned by
a single run, so two simulations can never see each other's doses.
"""
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from numpy import log as ln
from . import params as P
from ._solver import solve as solve_ivp

# Magic constants from the original scheduler, named so they can be seen and tuned.
FIRST_DOSE = 40.0          # day-0 dose, independent of the patient
SOFT_CAP = 43.0            # above this, the day-k dose is shrunk once by SHRINK
SOFT_CAP_NEXT_DAY = 43.7   # the day-(k+1) equivalent; inconsistent with SOFT_CAP
SHRINK = 0.95

#: The Gompertz growth term contains ln(rho / N), which is undefined once the tumour
#: is driven to zero. Short cycle intervals do exactly that, and the resulting NaN
#: propagated through the rest of the solve. A cell count is physically non-negative,
#: so the growth term is evaluated at a tiny positive floor instead. This is far below
#: one cell and cannot affect any trajectory that stays in a meaningful range.
N_FLOOR = 1e-9


@dataclass
class TumorConfig:
    interval_days: int = 14
    horizon_days: int = 120
    n0: float = float(P.N_0)
    bsa_scaled: float = 179.375   # BSA * 100, the FES-2 'calculated_dose' input
    log_tumor: bool = True


@dataclass
class TumorResult:
    t: np.ndarray
    C: np.ndarray
    N: np.ndarray
    T: np.ndarray
    dose_by_day: np.ndarray       # dose administered on each whole day
    cycle_doses: List[float] = field(default_factory=list)   # one per dose cycle

    @property
    def log10N(self):
        return np.log10(np.clip(self.N, 1e-12, None))


class _Schedule:
    """Run-local dose bookkeeping shared by the adaptive and replay right-hand sides."""

    def __init__(self, horizon_days):
        self.dose = [0.0] * (horizon_days + 1)
        self.horizon = horizon_days

    def day(self, t):
        return min(int(t), self.horizon)


class AdaptiveSchedule(_Schedule):
    """Doses decided on the fly by FES-1 and FES-2 (replaces ``dSdt``)."""

    def __init__(self, controller, cfg: TumorConfig):
        super().__init__(cfg.horizon_days)
        self.controller = controller
        self.cfg = cfg
        self.decisions = {}   # day -> dict, for the treatment calendar / audit trail

    def dose_at(self, t, N_t, T_t):
        d = self.day(t)
        interval = self.cfg.interval_days

        if d == 0:
            self.dose[0] = FIRST_DOSE
            return FIRST_DOSE
        if d == 1:
            self.dose[1] = 0.0
            return 0.0

        if d % interval == 0:
            # Memoised for this run only: LSODA evaluates the RHS many times per day,
            # and the fuzzy decision must be identical across those evaluations.
            if self.dose[d] != 0:
                return self.dose[d]
            base, why1 = self.controller.fes1_explained(N_t, T_t)
            pct, why2 = self.controller.fes2_explained(N_t, T_t, self.cfg.bsa_scaled)
            value = base + base * pct
            capped = value > SOFT_CAP
            if capped:
                value = value * SHRINK
            self.dose[d] = value
            self.decisions[d] = {
                "day": d, "fes1_dose": base, "fes2_increase": pct,
                "raw": base + base * pct, "capped": capped,
                "soft_cap": SOFT_CAP, "shrink": SHRINK, "dose": value,
                "n_cells": float(N_t), "toxicity": float(T_t),
                "fes1_inputs": why1, "fes2_inputs": why2,
            }
            return value

        if d % interval == 1:
            # The original recomputed FES-1/FES-2 here and threw the result away.
            value = self.dose[d - 1]
            if value > SOFT_CAP_NEXT_DAY:
                value = value * SHRINK
            self.dose[d] = value
            return value

        self.dose[d] = 0.0
        return 0.0


class ReplaySchedule(_Schedule):
    """Doses read from a fixed per-cycle list (replaces ``dS3dt``)."""

    def __init__(self, cycle_doses, cfg: TumorConfig):
        super().__init__(cfg.horizon_days)
        self.cycle_doses = list(cycle_doses)
        self.cfg = cfg

    def _cycle(self, d):
        return min(d // self.cfg.interval_days, len(self.cycle_doses) - 1)

    def dose_at(self, t, N_t, T_t):
        d = self.day(t)
        interval = self.cfg.interval_days

        if d == 0:
            self.dose[0] = self.cycle_doses[0]
            return self.cycle_doses[0]
        if d == 1:
            return 0.0
        if d % interval in (0, 1):
            value = self.cycle_doses[self._cycle(d)]
            self.dose[d] = value
            return value

        self.dose[d] = 0.0
        return 0.0


def _rhs(t, S, schedule):
    C_t, N_t, T_t = S
    D_t = schedule.dose_at(t, N_t, T_t)

    C_eff = (C_t - P.C_th) if C_t >= P.C_th else 0.0
    N_safe = N_t if N_t > N_FLOOR else N_FLOOR

    return [
        D_t - P.lambdaa * C_t,
        ((1 / P.tau_g) * ln((ln(P.rho_g / P.N_0)) / ln(P.rho_g / (2 * P.N_0)))
         * N_safe * ln(P.rho_g / N_safe)) - (P.K_eff * C_eff * N_t),
        C_t - (P.eta * T_t),
    ]


def _integrate(schedule, cfg: TumorConfig) -> TumorResult:
    horizon = cfg.horizon_days
    t_eval = np.linspace(0, horizon, horizon + 1, dtype=int)
    sol = solve_ivp(_rhs, (0, horizon), y0=[P.C_0, cfg.n0, 0],
                    method="LSODA", t_eval=t_eval, args=(schedule,))
    dose_by_day = np.array(schedule.dose, dtype=float)
    cycles = [dose_by_day[d] for d in range(0, horizon + 1) if d % cfg.interval_days == 0]
    return TumorResult(t=sol.t, C=sol.y[0], N=sol.y[1], T=sol.y[2],
                       dose_by_day=dose_by_day, cycle_doses=cycles)


def simulate_adaptive(controller, cfg: TumorConfig):
    """First pass: let the fuzzy controllers choose the schedule."""
    schedule = AdaptiveSchedule(controller, cfg)
    return _integrate(schedule, cfg), schedule


def simulate_replay(cycle_doses, cfg: TumorConfig) -> TumorResult:
    """Second pass: re-simulate with the safety-corrected schedule."""
    return _integrate(ReplaySchedule(cycle_doses, cfg), cfg)
