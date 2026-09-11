"""FES-1 and FES-2: the two Mamdani fuzzy controllers that set the dose.

Membership functions and rule bases are unchanged from ``fes1.py`` / ``fes2.py``.
Two things are different:

*   The ``ControlSystem`` (antecedents, membership functions, rule graph) is built
    **once at import** instead of being reconstructed on every call. The original
    rebuilt 21 + 30 rules and a networkx graph per invocation, from inside an ODE
    right-hand side; that dominated the runtime of the whole application.
*   ``tumour_size`` is fed ``log10(N)`` rather than the raw cell count. Its universe
    is ``arange(-3, 12)``, which is the range of log10(cells) for 1e-3 .. 1e11.
    With a raw count of ~1e10 the value falls outside the universe and scikit-fuzzy's
    ``interp_membership`` (``zero_outside_x=True``) returns **zero membership for every
    tumour-size term**, so the tumour-size dimension carried no information at all and
    the controllers responded only to toxicity. Pass ``log_tumor=False`` to reproduce
    the original behaviour.
"""
import threading

import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl

# scikit-fuzzy's ControlSystemSimulation is not safe to share between threads, and
# it memoises per input set; the rule graph underneath it is immutable and cheap to
# share. Build the graph once, hand out a fresh simulation per run.
_LOCK = threading.Lock()


def _build_fes1():
    tumor = ctrl.Antecedent(np.arange(-3, 12, 1), "tumor_size")
    tox = ctrl.Antecedent(np.arange(0, 121, 1), "toxicity")
    dose = ctrl.Consequent(np.arange(10, 51, 1), "dose")

    tumor["VS"] = fuzz.gbellmf(tumor.universe, .8242, 3.278, -2.5)
    tumor["S"] = fuzz.gbellmf(tumor.universe, .8242, 3.278, 0)
    tumor["M"] = fuzz.gbellmf(tumor.universe, .8242, 3.278, 3)
    tumor["B"] = fuzz.gbellmf(tumor.universe, .8242, 3.278, 6.5)
    tumor["VB"] = fuzz.gbellmf(tumor.universe, .8242, 3.278, 10)

    tox["VL"] = fuzz.gaussmf(tox.universe, 0.15, 13.89)
    tox["L"] = fuzz.trimf(tox.universe, [0, 30, 60])
    tox["M"] = fuzz.trimf(tox.universe, [30, 60, 90])
    tox["H"] = fuzz.trimf(tox.universe, [60, 90, 120])
    tox["VH"] = fuzz.gaussmf(tox.universe, 117.2, 2.421)

    dose["VVL"] = fuzz.trimf(dose.universe, [3.333, 10, 16.67])
    dose["VL"] = fuzz.trimf(dose.universe, [10, 16.67, 23.33])
    dose["L"] = fuzz.trimf(dose.universe, [16.67, 23.33, 30])
    dose["M"] = fuzz.trimf(dose.universe, [23.33, 30, 36.67])
    dose["H"] = fuzz.trimf(dose.universe, [30, 36.67, 43.33])
    dose["VH"] = fuzz.trimf(dose.universe, [36.67, 43.33, 50])
    dose["VVH"] = fuzz.trimf(dose.universe, [43.33, 50, 56.67])

    rules = [
        ctrl.Rule(tumor["VS"], dose["VVL"]),
        ctrl.Rule(tumor["S"] & tox["VL"], dose["M"]),
        ctrl.Rule(tumor["S"] & tox["L"], dose["L"]),
        ctrl.Rule(tumor["S"] & tox["M"], dose["VL"]),
        ctrl.Rule(tumor["S"] & tox["H"], dose["VVL"]),
        ctrl.Rule(tumor["S"] & tox["VH"], dose["VVL"]),
        ctrl.Rule(tumor["M"] & tox["VL"], dose["H"]),
        ctrl.Rule(tumor["M"] & tox["L"], dose["L"]),
        ctrl.Rule(tumor["M"] & tox["M"], dose["L"]),
        ctrl.Rule(tumor["M"] & tox["H"], dose["VL"]),
        ctrl.Rule(tumor["M"] & tox["VH"], dose["VVL"]),
        ctrl.Rule(tumor["B"] & tox["VL"], dose["VH"]),
        ctrl.Rule(tumor["B"] & tox["L"], dose["H"]),
        ctrl.Rule(tumor["B"] & tox["M"], dose["M"]),
        ctrl.Rule(tumor["B"] & tox["H"], dose["L"]),
        ctrl.Rule(tumor["B"] & tox["VH"], dose["VL"]),
        ctrl.Rule(tumor["VB"] & tox["VL"], dose["VH"]),
        ctrl.Rule(tumor["VB"] & tox["L"], dose["VH"]),
        ctrl.Rule(tumor["VB"] & tox["M"], dose["H"]),
        ctrl.Rule(tumor["VB"] & tox["H"], dose["M"]),
        ctrl.Rule(tumor["VB"] & tox["VH"], dose["L"]),
    ]
    return ctrl.ControlSystem(rules)


def _build_fes2():
    tumor = ctrl.Antecedent(np.arange(-3, 12, 1), "tumor_size")
    tox = ctrl.Antecedent(np.arange(0, 121, 1), "toxicity")
    calc = ctrl.Antecedent(np.arange(0, 301, 10), "calculated_dose")
    inc = ctrl.Consequent(np.arange(0, .9, .1), "percent_dose_increase")

    tumor["VS"] = fuzz.gaussmf(tumor.universe, -1.6, 2.2)
    tumor["S"] = fuzz.trimf(tumor.universe, [-.5, 2.25, 5])
    tumor["B"] = fuzz.trimf(tumor.universe, [3, 5.75, 8.5])
    tumor["VB"] = fuzz.gaussmf(tumor.universe, 9.85, 1.2)

    tox["VL"] = fuzz.trapmf(tox.universe, [-27, -3, 3, 30])
    tox["L"] = fuzz.trimf(tox.universe, [0, 30, 60])
    tox["M"] = fuzz.trimf(tox.universe, [30, 60, 90])
    tox["H"] = fuzz.trimf(tox.universe, [60, 90, 120])
    tox["VH"] = fuzz.trapmf(tox.universe, [90, 117, 123, 147])

    calc["LD"] = fuzz.trapmf(calc.universe, [0, 0, 165, 180])
    calc["ND"] = fuzz.trimf(calc.universe, [165, 180, 202])
    calc["OD"] = fuzz.trimf(calc.universe, [184, 202, 216])
    calc["VOD"] = fuzz.trapmf(calc.universe, [200, 245, 300, 300])

    inc["normal"] = fuzz.trapmf(inc.universe, [0, 0, 0, 0.1])
    inc["low_increase"] = fuzz.trimf(inc.universe, [0.05, 0.15, 0.25])
    inc["increase"] = fuzz.trimf(inc.universe, [0.2, 0.3, 0.4])
    # 'very_increase' is defined in the original but referenced by no rule, so the
    # controller can never output more than ~0.3. Kept for fidelity.
    inc["very_increase"] = fuzz.trapmf(inc.universe, [0.35, 0.45, 0.8, 0.8])

    rules = [
        ctrl.Rule(tumor["VS"] & calc["LD"], inc["normal"]),
        ctrl.Rule(tumor["VS"] & calc["ND"], inc["low_increase"]),
        ctrl.Rule(tumor["VS"] & calc["OD"], inc["increase"]),
        ctrl.Rule(tumor["VS"] & calc["VOD"], inc["increase"]),
        ctrl.Rule(tumor["S"] & tox["VL"] & calc["LD"], inc["normal"]),
        ctrl.Rule(tumor["S"] & tox["VL"] & calc["ND"], inc["low_increase"]),
        ctrl.Rule(tumor["S"] & tox["VL"] & calc["OD"], inc["increase"]),
        ctrl.Rule(tumor["S"] & tox["L"] & calc["LD"], inc["normal"]),
        ctrl.Rule(tumor["S"] & tox["L"] & calc["ND"], inc["low_increase"]),
        ctrl.Rule(tumor["S"] & tox["L"] & calc["OD"], inc["increase"]),
        ctrl.Rule(tumor["S"] & tox["M"] & calc["LD"], inc["normal"]),
        ctrl.Rule(tumor["S"] & tox["M"] & calc["ND"], inc["normal"]),
        ctrl.Rule(tumor["S"] & tox["M"] & calc["OD"], inc["increase"]),
        ctrl.Rule(tumor["S"] & tox["H"], inc["normal"]),
        ctrl.Rule(tumor["S"] & tox["VH"], inc["normal"]),
        ctrl.Rule(tumor["S"] & tox["VL"] & calc["VOD"], inc["increase"]),
        ctrl.Rule(tumor["S"] & tox["L"] & calc["VOD"], inc["increase"]),
        ctrl.Rule(tumor["S"] & tox["M"] & calc["VOD"], inc["increase"]),
        ctrl.Rule(tumor["B"] & tox["VL"], inc["normal"]),
        ctrl.Rule(tumor["B"] & tox["L"], inc["normal"]),
        ctrl.Rule(tumor["B"] & tox["M"], inc["normal"]),
        ctrl.Rule(tumor["B"] & tox["H"], inc["normal"]),
        ctrl.Rule(tumor["B"] & tox["VH"], inc["normal"]),
        ctrl.Rule(tumor["VB"] & tox["VL"], inc["normal"]),
        ctrl.Rule(tumor["VB"] & tox["L"], inc["normal"]),
        ctrl.Rule(tumor["VB"] & tox["M"] & calc["LD"], inc["normal"]),
        ctrl.Rule(tumor["VB"] & tox["M"] & calc["ND"], inc["normal"]),
        ctrl.Rule(tumor["VB"] & tox["M"] & calc["VOD"], inc["increase"]),
        ctrl.Rule(tumor["VB"] & tox["H"], inc["normal"]),
        ctrl.Rule(tumor["VB"] & tox["VH"], inc["normal"]),
    ]
    return ctrl.ControlSystem(rules)


FES1_SYSTEM = _build_fes1()
FES2_SYSTEM = _build_fes2()


#: Linguistic terms, for reporting which ones a given input activates.
FES1_TUMOR_TERMS = {"VS": "Very small", "S": "Small", "M": "Medium", "B": "Big", "VB": "Very big"}
TOX_TERMS = {"VL": "Very low", "L": "Low", "M": "Medium", "H": "High", "VH": "Very high"}
DOSE_TERMS = {"VVL": "Very very low", "VL": "Very low", "L": "Low", "M": "Medium",
              "H": "High", "VH": "Very high", "VVH": "Very very high"}
FES2_TUMOR_TERMS = {"VS": "Very small", "S": "Small", "B": "Big", "VB": "Very big"}
BSA_TERMS = {"LD": "Low dose", "ND": "Normal dose", "OD": "Over dose", "VOD": "Very over dose"}


def _activations(sim, variable, labels):
    """Membership of each term of ``variable`` for the value just computed.

    scikit-fuzzy stores the fuzzified membership on each Term keyed by the simulation,
    so this reads the actual degrees that drove the last inference rather than
    recomputing them.
    """
    out = []
    for key, term in variable.terms.items():
        try:
            mu = float(term.membership_value[sim])
        except (KeyError, TypeError):
            continue
        if mu > 0.001:
            out.append({"term": key, "label": labels.get(key, key), "mu": round(mu, 3)})
    return sorted(out, key=lambda d: -d["mu"])


def tumor_input(n_cells, log_tumor=True):
    """Map a tumour cell count onto the controllers' ``tumor_size`` universe."""
    if not log_tumor:
        return n_cells
    return float(np.log10(max(float(n_cells), 1e-3)))


class FuzzyController:
    """Per-run wrapper holding fresh simulation objects over the shared rule graphs."""

    def __init__(self, log_tumor=True):
        self.log_tumor = log_tumor
        with _LOCK:
            self._fes1 = ctrl.ControlSystemSimulation(FES1_SYSTEM)
            self._fes2 = ctrl.ControlSystemSimulation(FES2_SYSTEM)

    def fes1(self, n_cells, toxicity):
        """Base dose from tumour burden and cumulative toxicity."""
        self._fes1.input["tumor_size"] = tumor_input(n_cells, self.log_tumor)
        self._fes1.input["toxicity"] = toxicity
        self._fes1.compute()
        return float(self._fes1.output["dose"])

    def fes1_explained(self, n_cells, toxicity):
        """Same as :meth:`fes1`, plus the term memberships that produced it."""
        dose = self.fes1(n_cells, toxicity)
        ants = FES1_SYSTEM.antecedents
        by_label = {a.label: a for a in ants}
        return dose, {
            "tumor": _activations(self._fes1, by_label["tumor_size"], FES1_TUMOR_TERMS),
            "toxicity": _activations(self._fes1, by_label["toxicity"], TOX_TERMS),
        }

    def fes2(self, n_cells, toxicity, calculated_dose):
        """Fractional dose increase from tumour burden, toxicity and scaled BSA."""
        self._fes2.input["tumor_size"] = tumor_input(n_cells, self.log_tumor)
        self._fes2.input["toxicity"] = toxicity
        self._fes2.input["calculated_dose"] = calculated_dose
        self._fes2.compute()
        return float(self._fes2.output["percent_dose_increase"])

    def fes2_explained(self, n_cells, toxicity, calculated_dose):
        """Same as :meth:`fes2`, plus the term memberships that produced it."""
        pct = self.fes2(n_cells, toxicity, calculated_dose)
        by_label = {a.label: a for a in FES2_SYSTEM.antecedents}
        return pct, {
            "tumor": _activations(self._fes2, by_label["tumor_size"], FES2_TUMOR_TERMS),
            "toxicity": _activations(self._fes2, by_label["toxicity"], TOX_TERMS),
            "bsa": _activations(self._fes2, by_label["calculated_dose"], BSA_TERMS),
        }
