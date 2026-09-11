"""Initial tumour burden presets.

The original hardcoded N_0 = 1e10 for every patient, which is arguably the single most
important patient-specific input in the whole model. These presets expose it in terms a
clinician reads rather than a bare cell count.

The mass equivalents use the standard teaching convention **1 g of tumour ~ 10^9 cells**
(~1 cm^3). They are order-of-magnitude descriptors, not staging criteria - formal staging
is TNM and is not a function of cell number.
"""
import math

import numpy as np
import skfuzzy as fuzz

#: FES-1's tumour-size membership functions, for reporting which linguistic term a
#: given burden activates. Centres match fuzzy.py.
_FES1_TERMS = [
    ("VS", "Very small", -2.5),
    ("S", "Small", 0.0),
    ("M", "Medium", 3.0),
    ("B", "Big", 6.5),
    ("VB", "Very big", 10.0),
]
_UNIVERSE = np.arange(-3, 12, 1)

#: The model constrains the usable range at both ends:
#:  * FES-1/FES-2 read tumour size on a log10 universe of -3..11, and scikit-fuzzy gives
#:    ZERO membership outside it - so a burden above 1e11 makes the fuzzy controllers
#:    blind to tumour size, which is the defect this rewrite fixed.
#:  * The Gompertz carrying capacity rho_g is 1e12. At N = rho_g growth is exactly zero,
#:    and above it the tumour spontaneously shrinks. So 1e12 is an asymptote the model
#:    approaches, not a burden it can start from.
#: Hence 1e11 (~100 g) is the largest burden this tool will accept.
PRESETS = [
    {
        "key": "very_bulky", "cells": 1e11, "label": "Very bulky / widely metastatic",
        "detail": "~100 g across multiple sites. The largest burden the model supports.",
    },
    {
        "key": "bulky", "cells": 1e10, "label": "Bulky disease",
        "detail": "~10 g of tumour. This is the burden the published model assumes.",
    },
    {
        "key": "detectable", "cells": 1e9, "label": "Clinically detectable",
        "detail": "~1 g, about 1 cm³ - roughly the smallest lesion standard imaging resolves.",
    },
    {
        "key": "subclinical", "cells": 1e8, "label": "Subclinical",
        "detail": "~100 mg; below the imaging threshold but well above molecular detection.",
    },
    {
        "key": "mrd", "cells": 1e6, "label": "Minimal residual disease",
        "detail": "~1 mg; the range flow cytometry and PCR assays are built to detect.",
    },
    {
        "key": "molecular", "cells": 1e3, "label": "Molecular remission",
        "detail": "Trace disease, near the limit of any current assay.",
    },
]

PRESET_BY_KEY = {p["key"]: p for p in PRESETS}
DEFAULT_PRESET = "bulky"

MIN_CELLS = 1e2
MAX_CELLS = 1e11        # fuzzy universe ceiling; rho_g = 1e12 is the growth asymptote
FUZZY_LOG_MIN, FUZZY_LOG_MAX = -3.0, 11.0


def describe(cells):
    """Human-readable descriptors for a tumour cell count."""
    cells = float(cells)
    log10 = math.log10(cells) if cells > 0 else float("-inf")

    grams = cells / 1e9  # 1 g ~ 1e9 cells
    if grams >= 1000:
        mass = f"~{grams / 1000:.3g} kg"
    elif grams >= 1:
        mass = f"~{grams:.3g} g"
    elif grams >= 1e-3:
        mass = f"~{grams * 1e3:.3g} mg"
    else:
        mass = "trace"

    # Which FES-1 linguistic term this burden most activates.
    best, best_mu = None, 0.0
    for key, label, centre in _FES1_TERMS:
        mf = fuzz.gbellmf(_UNIVERSE, .8242, 3.278, centre)
        mu = float(fuzz.interp_membership(_UNIVERSE, mf, log10))
        if mu > best_mu:
            best, best_mu = label, mu

    return {
        "cells": cells,
        "log10": round(log10, 3),
        "mass": mass,
        "fuzzy_term": best,
        "fuzzy_membership": round(best_mu, 3),
        "in_fuzzy_range": FUZZY_LOG_MIN <= log10 <= FUZZY_LOG_MAX,
    }


def clamp(cells):
    """Hold a burden inside the range the model can actually represent."""
    return max(MIN_CELLS, min(MAX_CELLS, float(cells)))


def nearest_preset(cells):
    """The preset whose burden is closest in log space, or None if nothing is near."""
    if cells <= 0:
        return None
    target = math.log10(cells)
    best = min(PRESETS, key=lambda p: abs(math.log10(p["cells"]) - target))
    return best if abs(math.log10(best["cells"]) - target) < 0.05 else None
