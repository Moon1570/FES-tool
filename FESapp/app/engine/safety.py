"""Per-organ peak-concentration safety constraints.

Replaces the original ``feedback.py``. Two things change:

1.  Every organ's state indices are looked up from :data:`pbpk.STATE_NAMES` by name
    instead of being hard-coded. The original hard-coded a natural ordering that did
    not match the order ``dS2dt`` actually returns (heart is returned after tumour),
    so heart/kidney/muscle/fat/tumour were each compared against a neighbour's limit.
2.  The check returns a structured report instead of a bare ``True``/``False``, so the
    UI can show *which* organ limited the dose and by how much.
"""
import math
from dataclasses import dataclass, field
from typing import List, Optional

from .pbpk import STATE_INDEX


@dataclass(frozen=True)
class Organ:
    key: str
    label: str
    states: tuple          # names in STATE_NAMES that make up this organ

    @property
    def indices(self):
        return [STATE_INDEX[s] for s in self.states]


#: The 13 constrained organs, in the order the input form presents them.
ORGANS: List[Organ] = [
    Organ("ven_blood",  "Venous blood",   ("C_v", "C_rbcv")),
    Organ("lung",       "Lung",           ("C_lv", "C_le", "C_lb")),
    Organ("art_blood",  "Arterial blood", ("C_art", "C_rbca")),
    Organ("gut",        "Gut",            ("C_gv",)),
    Organ("brain",      "Brain",          ("C_bv", "C_be", "C_bb")),
    Organ("spleen",     "Spleen",         ("C_sv", "C_se", "C_sb")),
    Organ("liver",      "Liver",          ("C_liv", "C_lie", "C_lib")),
    Organ("heart",      "Heart",          ("C_hv", "C_he", "C_hb")),
    Organ("kidney",     "Kidney",         ("C_kv", "C_ke", "C_kb")),
    Organ("muscle",     "Muscle",         ("C_mv", "C_me", "C_mb")),
    Organ("fat",        "Fat",            ("C_fv", "C_fe", "C_fb")),
    Organ("tumor",      "Tumour",         ("C_tv", "C_te", "C_tb")),
    Organ("others",     "Other tissue",   ("C_ov", "C_oe", "C_ob")),
]

ORGAN_KEYS = [o.key for o in ORGANS]
ORGAN_BY_KEY = {o.key: o for o in ORGANS}

#: Matches the original hard-coded placeholder so behaviour is unchanged by default.
LEGACY_DEFAULT_LIMIT = 50.0


@dataclass
class OrganReading:
    key: str
    label: str
    peak: float           # max over the 48 h window of the summed sub-compartments
    limit: float
    breached: bool
    utilisation: float    # peak / limit; >1 means over the ceiling

    def as_dict(self):
        def num(v):
            return round(v, 5) if math.isfinite(v) else None
        return {
            "key": self.key, "label": self.label, "peak": num(self.peak),
            "limit": num(self.limit), "breached": self.breached,
            "utilisation": num(self.utilisation),
        }


@dataclass
class SafetyReport:
    readings: List[OrganReading] = field(default_factory=list)

    @property
    def breached(self) -> bool:
        return any(r.breached for r in self.readings)

    @property
    def limiting(self) -> Optional[OrganReading]:
        """The organ furthest over its ceiling (the one that forces a dose cut)."""
        over = [r for r in self.readings if r.breached]
        return max(over, key=lambda r: r.utilisation) if over else None

    def as_dict(self):
        lim = self.limiting
        return {
            "breached": self.breached,
            "limiting": lim.key if lim else None,
            "limiting_label": lim.label if lim else None,
            "readings": [r.as_dict() for r in self.readings],
        }


def organ_peak(y, organ: Organ) -> float:
    """Peak summed concentration for one organ across the whole simulated window."""
    return float(y[organ.indices].sum(axis=0).max())


def evaluate(result, limits) -> SafetyReport:
    """Build a :class:`SafetyReport` for one :class:`~.pbpk.PBPKResult`.

    ``limits`` maps organ key -> ceiling. Missing keys fall back to the legacy default.
    """
    readings = []
    for organ in ORGANS:
        limit = float(limits.get(organ.key, LEGACY_DEFAULT_LIMIT))
        peak = organ_peak(result.y, organ)
        readings.append(OrganReading(
            key=organ.key, label=organ.label, peak=peak, limit=limit,
            breached=peak > limit,
            utilisation=(peak / limit) if limit > 0 else float("inf"),
        ))
    return SafetyReport(readings=readings)


# ---------------------------------------------------------------------------
# Limit presets
# ---------------------------------------------------------------------------
# IMPORTANT: the underlying model's concentration units are not calibrated against
# clinical exposure limits, and the source publication defines no per-organ ceilings.
# Everything below 'legacy' is ILLUSTRATIVE - chosen so the safety feedback loop is
# actually exercised - and must be presented as such, never as clinical thresholds.
#
# Reference peaks (summed sub-compartments, 48 h window) across the default schedule
# (doses 14.3 .. 42.5), after the day-1 dose application was made uniform:
#   ven_blood  5.1-15.0   lung   5.1-15.1   art_blood 4.8-14.3   gut    2.6-7.8
#   brain      4.7-14.0   spleen 1.9-5.6    liver     0.8-2.5    heart  4.8-14.2
#   kidney     4.9-14.5   muscle 2.3-6.7    fat       3.5-10.4   tumour 4.1-12.3
#   other      1.7-4.9

LEGACY_LIMITS = {o.key: LEGACY_DEFAULT_LIMIT for o in ORGANS}

STANDARD_LIMITS = {
    "ven_blood": 13.0, "lung": 13.0, "art_blood": 12.5, "gut": 7.0,
    "brain": 12.5, "spleen": 5.0, "liver": 4.0, "heart": 12.5,
    "kidney": 13.0, "muscle": 6.0, "fat": 9.0, "tumor": 25.0, "others": 4.5,
}

REDUCED_RENAL_LIMITS = dict(STANDARD_LIMITS, kidney=11.0, liver=3.0)
REDUCED_CARDIAC_LIMITS = dict(STANDARD_LIMITS, heart=11.0)

PRESETS = {
    "legacy": {
        "label": "Published configuration",
        "description": "All ceilings at 50.0, as in the original tool. No dose is ever adjusted.",
        "limits": LEGACY_LIMITS,
    },
    "standard": {
        "label": "Standard adult",
        "description": "Illustrative ceilings that exercise the safety feedback loop.",
        "limits": STANDARD_LIMITS,
    },
    "renal": {
        "label": "Reduced renal tolerance",
        "description": "Lower kidney and liver ceilings, as for impaired clearance.",
        "limits": REDUCED_RENAL_LIMITS,
    },
    "cardiac": {
        "label": "Reduced cardiac tolerance",
        "description": "Lower heart ceiling, as for anthracycline cardiotoxicity risk.",
        "limits": REDUCED_CARDIAC_LIMITS,
    },
}

DEFAULT_PRESET = "standard"
