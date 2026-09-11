"""35-state physiologically-based pharmacokinetic (PBPK) model.

The state vector order defined by :data:`STATE_NAMES` is the single source of truth
for the whole package. In the original code the ordering lived only in the tuple
unpacking inside ``dS2dt``; heart was moved out of numeric order (it carries the
``#EQ 18/19/20`` comments but is returned last-but-one), and ``feedback.py``
independently hard-coded a *different*, natural ordering. The result was that five
organs were checked against the wrong safety limit. Deriving every index from
:data:`STATE_NAMES` makes that class of bug unrepresentable.
"""
from dataclasses import dataclass

import numpy as np
from . import params as P
from ._solver import solve as solve_ivp

#: Names of the 35 PBPK states, in the exact order the RHS returns them.
#: Suffixes: ``v`` vascular, ``e`` extravascular, ``b`` bound.
STATE_NAMES = [
    "C_v", "C_rbcv",                  # venous blood plasma, venous RBC
    "C_lv", "C_le", "C_lb",           # lung
    "C_art", "C_rbca",                # arterial blood plasma, arterial RBC
    "C_gv",                           # gut (vascular only)
    "C_bv", "C_be", "C_bb",           # brain
    "C_sv", "C_se", "C_sb",           # spleen
    "C_liv", "C_lie", "C_lib",        # liver
    "C_kv", "C_ke", "C_kb",           # kidney
    "C_mv", "C_me", "C_mb",           # muscle
    "C_fv", "C_fe", "C_fb",           # fat
    "C_tv", "C_te", "C_tb",           # tumour
    "C_hv", "C_he", "C_hb",           # heart  (returned here, not in EQ order)
    "C_ov", "C_oe", "C_ob",           # other
]

#: name -> index, for callers that want a state by name rather than position.
STATE_INDEX = {name: i for i, name in enumerate(STATE_NAMES)}

# Integration window: 48 hours expressed in days, sampled hourly.
T_RANGE = (0, 2)
N_SAMPLES = 49
ATOL = RTOL = 1.49e-8

#: The published parameter set describes a 70 kg reference human. Without scaling, a
#: 50 kg and a 120 kg patient get identical predicted organ exposures - body weight only
#: reached FES-2 through BSA. Standard allometry: compartment volumes scale linearly with
#: body weight, blood flows with BW^0.75. First-order transfer rate constants are left
#: alone, which is what makes elimination slow with BW^-0.25 as observed.
REFERENCE_WEIGHT_KG = 70.0
FLOW_EXPONENT = 0.75

_VOLUME_NAMES = [
    "V_ven", "V_art", "V_lv", "V_le", "V_gv", "V_bv", "V_be", "V_sv", "V_se",
    "V_liv", "V_lie", "V_kv", "V_ke", "V_mv", "V_me", "V_fv", "V_fe", "V_tv",
    "V_te", "V_hv", "V_he", "V_ov", "V_oe",
]
_FLOW_NAMES = [
    "F_li", "F_l", "F_t", "F_g", "F_m", "F_s", "F_h", "F_f", "F_k", "F_b",
    "F_o", "F_tot",
]


class Physiology:
    """Body-weight-scaled copy of the PBPK parameter set.

    At the reference weight every factor is exactly 1.0, so a 70 kg patient reproduces
    the published numbers bit for bit.
    """

    __slots__ = _VOLUME_NAMES + _FLOW_NAMES + ["weight_kg", "volume_factor", "flow_factor"]

    def __init__(self, weight_kg=REFERENCE_WEIGHT_KG, enabled=True):
        ratio = (float(weight_kg) / REFERENCE_WEIGHT_KG) if enabled else 1.0
        self.weight_kg = float(weight_kg)
        self.volume_factor = ratio
        self.flow_factor = ratio ** FLOW_EXPONENT
        for name in _VOLUME_NAMES:
            setattr(self, name, getattr(P, name) * self.volume_factor)
        for name in _FLOW_NAMES:
            setattr(self, name, getattr(P, name) * self.flow_factor)


REFERENCE_PHYSIOLOGY = Physiology()


@dataclass
class PBPKResult:
    """Concentration-time profiles for one administered dose."""

    dose: float
    t_hours: np.ndarray
    y: np.ndarray  # shape (35, N_SAMPLES), rows ordered as STATE_NAMES

    def state(self, name):
        return self.y[STATE_INDEX[name]]


def _rhs(t, S, dose_rate, legacy_dose_gate=False, phys=REFERENCE_PHYSIOLOGY):
    (C_v, C_rbcv, C_lv, C_le, C_lb, C_art, C_rbca, C_gv, C_bv, C_be, C_bb,
     C_sv, C_se, C_sb, C_liv, C_lie, C_lib, C_kv, C_ke, C_kb, C_mv, C_me,
     C_mb, C_fv, C_fe, C_fb, C_tv, C_te, C_tb, C_hv, C_he, C_hb, C_ov,
     C_oe, C_ob) = S

    # The dose enters as a constant rate and is applied again on the second day.
    # The original gated this on `int(dose) != 40`, which - because int() truncates -
    # silently exempted every dose in [40, 41). That produced a discontinuity where
    # 39.99 gave roughly twice the exposure of 40.00. The gate was a leftover from
    # when the day-0 dose was hardcoded to 40; it is applied uniformly now.
    dose = dose_rate
    if int(t) == 1 and not (legacy_dose_gate and int(dose_rate) == 40):
        dose = dose + dose_rate

    return [
        ((((phys.F_l*C_le)+(phys.F_b*C_be)+(phys.F_s*C_se)+(phys.F_li*C_lie)+(phys.F_h*C_he)+(phys.F_k*C_ke)+(phys.F_m*C_me)+(phys.F_f*C_fe)+(phys.F_t*C_te)+(phys.F_o*C_oe)) - (phys.F_tot * C_v))/(phys.V_ven * P.one_sub_f_hem))+(dose/(phys.V_ven*P.one_sub_f_hem))+(P.one_by_one_sub_f_hem*P.f_hem*P.k_rbcplas*C_rbcv)-(P.k_plasrbc*P.f_unb*C_v),  # EQ 1
        (P.one_sub_f_hem*P.k_plasrbc*P.f_unb*C_v/P.f_hem)-(P.k_rbcplas*C_rbcv),  # EQ 2
        ((phys.F_l/phys.V_lv)*(C_v-C_lv)) - (P.k_lve*P.f_unb*C_lv) + ((phys.V_le*P.k_lev*C_le)/phys.V_lv),  # EQ 3
        (phys.V_lv*P.k_lve*P.f_unb*C_lv/phys.V_le)-(P.k_lev*C_le)+(P.k_bind_out*C_lb)-(P.k_bind_in*C_le),  # EQ 4
        (P.k_bind_in*C_le)-(P.k_bind_out*C_lb),  # EQ 5
        (((phys.F_l*C_lv)-(phys.F_tot*C_art))/(phys.V_art*P.one_sub_f_hem))+(P.f_hem*P.k_rbcplas*C_rbca/P.one_sub_f_hem)-(P.k_plasrbc*P.f_unb*C_art),  # EQ 6
        (P.one_sub_f_hem*P.k_rbcplas*P.f_unb*C_art/P.f_hem) - (P.k_rbcplas*C_rbca),  # EQ 7
        ((phys.F_g/phys.V_gv)*(C_art - C_gv)),  # EQ 8
        ((phys.F_b/phys.V_bv)*(C_art - C_bv)) - (P.k_bve*P.f_unb*C_bv) + ((phys.V_be*P.k_bev*C_be)/phys.V_bv),  # EQ 9
        (phys.V_bv*P.k_bve*P.f_unb*C_bv/phys.V_be)-(P.k_bev*C_be)+(P.k_bind_out*C_bb)-(P.k_bind_in*C_be),  # EQ 10
        (P.k_bind_in*C_be)-(P.k_bind_out*C_bb),  # EQ 11
        ((phys.F_s/phys.V_sv)*(C_art - C_sv)) - (P.k_sve*P.f_unb*C_sv) + ((phys.V_se*P.k_sev*C_se)/phys.V_sv),  # EQ 12
        (phys.V_sv*P.k_sve*P.f_unb*C_sv/phys.V_se)-(P.k_sev*C_se)+(P.k_bind_out*C_sb)-(P.k_bind_in*C_se),  # EQ 13
        (P.k_bind_in*C_se)-(P.k_bind_out*C_sb),  # EQ 14
        ((phys.F_li*C_art)+(phys.F_g*C_gv)+(phys.F_s*C_sv)-(C_liv*(phys.F_g+phys.F_s+phys.F_li))/phys.V_liv)-(P.k_live*P.f_unb*C_liv)+(phys.V_lie*P.k_liev*C_lie/phys.V_liv),  # EQ 15
        (phys.V_liv*P.k_live*P.f_unb*C_liv/phys.V_lie)-(P.k_liev*C_lie)+(P.k_bind_out*C_lib)-(P.k_bind_in*C_lie)-(P.k_clli*C_lie),  # EQ 16
        (P.k_bind_in*C_lie)-(P.k_bind_out*C_lib)-(P.k_clli*C_lib),  # EQ 17
        ((phys.F_k/phys.V_kv)*(C_art - C_kv)) - (P.k_kve*P.f_unb*C_kv) + ((phys.V_ke*P.k_kev*C_ke)/phys.V_kv),  # EQ 21
        (phys.V_kv*P.k_kve*P.f_unb*C_kv/phys.V_ke)-(P.k_kev*C_ke)+(P.k_bind_out*C_kb)-(P.k_bind_in*C_ke),  # EQ 22
        (P.k_bind_in*C_ke)-(P.k_bind_out*C_kb),  # EQ 23
        ((phys.F_m/phys.V_mv)*(C_art - C_mv)) - (P.k_mve*P.f_unb*C_mv) + ((phys.V_me*P.k_mev*C_me)/phys.V_mv),  # EQ 24
        (phys.V_mv*P.k_mve*P.f_unb*C_mv/phys.V_me)-(P.k_mev*C_me)+(P.k_bind_out*C_mb)-(P.k_bind_in*C_me),  # EQ 25
        (P.k_bind_in*C_me)-(P.k_bind_out*C_mb),  # EQ 26
        (phys.F_f*(C_art-C_fv)/phys.V_fv)-(P.k_fve*P.f_unb*C_fv)+(phys.V_fe*P.k_fev*C_fe/phys.V_fv),  # EQ 27
        (phys.V_fv*P.k_fve*P.f_unb*C_fv/phys.V_fe)-(P.k_fev*C_fe)+(P.k_bind_out*C_fb)-(P.k_bind_in*C_fe),  # EQ 28
        (P.k_bind_in*C_fe)-(P.k_bind_out*C_fb),  # EQ 29
        (phys.F_t*(C_art-C_tv)/phys.V_tv)-(P.k_tve*P.f_unb*C_tv)+(phys.V_te*P.k_tev*C_te/phys.V_tv),  # EQ 30
        (phys.V_tv*P.k_tve*P.f_unb*C_tv/phys.V_te)-(P.k_tev*C_te)+(P.k_bind_out*C_tb)-(P.k_bind_in*C_te),  # EQ 31
        (P.k_bind_in*C_te)-(P.k_bind_out*C_tb),  # EQ 32
        (phys.F_h*(C_art-C_hv)/phys.V_hv)-(P.k_hve*P.f_unb*C_hv)+(phys.V_he*P.k_hev*C_he/phys.V_hv),  # EQ 18
        (phys.V_hv*P.k_hve*P.f_unb*C_hv/phys.V_he)-(P.k_hev*C_he)+(P.k_bind_out*C_hb)-(P.k_bind_in*C_he),  # EQ 19
        (P.k_bind_in*C_he)-(P.k_bind_out*C_hb),  # EQ 20
        (phys.F_o*(C_art-C_ov)/phys.V_ov)-(P.k_ove*P.f_unb*C_ov)+(phys.V_oe*P.k_oev*C_oe/phys.V_ov),  # EQ 33
        (phys.V_ov*P.k_ove*P.f_unb*C_ov/phys.V_oe)-(P.k_oev*C_oe)+(P.k_bind_out*C_ob)-(P.k_bind_in*C_oe),  # EQ 34
        (P.k_bind_in*C_oe)-(P.k_bind_out*C_ob),  # EQ 35
    ]


def simulate_pbpk(dose, legacy_dose_gate=False, phys=REFERENCE_PHYSIOLOGY):
    """Integrate the 35-compartment model for a single administered ``dose``.

    ``legacy_dose_gate=True`` restores the original ``int(dose) != 40`` exemption,
    used only to prove the port reproduces the published numbers.
    """
    t_eval = np.linspace(T_RANGE[0], T_RANGE[1], N_SAMPLES, dtype=float)
    sol = solve_ivp(
        _rhs, T_RANGE, y0=[0.0] * len(STATE_NAMES), method="LSODA",
        t_eval=t_eval, atol=ATOL, rtol=RTOL, args=(dose, legacy_dose_gate, phys),
    )
    return PBPKResult(dose=float(dose), t_hours=sol.t * 24.0, y=sol.y)
