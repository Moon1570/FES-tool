"""Model constants for the FES chemotherapy scheduling tool.

Tumour/PK constants (Martin's model) and the ~90 physiological PBPK constants are
reproduced verbatim from the original ``views.py`` so numerical results are unchanged.
"""

# ---------------------------------------------------------------------------
# Martin's tumour / pharmacokinetic / toxicity model
# ---------------------------------------------------------------------------
C_0 = 0
N_0 = 10000000000        # initial tumour cell count
tau_g = 150              # Gompertz growth time constant (days)
rho_g = 1000000000000    # Gompertz carrying capacity
K_eff = 2.7 * 10 ** -2   # cell-kill coefficient
lambdaa = .27            # drug elimination rate (1/day)
eta = .4                 # toxicity clearance rate (1/day)
C_th = 10                # concentration threshold below which there is no kill

# Declared in the original model but never enforced. Kept as documented, named
# reference values rather than dead code; see safety.py for the constraints that
# ARE enforced (per-organ peak concentration limits).
D_th = 10
D_max = 50
T_max = 100

# ---------------------------------------------------------------------------
# PBPK model
# ---------------------------------------------------------------------------
k_live = 10.251
k_clli = 0.1023
k_lev = 0.0365
k_tev = 0.0006
k_mev = 0.0158
k_sev = 0.0445
k_hev = 0.0495
k_fev = 0.0079
k_kev = 0.1859
k_bev = 0.0573
k_oev = .0099
k_liev = 0.0965
k_lve = 0.2662
k_tve = 0.110
k_mve = 0.5952
k_sve = 1.8667
k_hve = 2.246
k_fve = 0.2162
k_kve = 2.924
k_bve = 0.0547
k_ove = 0.7451
k_rbcplas = 0.00128
k_plasrbc = 0.000348
k_bind_in = 0.001015
k_bind_out = 0.000895

f_unb = 0.05
f_hem = 0.45
one_sub_f_hem = 1 - f_hem
one_by_f_hem = 1/f_hem
one_by_one_sub_f_hem = 1/one_sub_f_hem

F_li = 0.45
V_li = 1.80
f_li = 0.16
F_l = 5.60
V_l = 0.53
f_l = 0.30
F_t = 0.03
V_t = 0.2
f_t = 0.05
F_g = 1.13
V_g = 1.13
F_m = 0.59
V_m = 28.0
f_m = 0.03
F_s = 0.02
V_s = 0.18
f_s = 0.20
F_h = 0.26
V_h = 0.33
f_h = 0.02
F_f = 0.74
V_f = 15.0
f_f = 0.03
F_k = 1.24
V_k = 0.31
f_k = 0.24
F_b = 0.78
V_b = 1.40
f_b = 0.04
F_o = 0.36
V_o = 15.8
f_o = 0.05
F_tot = 5.60
V_ven = 3.318
V_art = 2.212

V_lv = 0.159
V_le = 0.371
V_gv = 1.27
V_bv = 0.056
V_be = 1.344
V_sv = 0.036
V_se = 0.144
V_liv = 0.288
V_lie = 1.512
V_kv = 0.0744
V_ke = 0.2356
V_mv = 0.84
V_me = 27.16
V_fv = 0.45
V_fe = 14.55
V_tv = 0.01
V_te = 0.19
V_hv = 0.0066
V_he = 0.3234
V_ov = 0.79
V_oe = 15.01
