from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from app import functions, fes1, fes2
from app import generateFigure

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from numpy import log as ln
import skfuzzy as fuzz
from skfuzzy import control as ctrl
#declear a list
dose = [0] * 121


C_0 = 0
N_0 = 10000000000
tau_g = 150
rho_g = 1000000000000
K_eff = 2.7 * 10 ** -2
lambdaa = .27
eta = .4
D_th = 10
D_max = 50
#C_th = 4.1 * 10 **3 #Is this the Concentration Threshold C_cum?
T_max = 100
C_th = 10
t = int(0)
temp = 0.0



#Parameters for PBPK model
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


def dS2dt(t, S):
    C_v, C_rbcv, C_lv, C_le, C_lb, C_art, C_rbca, C_gv, C_bv, C_be, C_bb, C_sv, C_se, C_sb, C_liv, C_lie, C_lib, C_kv, C_ke, C_kb, C_mv, C_me, C_mb, C_fv, C_fe, C_fb, C_tv, C_te, C_tb, C_hv, C_he, C_hb, C_ov, C_oe, C_ob = S
    dose = adminstrated_dose
        
    if(int(dose) != 45 and int(t) == 1):
        dose = dose + adminstrated_dose
    


    return [((((F_l*C_le)+(F_b*C_be)+(F_s*C_se)+(F_li*C_lie)+(F_h*C_he)+(F_k*C_ke)+(F_m*C_me)+(F_f*C_fe)+(F_t*C_te)+(F_o*C_oe)) - (F_tot * C_v))/(V_ven * one_sub_f_hem))+(dose/(V_ven*one_sub_f_hem))+(one_by_one_sub_f_hem*f_hem*k_rbcplas*C_rbcv)-(k_plasrbc*f_unb*C_v), #EQ 1
            (one_sub_f_hem*k_plasrbc*f_unb*C_v/f_hem)-(k_rbcplas*C_rbcv), #EQ 2
            ((F_l/V_lv)*(C_v-C_lv)) - (k_lve*f_unb*C_lv) + ((V_le*k_lev*C_le)/V_lv), #EQ 3
            (V_lv*k_lve*f_unb*C_lv/V_le)-(k_lev*C_le)+(k_bind_out*C_lb)-(k_bind_in*C_le), #EQ 4
            (k_bind_in*C_le)-(k_bind_out*C_lb), #EQ 5
            (((F_l*C_lv)-(F_tot*C_art))/(V_art*one_sub_f_hem))+(f_hem*k_rbcplas*C_rbca/one_sub_f_hem)-(k_plasrbc*f_unb*C_art), #EQ 6
            (one_sub_f_hem*k_rbcplas*f_unb*C_art/f_hem) - (k_rbcplas*C_rbca), #EQ 7
            ((F_g/V_gv)*(C_art - C_gv)) , #EQ 8
            ((F_b/V_bv)*(C_art - C_bv)) - (k_bve*f_unb*C_bv) + ((V_be*k_bev*C_be)/V_bv), #EQ 9
            (V_bv*k_bve*f_unb*C_bv/V_be)-(k_bev*C_be)+(k_bind_out*C_bb)-(k_bind_in*C_be), #EQ 10\
            (k_bind_in*C_be)-(k_bind_out*C_bb), #EQ 11
            ((F_s/V_sv)*(C_art - C_sv)) - (k_sve*f_unb*C_sv) + ((V_se*k_sev*C_se)/V_sv), #EQ 12
            (V_sv*k_sve*f_unb*C_sv/V_se)-(k_sev*C_se)+(k_bind_out*C_sb)-(k_bind_in*C_se), #EQ 13
            (k_bind_in*C_se)-(k_bind_out*C_sb), #EQ 14
            ((F_li*C_art)+(F_g*C_gv)+(F_s*C_sv)-(C_liv*(F_g+F_s+F_li))/V_liv)-(k_live*f_unb*C_liv)+(V_lie*k_liev*C_lie/V_liv), #EQ 15
            (V_liv*k_live*f_unb*C_liv/V_lie)-(k_liev*C_lie)+(k_bind_out*C_lib)-(k_bind_in*C_lie)-(k_clli*C_lie), #EQ 16
            (k_bind_in*C_lie)-(k_bind_out*C_lib)-(k_clli*C_lib), #EQ 17
            (F_h*(C_art-C_hv)/V_hv)-(k_hve*f_unb*C_hv)+(V_he*k_hev*C_he/V_hv), #EQ 18
            (V_hv*k_hve*f_unb*C_hv/V_he)-(k_hev*C_he)+(k_bind_out*C_hb)-(k_bind_in*C_he), #EQ 19
            (k_bind_in*C_he)-(k_bind_out*C_hb), #EQ 20
            ((F_k/V_kv)*(C_art - C_kv)) - (k_kve*f_unb*C_kv) + ((V_ke*k_kev*C_ke)/V_kv), #EQ 21
            (V_kv*k_kve*f_unb*C_kv/V_ke)-(k_kev*C_ke)+(k_bind_out*C_kb)-(k_bind_in*C_ke), #EQ 22
            (k_bind_in*C_ke)-(k_bind_out*C_kb), #EQ 23
            ((F_m/V_mv)*(C_art - C_mv)) - (k_mve*f_unb*C_mv) + ((V_me*k_mev*C_me)/V_mv), #EQ 24
            (V_mv*k_mve*f_unb*C_mv/V_me)-(k_mev*C_me)+(k_bind_out*C_mb)-(k_bind_in*C_me), #EQ 25
            (k_bind_in*C_me)-(k_bind_out*C_mb), #EQ 26
            (F_f*(C_art-C_fv)/V_fv)-(k_fve*f_unb*C_fv)+(V_fe*k_fev*C_fe/V_fv), #EQ 27
            (V_fv*k_fve*f_unb*C_fv/V_fe)-(k_fev*C_fe)+(k_bind_out*C_fb)-(k_bind_in*C_fe), #EQ 28
            (k_bind_in*C_fe)-(k_bind_out*C_fb), #EQ 29
            (F_t*(C_art-C_tv)/V_tv)-(k_tve*f_unb*C_tv)+(V_te*k_tev*C_te/V_tv), #EQ 30
            (V_tv*k_tve*f_unb*C_tv/V_te)-(k_tev*C_te)+(k_bind_out*C_tb)-(k_bind_in*C_te), #EQ 31
            (k_bind_in*C_te)-(k_bind_out*C_tb), #EQ 32
            (F_o*(C_art-C_ov)/V_ov)-(k_ove*f_unb*C_ov)+(V_oe*k_oev*C_oe/V_ov), #EQ 33
            (V_ov*k_ove*f_unb*C_ov/V_oe)-(k_oev*C_oe)+(k_bind_out*C_ob)-(k_bind_in*C_oe), #EQ 34
            (k_bind_in*C_oe)-(k_bind_out*C_ob) #EQ 35
            ]


def dSdt(t, S):
    global temp
    #temp = 0.0
    C_t, N_t, T_t = S
    D_t =0    
    flag = {}

    
    if int(t) == 0:
        D_t = 45.0
        dose[int(t)] = D_t

    elif int(t)==1:
        D_t = 0.0
        dose[int(t)] = D_t

    elif int(t)%intval_time == 0 :
        if dose[int(t)] != 0:
            D_t = dose[int(t)]

        else:
            fes1Dose = fes1.fes1(N_t, T_t) 

            percentIncrease = fes2.fes2(N_t, T_t, BSA)

            D_t = fes1Dose + fes1Dose * percentIncrease
            if D_t>43:
                D_t = D_t*.95

            dose[int(t)] = D_t
            temp = D_t


    elif int(t)%intval_time == 1:
        fes1Dose = fes1.fes1(N_t, T_t) 
        percentIncrease = fes2.fes2(N_t, T_t, BSA)

        D_t = dose[int(t-1)]

        if D_t>43:
            D_t = D_t*.95

        dose[int(t)] = D_t



        #D_t = temp
        #dose[int(t)] = D_t
       # print(D_t, t, int(t), "Cooled", temp)

    else:
    #  print("false", t, int(t))
        dose[int(t)] = 0
        D_t = dose[int(t)]
        #print(t, int(t), C_t, N_t, T_t, D_t)


    if C_t>=C_th:
        H_ct_cth = 1
    else:
        H_ct_cth = 0
        
    C_eff = (C_t-C_th)*H_ct_cth

   # print(D_t, t)
    
   # return [ 40-.27*C_t,
   #         ((1/tau_g)*ln((ln(rho_g/N_0))/ln(rho_g/(2*N_0)))*N_t*ln(rho_g/N_t)) - (K_eff*C_eff* N_t),
   #         C_t-(eta*T_t)]

    #print(K_eff, C_eff, N_t, t)
    return [ D_t-lambdaa*C_t,
            ((1/tau_g)*ln((ln(rho_g/N_0))/ln(rho_g/(2*N_0)))*N_t*ln(rho_g/N_t)) - (K_eff*C_eff* N_t),
            C_t-(eta*T_t)]


# Create your views here.
def home(request):
    return render(request, 'home.html')

def calc(request):
    
    global intval_time
    #get name from request
    name = request.POST.get('name')
    weight = int(request.POST.get('weight'))
    intval_time = int(request.POST.get('intervalTime'))
    print(intval_time)
    #print to console

    global BSA
    BSA = functions.calcBSA(weight)


    
    C_t_0 = 0
    N_t_0 = N_0
    T_t_0 = 0


    C_v_0, C_rbcv_0, C_lv_0, C_le_0, C_lb_0, C_art_0, C_rbca_0, C_gv_0, C_bv_0, C_be_0, C_bb_0, C_sv_0, C_se_0, C_sb_0, C_liv_0, C_lie_0, C_lib_0, C_kv_0, C_ke_0, C_kb_0, C_mv_0, C_me_0, C_mb_0, C_fv_0, C_fe_0, C_fb_0, C_tv_0, C_te_0, C_tb_0, C_hv_0, C_he_0, C_hb_0, C_ov_0, C_oe_0, C_ob_0 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


    t_eval  = np.linspace(0, 120, 121, dtype=int)
    t_range = (0, 120)

    


    sol = solve_ivp(dSdt, t_range, y0=[C_t_0, N_t_0, T_t_0],method='LSODA',t_eval=t_eval)
    t = sol.t
    C_t = sol.y[0]
    N_t = sol.y[1]
    T_t = sol.y[2]

    logNT = np.log10(N_t)
    #print(dose)
    #print(logNT)
    toxPlot = generateFigure.get_plot(T_t, "Toxicity Vs Days", "Day", "Toxicity")

    noCellPlot = generateFigure.get_plot(logNT, "No of cells Vs Days", "Day", "Cells")

    #print(N_t[84])
    print(T_t)
    dosePlot = generateFigure.get_plot(dose, "Dose Vs Days", "Day", "Dose")

    uniqueDose = []
    for i in range(0, 121):
        if i%intval_time == 0:
            uniqueDose.append(dose[i])

    print(uniqueDose)

    

    t_eval2  = np.linspace(0, 2, 48, dtype=float)
    t_range2 = (0, 2)
    
    global adminstrated_dose
    pbpk = []
    for dt in uniqueDose:
        adminstrated_dose = dt
        sol2 = solve_ivp(dS2dt, t_range2, y0=[C_v_0, C_rbcv_0, C_lv_0, C_le_0, C_lb_0, C_art_0, C_rbca_0, C_gv_0, C_bv_0, C_be_0, C_bb_0, C_sv_0, C_se_0, C_sb_0, C_liv_0, C_lie_0, C_lib_0, C_kv_0, C_ke_0, C_kb_0, C_mv_0, C_me_0, C_mb_0, C_fv_0, C_fe_0, C_fb_0, C_tv_0, C_te_0, C_tb_0, C_hv_0, C_he_0, C_hb_0, C_ov_0, C_oe_0, C_ob_0],method = 'LSODA',t_eval=t_eval2, atol=1.49e-8, rtol=1.49e-8)
        pbpk.append(sol2)
       
    #sol2 = solve_ivp(dS2dt, t_range2, y0=[C_v_0, C_rbcv_0, C_lv_0, C_le_0, C_lb_0, C_art_0, C_rbca_0, C_gv_0, C_bv_0, C_be_0, C_bb_0, C_sv_0, C_se_0, C_sb_0, C_liv_0, C_lie_0, C_lib_0, C_kv_0, C_ke_0, C_kb_0, C_mv_0, C_me_0, C_mb_0, C_fv_0, C_fe_0, C_fb_0, C_tv_0, C_te_0, C_tb_0, C_hv_0, C_he_0, C_hb_0, C_ov_0, C_oe_0, C_ob_0],method = 'LSODA',t_eval=t_eval2, atol=1.49e-8, rtol=1.49e-8)
    t = sol2.t
    #N_t = sol2.y[1]
    #T_t = sol2.y[2]


    dose1 = pbpk[0]
    dose2 = pbpk[1]
    dose3 = pbpk[2]
    dose4 = pbpk[3]
    dose5 = pbpk[4]
    dose6 = pbpk[5]
    dose7 = pbpk[6]
    dose8 = pbpk[7]
    dose9 = pbpk[8]

    


    venousBloodPlot1 = generateFigure.get_plot(dose1.y[0], "1st Dose Venous Blood Vs Time", "Time", "Venous Blood")
    venousBloodPlot2 = generateFigure.get_plot(dose2.y[0], "2nd Dose Venous Blood Vs Time", "Time", "Venous Blood")
    venousBloodPlot3 = generateFigure.get_plot(dose3.y[0], "3rd Dose Venous Blood Vs Time", "Time", "Venous Blood")

    return render(request, 'result.html', {'name':name, 'weight':weight, 'BSA':BSA, 'toxPlot':toxPlot, 'noCellPlot':noCellPlot, 'dosePlot':dosePlot, 'venousBloodPlot1':venousBloodPlot1, 'venousBloodPlot2':venousBloodPlot2, 'venousBloodPlot3':venousBloodPlot3})

