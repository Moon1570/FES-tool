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


def dS2dt(t, S):
    C_v, C_rbcv, C_lv, C_le, C_lb, C_art, C_rbca, C_gv, C_bv, C_be, C_bb, C_sv, C_se, C_sb, C_liv, C_lie, C_lib, C_kv, C_ke, C_kb, C_mv, C_me, C_mb, C_fv, C_fe, C_fb, C_tv, C_te, C_tb, C_hv, C_he, C_hb, C_ov, C_oe, C_ob = S
    
    return [ one_by_one_sub_f_hem * (1),
            ((1/tau_g)*ln((ln(rho_g/N_0))/ln(rho_g/(2*N_0)))*N_t*ln(rho_g/N_t)) - (K_eff*C_eff* N_t),
            C_t-(eta*T_t)]


def dSdt(t, S):
    global temp
    #temp = 0.0
    C_t, N_t, T_t = S
    D_t =0    
    flag = {}

    #print(t, int(t), D_t)
    if int(t) == 0:
        D_t = 45.0
        dose[int(t)] = D_t

    elif int(t)==1:
        D_t = 0.0
        dose[int(t)] = D_t

    elif int(t)%14 == 0 :
        if dose[int(t)] != 0:
            D_t = dose[int(t)]
            print(dose[int(t)], t)

        else:
            fes1Dose = fes1.fes1(N_t, T_t) 

            percentIncrease = fes2.fes2(N_t, T_t, BSA)

            D_t = fes1Dose + fes1Dose * percentIncrease
            if D_t>43:
                D_t = D_t*.95

            dose[int(t)] = D_t
            temp = D_t


    elif int(t)%14 == 1:
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
    return render(request, 'home.html', {'name':'Moon'})

def calc(request):
    #get name from request
    name = request.POST.get('name')
    weight = int(request.POST.get('weight'))
    intervalTime = int(request.POST.get('intervalTime'))

    #print to console
    print("Printing - " + name, weight, intervalTime)

    global BSA
    BSA = functions.calcBSA(weight)


    
    C_t_0 = 0
    N_t_0 = N_0
    T_t_0 = 0

    t_eval  = np.linspace(0, 120, 121, dtype=int)
    t_range = (0, 120)


    sol = solve_ivp(dSdt, t_range, y0=[C_t_0, N_t_0, T_t_0],method='LSODA',t_eval=t_eval)
    t = sol.t
    C_t = sol.y[0]
    N_t = sol.y[1]
    T_t = sol.y[2]

    logNT = np.log10(N_t)
    print(dose)
    print(logNT)
    toxPlot = generateFigure.get_plot(T_t, "Toxicity Vs Days", "Day", "Toxicity")

    noCellPlot = generateFigure.get_plot(logNT, "No of cells Vs Days", "Day", "Cells")

    print(N_t[84])
    print(T_t)
    dosePlot = generateFigure.get_plot(dose, "Dose Vs Days", "Day", "Dose")

    


    return render(request, 'result.html', {'name':name, 'weight':weight, 'BSA':BSA, 'toxPlot':toxPlot, 'noCellPlot':noCellPlot, 'dosePlot':dosePlot})

