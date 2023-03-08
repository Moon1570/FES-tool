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

