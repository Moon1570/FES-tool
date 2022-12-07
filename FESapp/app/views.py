from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from app import functions
from app import generateFigure

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from numpy import log as ln

C_0 = 0
N_0 = 10** 10
tau_g = 150
rho_g = 10 ** 12
K_eff = 2.7 * 10 ** -2
lambdaa = .27
eta = .4
D_th = 10
D_max = 50
#C_th = 4.1 * 10 **3 #Is this the Concentration Threshold C_cum?
T_max = 100
C_th = 10

def dSdt(t, S):
    C_t, N_t, T_t = S
        

    if int(t)%15 == 0 or (int(t)+1)%15 ==0 or t==0:
    # print("true", t, int(t))
        D_t = 40
    else:
    #  print("false", t, int(t))
        D_t = 0
        
    if C_t>=C_th:
        H_ct_cth = 1
    else:
        H_ct_cth = 0
        
    C_eff = (C_t-C_th)*H_ct_cth

    
   # return [ 40-.27*C_t,
   #         ((1/tau_g)*ln((ln(rho_g/N_0))/ln(rho_g/(2*N_0)))*N_t*ln(rho_g/N_t)) - (K_eff*C_eff* N_t),
   #         C_t-(eta*T_t)]

    return [ D_t-lambdaa*C_t,
            ((1/tau_g)*ln((ln(rho_g/N_0))/ln(rho_g/(2*N_0)))*N_t*ln(rho_g/N_t)) - (.27*1* N_t),
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

    print(t)
    chart = generateFigure.get_plot(T_t)

    


    return render(request, 'result.html', {'name':name, 'weight':weight, 'BSA':BSA, 'chart':chart})

##Martins Model
    #Initializing Variables
    


    #Packaging Equations to a function
    
