"""
Fuzzy logic one to get the dose based on toxicity and number of cells only which will be modified by FES2 later on based on patient's BSA
"""

def fes2(no_of_cells, toxicity, calculated_dose):

    percent_increase = 0.0

    import numpy as np
    from numpy import log as ln
    import skfuzzy as fuzz
    from skfuzzy import control as ctrl

    #defining Range of Fuzzy inputs and outputs
    fuzzTumorSize = ctrl.Antecedent(np.arange(-3, 12, 1), 'tumor_size')
    fuzzToxicity = ctrl.Antecedent(np.arange(0, 121, 1), 'toxicity')
    fuzzCalculatedDose = ctrl.Antecedent(np.arange(0, 301, 10), 'calculated_dose')
    fuzzPercentDoseIncrease = ctrl.Consequent(np.arange(0, .9, .1), 'percent_dose_increase')

    #defining membership functions for fuzzy inputs
    #For Tumor Size
    fuzzTumorSize['VS'] = fuzz.gaussmf(fuzzTumorSize.universe, -1.6, 2.2)
    fuzzTumorSize['S'] = fuzz.trimf(fuzzTumorSize.universe, [-.5, 2.25, 5])
    fuzzTumorSize['B'] = fuzz.trimf(fuzzTumorSize.universe, [3, 5.75, 8.5])
    fuzzTumorSize['VB'] = fuzz.gaussmf(fuzzTumorSize.universe, 9.85, 1.2)

    #For Toxicity
    fuzzToxicity['VL'] = fuzz.trapmf(fuzzToxicity.universe, [-27, -3, 3, 30])
    fuzzToxicity['L']= fuzz.trimf(fuzzToxicity.universe, [0, 30, 60])
    fuzzToxicity['M']= fuzz.trimf(fuzzToxicity.universe, [30, 60, 90])
    fuzzToxicity['H']= fuzz.trimf(fuzzToxicity.universe, [60, 90, 120])
    fuzzToxicity['VH']= fuzz.trapmf(fuzzToxicity.universe, [90, 117, 123, 147])

    #For Calculated Dose
    fuzzCalculatedDose['LD'] = fuzz.trapmf(fuzzCalculatedDose.universe, [0, 0, 165, 180])
    fuzzCalculatedDose['ND'] = fuzz.trimf(fuzzCalculatedDose.universe, [165, 180, 202])
    fuzzCalculatedDose['OD'] = fuzz.trimf(fuzzCalculatedDose.universe, [184, 202, 216])
    fuzzCalculatedDose['VOD'] = fuzz.trapmf(fuzzCalculatedDose.universe, [200, 245, 300, 300])

    #For Percent Dose Increase
    fuzzPercentDoseIncrease['normal'] = fuzz.trapmf(fuzzPercentDoseIncrease.universe,  [0, 0, 0, 0.1])
    fuzzPercentDoseIncrease['low_increase'] = fuzz.trimf(fuzzPercentDoseIncrease.universe, [0.05, 0.15, 0.25])
    fuzzPercentDoseIncrease['increase'] = fuzz.trimf(fuzzPercentDoseIncrease.universe, [0.2, 0.3, 0.4])
    fuzzPercentDoseIncrease['very_increase'] = fuzz.trapmf(fuzzPercentDoseIncrease.universe, [0.35, 0.45, 0.8, 0.8])

    #Defining rules for fuzzy logic
    #For tumor size = Very Small
    rule1 = ctrl.Rule(fuzzTumorSize['VS'] & fuzzCalculatedDose['LD'], fuzzPercentDoseIncrease['normal'])
    rule2 = ctrl.Rule(fuzzTumorSize['VS'] & fuzzCalculatedDose['ND'], fuzzPercentDoseIncrease['low_increase'])
    rule3 = ctrl.Rule(fuzzTumorSize['VS'] & fuzzCalculatedDose['OD'], fuzzPercentDoseIncrease['increase'])
    rule4 = ctrl.Rule(fuzzTumorSize['VS'] & fuzzCalculatedDose['VOD'], fuzzPercentDoseIncrease['increase'])

    #For tumor size =  Small
    rule5 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VL'] & fuzzCalculatedDose['LD'], fuzzPercentDoseIncrease['normal'])
    rule6 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VL'] & fuzzCalculatedDose['ND'], fuzzPercentDoseIncrease['low_increase'])
    rule7 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VL'] & fuzzCalculatedDose['OD'], fuzzPercentDoseIncrease['increase'])
    rule8 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['L'] & fuzzCalculatedDose['LD'], fuzzPercentDoseIncrease['normal'])
    rule9 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['L'] & fuzzCalculatedDose['ND'], fuzzPercentDoseIncrease['low_increase'])
    rule10 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['L'] & fuzzCalculatedDose['OD'], fuzzPercentDoseIncrease['increase'])
    rule11 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['M'] & fuzzCalculatedDose['LD'], fuzzPercentDoseIncrease['normal'])
    rule12 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['M'] & fuzzCalculatedDose['ND'], fuzzPercentDoseIncrease['normal'])
    rule13 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['M'] & fuzzCalculatedDose['OD'], fuzzPercentDoseIncrease['increase'])
    rule14 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['H'], fuzzPercentDoseIncrease['normal'])
    rule15 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VH'], fuzzPercentDoseIncrease['normal'])
    rule16 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VL'] & fuzzCalculatedDose['VOD'], fuzzPercentDoseIncrease['increase'])
    rule17 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['L'] & fuzzCalculatedDose['VOD'], fuzzPercentDoseIncrease['increase'])
    rule18 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['M'] & fuzzCalculatedDose['VOD'], fuzzPercentDoseIncrease['increase'])

    #For tumor size =  Big
    rule19 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['VL'], fuzzPercentDoseIncrease['normal'])
    rule20 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['L'], fuzzPercentDoseIncrease['normal'])
    rule21 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['M'], fuzzPercentDoseIncrease['normal'])
    rule22 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['H'], fuzzPercentDoseIncrease['normal'])
    rule23 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['VH'], fuzzPercentDoseIncrease['normal'])

    #For tumor size =  Very Big
    rule24 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['VL'], fuzzPercentDoseIncrease['normal'])
    rule25 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['L'], fuzzPercentDoseIncrease['normal'])
    rule26 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['M'] & fuzzCalculatedDose['LD'], fuzzPercentDoseIncrease['normal'])
    rule27 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['M'] & fuzzCalculatedDose['ND'], fuzzPercentDoseIncrease['normal'])
    rule28 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['M'] & fuzzCalculatedDose['VOD'], fuzzPercentDoseIncrease['increase'])
    rule29 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['H'], fuzzPercentDoseIncrease['normal'])
    rule30 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['VH'], fuzzPercentDoseIncrease['normal'])
    #rule31 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['M'] & fuzzCalculatedDose['OD'], fuzzPercentDoseIncrease['increase'])


    #Simulation
    percent_dose_increase_ctrl = ctrl.ControlSystem([rule1, rule2, rule3, rule4, rule5,
                                rule6, rule7, rule8, rule9, rule10,
                                rule11, rule12, rule13, rule14, rule15,
                                rule16, rule17, rule18, rule19, rule20, 
                                rule21, rule22, rule23, rule24, rule25,
                                rule26, rule27, rule28, rule29, rule30])
    increasing = ctrl.ControlSystemSimulation(percent_dose_increase_ctrl)

    increasing.input['tumor_size'] = no_of_cells
    increasing.input['toxicity'] = toxicity
    increasing.input['calculated_dose'] = calculated_dose
    increasing.compute()
    percent_increase = increasing.output['percent_dose_increase']

    return percent_increase