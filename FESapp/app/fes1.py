"""
Fuzzy logic one to get the dose based on toxicity and number of cells only which will be modified by FES2 later on based on patient's BSA
"""

def fes1(no_of_cells, toxicity):

    dose = 10

    import numpy as np
    from numpy import log as ln
    import skfuzzy as fuzz
    from skfuzzy import control as ctrl

    #defining Range of Fuzzy inputs and outputs
    fuzzTumorSize = ctrl.Antecedent(np.arange(-3, 12, 1), 'tumor_size')
    fuzzToxicity = ctrl.Antecedent(np.arange(0, 121, 1), 'toxicity')
    fuzzDose = ctrl.Consequent(np.arange(10, 51, 1), 'dose')

    #defining membership functions for fuzzy inputs
    #For Tumor Size
    fuzzTumorSize['VS'] = fuzz.gbellmf(fuzzTumorSize.universe, .8242, 3.278, -2.5)
    fuzzTumorSize['S'] = fuzz.gbellmf(fuzzTumorSize.universe, .8242, 3.278, 0)
    fuzzTumorSize['M'] = fuzz.gbellmf(fuzzTumorSize.universe, .8242, 3.278, 3)
    fuzzTumorSize['B'] = fuzz.gbellmf(fuzzTumorSize.universe, .8242, 3.278, 6.5)
    fuzzTumorSize['VB'] = fuzz.gbellmf(fuzzTumorSize.universe, .8242, 3.278, 10)

    #For Toxicity
    fuzzToxicity['VL'] = fuzz.gaussmf(fuzzToxicity.universe, 0.15, 13.89)
    fuzzToxicity['L']= fuzz.trimf(fuzzToxicity.universe, [0, 30, 60])
    fuzzToxicity['M']= fuzz.trimf(fuzzToxicity.universe, [30, 60, 90])
    fuzzToxicity['H']= fuzz.trimf(fuzzToxicity.universe, [60, 90, 120])
    fuzzToxicity['VH']= fuzz.gaussmf(fuzzToxicity.universe, 117.2, 2.421)

    #For Dose
    fuzzDose['VVL'] = fuzz.trimf(fuzzDose.universe, [3.333, 10, 16.67])
    fuzzDose['VL'] = fuzz.trimf(fuzzDose.universe, [10, 16.67, 23.33])
    fuzzDose['L'] = fuzz.trimf(fuzzDose.universe, [16.67, 23.33, 30])
    fuzzDose['M'] = fuzz.trimf(fuzzDose.universe, [23.33, 30, 36.67])
    fuzzDose['H'] = fuzz.trimf(fuzzDose.universe, [30, 36.67, 43.33])
    fuzzDose['VH'] = fuzz.trimf(fuzzDose.universe, [36.67, 43.33, 50])
    fuzzDose['VVH'] = fuzz.trimf(fuzzDose.universe, [43.33, 50, 56.67])

    #Defining rules for fuzzy logic
    #For tumor size = Very Small
    rule1 = ctrl.Rule(fuzzTumorSize['VS'], fuzzDose['VVL'])

    #For tumor size = Small
    rule2 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VL'], fuzzDose['M'])
    rule3 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['L'], fuzzDose['L'])
    rule4 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['M'], fuzzDose['VL'])
    rule5 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['H'], fuzzDose['VVL'])
    rule6 = ctrl.Rule(fuzzTumorSize['S'] & fuzzToxicity['VH'], fuzzDose['VVL'])

    #For tumor size = Medium
    rule7 = ctrl.Rule(fuzzTumorSize['M'] & fuzzToxicity['VL'], fuzzDose['H'])
    rule8 = ctrl.Rule(fuzzTumorSize['M'] & fuzzToxicity['L'], fuzzDose['L'])
    rule9 = ctrl.Rule(fuzzTumorSize['M'] & fuzzToxicity['M'], fuzzDose['L'])
    rule10 = ctrl.Rule(fuzzTumorSize['M'] & fuzzToxicity['H'], fuzzDose['VL'])
    rule11 = ctrl.Rule(fuzzTumorSize['M'] & fuzzToxicity['VH'], fuzzDose['VVL'])

    #For tumor size = Big
    rule12 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['VL'], fuzzDose['VH'])
    rule13 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['L'], fuzzDose['H'])
    rule14 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['M'], fuzzDose['M'])
    rule15 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['H'], fuzzDose['L'])
    rule16 = ctrl.Rule(fuzzTumorSize['B'] & fuzzToxicity['VH'], fuzzDose['VL'])

    #For tumor size = Very Big
    rule17 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['VL'], fuzzDose['VH'])
    rule18 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['L'], fuzzDose['VH'])
    rule19 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['M'], fuzzDose['H'])
    rule20 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['H'], fuzzDose['M'])
    rule21 = ctrl.Rule(fuzzTumorSize['VB'] & fuzzToxicity['VH'], fuzzDose['L'])

    #Simulation
    dose_ctrl = ctrl.ControlSystem([rule1, rule2, rule3, rule4, rule5,
                                rule6, rule7, rule8, rule9, rule10,
                                rule11, rule12, rule13, rule14, rule15,
                                rule16, rule17, rule18, rule19, rule20, rule21])
    dosing = ctrl.ControlSystemSimulation(dose_ctrl)

    #Calculating & returnng dose
    dosing.input['tumor_size'] = no_of_cells
    dosing.input['toxicity'] = toxicity
    dosing.compute()
    dose = dosing.output['dose']
    #print("Dose: ", dose)

    return dose