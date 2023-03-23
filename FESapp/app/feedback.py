


def checkVenousBlood(C_ven, C_rbcv, max_dose):
    
    total = []
    for i in range(0, 48):
        total.append(C_ven[i] + C_rbcv[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkLung(C_lv, C_le, C_lb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_lv[i] + C_le[i] + C_lb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkArterialBlood(C_art, C_rbca, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_art[i] + C_rbca[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkGut(C_gv, max_dose):
    if max(C_gv) > max_dose:
        return True
    else:
        return False

def checkBrain(C_bv, C_be, C_bb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_bv[i] + C_be[i] + C_bb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkSpleen(C_sv, C_se, C_sb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_sv[i] + C_se[i] + C_sb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkLiver(C_liv, C_lie, C_lib, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_liv[i] + C_lie[i] + C_lib[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkHeart(C_hv, C_he, C_hb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_hv[i] + C_he[i] + C_hb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkKidney(C_kv, C_ke, C_kb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_kv[i] + C_ke[i] + C_kb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkMuscle(C_mv, C_me, C_mb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_mv[i] + C_me[i] + C_mb[i])

    print("Max: ", max(total), "Max Dose: ", max_dose)

    if max(total) > max_dose:
        return True
    else:
        return False

def checkFat(C_fv, C_fe, C_fb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_fv[i] + C_fe[i] + C_fb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkTumor(C_tv, C_te, C_tb, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_tv[i] + C_te[i] + C_tb[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def checkOther(C_ov, C_oe, C_ob, max_dose):
    total = []
    for i in range(0, 48):
        total.append(C_ov[i] + C_oe[i] + C_ob[i])

    if max(total) > max_dose:
        return True
    else:
        return False

def check(dose, max_dose):
    venBloodOverDose = checkVenousBlood(dose.y[0], dose.y[1], max_dose[0])
    lungOverDose = checkLung(dose.y[2], dose.y[3], dose.y[4], max_dose[1])
    artBloodOverDose = checkArterialBlood(dose.y[5], dose.y[6], max_dose[2])
    gutOverDose = checkGut(dose.y[7], max_dose[3])
    brainOverDose = checkBrain(dose.y[8], dose.y[9], dose.y[10], max_dose[4])
    spleenOverDose = checkSpleen(dose.y[11], dose.y[12], dose.y[13], max_dose[5])
    liverOverDose = checkLiver(dose.y[14], dose.y[15], dose.y[16], max_dose[6])
    heartOverDose = checkHeart(dose.y[17], dose.y[18], dose.y[19], max_dose[7])
    kidneyOverDose = checkKidney(dose.y[20], dose.y[21], dose.y[22], max_dose[8])
    muscleOverDose = checkMuscle(dose.y[23], dose.y[24], dose.y[25], max_dose[9])
    fatOverDose = checkFat(dose.y[26], dose.y[27], dose.y[28], max_dose[10])
    tumorOverDose = checkTumor(dose.y[29], dose.y[30], dose.y[31], max_dose[11])
    otherOverDose = checkOther(dose.y[32], dose.y[33], dose.y[34], max_dose[12])


    if venBloodOverDose == True or lungOverDose == True or artBloodOverDose == True or gutOverDose == True or brainOverDose == True or spleenOverDose == True or liverOverDose == True or heartOverDose == True or kidneyOverDose == True or muscleOverDose == True or fatOverDose == True or tumorOverDose == True or otherOverDose == True:
        return True
    else:
        return False