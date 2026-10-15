#!python
#
import numpy as np
import openpyxl
import symdata
import symmetry
import symmetry_stress as symstr
import elastic3rd.esutils as esutils

def write_strainmode(xlsname, CrystalType, Ord, flag_stress = 0, StrainIn=None):
    #CijkInd = symdata.coef_ind(CrystalType, Ord)
    #CoefCoef = symdata.coef_crystal(CrystalType, Ord)  
    if flag_stress:
        if StrainIn is None:
            StrainModeCoefAll, StrainMode = symstr.gen_strain_mode(CrystalType, Ord)
        else:
            (Cijk, StrainModeCoefAll, StrainMode) = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
        ExtName = "-StressMode.xlsx"
    else:
        if StrainIn is None:
            StrainModeCoefAll, StrainMode = symmetry.gen_strain_mode(CrystalType, Ord)
        else:
            (Cijk, StrainModeCoefAll, StrainMode) = symmetry.CoefForSingleMode(CrystalType, Ord, StrainIn)
        ExtName = "-EnergyMode.xlsx"
    dest_filename = xlsname.upper() + "-Ord" + str(Ord) + ExtName

    wb = openpyxl.Workbook()
    ws1 = wb.active
    ws1.title = "StrainMode"
    for i in range(0, len(StrainMode)):
        for j in range(0, 3):
            StrainModeij = StrainMode[i][j]
            #print StrainModeij
            ws1.append(StrainModeij.tolist())
        ws1.append([])

    for i in range(2 ,Ord + 1):
        CijInd = symdata.coef_ind(CrystalType, i)
        CijUniq = np.array(symdata.get_unique(CijInd))
        Cijk = symmetry.num2cijk(CijUniq, i)
        StrainModeCoef = eval("StrainModeCoefAll.coef" + str(i))
        SheetName = "StrainModeCoef" + str(i)
        wsi = wb.create_sheet(title=SheetName)
        wsi.append(Cijk.tolist())
        wsi.append([])

        for i in range(0, len(StrainModeCoef)):
            StrainModeCoeflst = StrainModeCoef[i].tolist()
            SMCoefTot = ""
            for k in range(0, len(Cijk)):
                if StrainModeCoef[i][k] != 0:
                    if StrainModeCoef[i][k] == 1:
                        CoefBfCijk = ""
                    else:
                        CoefBfCijk = str(StrainModeCoef[i][k])
                    if SMCoefTot == "" or StrainModeCoef[i][k] < 0:
                        SMCoefTot = SMCoefTot + CoefBfCijk + "C" + str(Cijk[k])
                    else:
                        SMCoefTot = SMCoefTot + "+" + CoefBfCijk + "C" + str(Cijk[k])
            StrainModeCoeflst.append(SMCoefTot)
            wsi.append(StrainModeCoeflst)
            wsi.append([])

    wb.save(filename = dest_filename)


#for test
#strain = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])
#strain = np.eye(3)
CrystalType = "c1"
Ord = 5
flag_stress = 1

#StrainAll = symmetry.gen_strain(1)
#print StrainAll
#print len(StrainAll)
#StrainIn = esutils.read_strainmode('STRAINMODE')
StrainIn = np.array([[0.90, -0.20,  0.04,  0.40,  0.50,  0.00], [0.30,  0.80, -0.12,  0.40,  0.50,  0.00],
    [0.40,  0.20,  0.03,  1.00,  0.10,  0.00]])
write_strainmode(CrystalType, CrystalType, Ord, flag_stress, StrainIn=StrainIn)