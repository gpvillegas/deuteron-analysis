import pickle as P
import R_Barlow as RB
import LT.box as B

"""
CALCULATE SYSTEMATICS: this code uses the class R_Barlow to perform the 
Barlow test on given systematics.

CUT SYSTEMATICS
1. Nominal Cuts
 Window cut: hms_delta; xmin = -8.0, xmax = 8.0
 Window cut: shms_delta; xmin = -10.0, xmax = 22.0
 Window cut: shms_calPID; xmin = 0.7, xmax = 1.3
 Window cut: Em_cut; xmin = -0.05, xmax = 0.05
 Window cut: Q2_cut; xmin = 4.0, xmax = 5.0
 Value cut: noEDTM; value = 0.0
 HMS collimator cut ON
 BCM1_cut
 ztar_diff_cut
 epCoinTime_cut

2. Varying Cuts *DONE
    *1 - HCOLL: HMS Collimator Cut - vary scale
        A) 1.1, 1.0 (Nominal),  B) 0.9,  C) 0.5
    *2 - shms_calPID: SHMS Calorimeter Cut
        ON (Nominal) / OFF  
    *3 - Em_cut: Missing Energy Cut
        A) (-0.08,0.08), B) (-0.07,0.07), C) (-0.06,0.06), (-0.05,0.05) (Nominal),
        D) (-0.04,0.04), E) (-0.03,0.03)
    *4 - ZT: Ztar Difference (HMS Ztar - SHMS Ztar) cut - cut around mean value
        A) (-2.5,2.5), (-2,2) (Nominal), B) (-1.5,1.5), C) (-1,1)
    *5 - epCoinTime_cut: e,p Coincidence Time Cut
        ON (Nominal) / OFF
    *6 - noEDTM: EDTM TDC Signal Cut 
        ON (Nominal) / OFF
    *7 - BCM1_cut: BCM1 Current Cut 
        ON (Nominal) / OFF
        
Other effects tested here:
    1. Radiative correction factor model: FSI/PWIA and JML/CDB dependence
    2. Vertex correction model: FSI/PWIA dependence

"""

def load_pickle(pcl_fnames):
    for i in range(len(pcl_fnames)):
        if i == 0:
            with open(pcl_fnames[i],'rb') as f:
                h = P.load(f)
        else:
            with open(pcl_fnames[i],'rb') as f:
                h.update(P.load(f))
    return h

#%% load pickle file with histograms
#Nominal Cuts
f_JML_120 = 'yield_n_xsec/pm_120_data_JML_histos.pcl'
f_JML_580 = 'yield_n_xsec/pm_580_data_JML_histos.pcl'
f_JML_800 = 'yield_n_xsec/pm_800_data_JML_histos.pcl'
f_JML_900 = 'yield_n_xsec/pm_900_data_JML_histos.pcl'

f_CDB_120 = 'yield_n_xsec/pm_120_data_CD-Bonn_histos.pcl'
f_CDB_580 = 'yield_n_xsec/pm_580_data_CD-Bonn_histos.pcl'
f_CDB_800 = 'yield_n_xsec/pm_800_data_CD-Bonn_histos.pcl'
f_CDB_900 = 'yield_n_xsec/pm_900_data_CD-Bonn_histos.pcl'

f_JML = [f_JML_120,f_JML_580,f_JML_800,f_JML_900]
f_CDB = [f_CDB_120,f_CDB_580,f_CDB_800,f_CDB_900]

h_JML = load_pickle(f_JML)
h_CDB = load_pickle(f_CDB)

#HCOLL Variation A
f_HCOLL_A_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_A/pm_120_data_CD-Bonn_histos_HCOLL_A.pcl'
f_HCOLL_A_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_A/pm_580_data_CD-Bonn_histos_HCOLL_A.pcl'
f_HCOLL_A_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_A/pm_800_data_CD-Bonn_histos_HCOLL_A.pcl'
f_HCOLL_A_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_A/pm_900_data_CD-Bonn_histos_HCOLL_A.pcl'

f_HCOLL_A = [f_HCOLL_A_pm120,f_HCOLL_A_pm580,f_HCOLL_A_pm800,f_HCOLL_A_pm900]

h_HCOLL_A = load_pickle(f_HCOLL_A)

#HCOLL Variation B
f_HCOLL_B_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_B/pm_120_data_CD-Bonn_histos_HCOLL_B.pcl'
f_HCOLL_B_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_B/pm_580_data_CD-Bonn_histos_HCOLL_B.pcl'
f_HCOLL_B_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_B/pm_800_data_CD-Bonn_histos_HCOLL_B.pcl'
f_HCOLL_B_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_B/pm_900_data_CD-Bonn_histos_HCOLL_B.pcl'

f_HCOLL_B = [f_HCOLL_B_pm120,f_HCOLL_B_pm580,f_HCOLL_B_pm800,f_HCOLL_B_pm900]

h_HCOLL_B = load_pickle(f_HCOLL_B)

#HCOLL Variation B
f_HCOLL_C_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_C/pm_120_data_CD-Bonn_histos_HCOLL_C.pcl'
f_HCOLL_C_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_C/pm_580_data_CD-Bonn_histos_HCOLL_C.pcl'
f_HCOLL_C_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_C/pm_800_data_CD-Bonn_histos_HCOLL_C.pcl'
f_HCOLL_C_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/HCOLL_C/pm_900_data_CD-Bonn_histos_HCOLL_C.pcl'

f_HCOLL_C = [f_HCOLL_C_pm120,f_HCOLL_C_pm580,f_HCOLL_C_pm800,f_HCOLL_C_pm900]

h_HCOLL_C = load_pickle(f_HCOLL_C)

#shms_calPID OFF
f_pcalPID_OFF_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_OFF/pm_120_data_CD-Bonn_histos_shms_calPID_OFF.pcl'
f_pcalPID_OFF_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_OFF/pm_580_data_CD-Bonn_histos_shms_calPID_OFF.pcl'
f_pcalPID_OFF_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_OFF/pm_800_data_CD-Bonn_histos_shms_calPID_OFF.pcl'
f_pcalPID_OFF_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_OFF/pm_900_data_CD-Bonn_histos_shms_calPID_OFF.pcl'

f_pcalPID_OFF = [f_pcalPID_OFF_pm120,f_pcalPID_OFF_pm580,f_pcalPID_OFF_pm800,
                 f_pcalPID_OFF_pm900]

h_pcalPID_OFF = load_pickle(f_pcalPID_OFF)

#shms_calPID Variation A
f_pcalPID_A_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_A/pm_120_data_CD-Bonn_histos_shms_calPID_A.pcl'
f_pcalPID_A_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_A/pm_580_data_CD-Bonn_histos_shms_calPID_A.pcl'
f_pcalPID_A_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_A/pm_800_data_CD-Bonn_histos_shms_calPID_A.pcl'
f_pcalPID_A_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_A/pm_900_data_CD-Bonn_histos_shms_calPID_A.pcl'

f_pcalPID_A = [f_pcalPID_A_pm120,f_pcalPID_A_pm580,f_pcalPID_A_pm800,
                 f_pcalPID_A_pm900]

h_pcalPID_A = load_pickle(f_pcalPID_A)

# #shms_calPID Variation B
f_pcalPID_B_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_B/pm_120_data_CD-Bonn_histos_shms_calPID_B.pcl'
f_pcalPID_B_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_B/pm_580_data_CD-Bonn_histos_shms_calPID_B.pcl'
f_pcalPID_B_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_B/pm_800_data_CD-Bonn_histos_shms_calPID_B.pcl'
f_pcalPID_B_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_B/pm_900_data_CD-Bonn_histos_shms_calPID_B.pcl'

f_pcalPID_B = [f_pcalPID_B_pm120,f_pcalPID_B_pm580,f_pcalPID_B_pm800,
                 f_pcalPID_B_pm900]

h_pcalPID_B = load_pickle(f_pcalPID_B)

#shms_calPID Variation C
f_pcalPID_C_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_C/pm_120_data_CD-Bonn_histos_shms_calPID_C.pcl'
f_pcalPID_C_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_C/pm_580_data_CD-Bonn_histos_shms_calPID_C.pcl'
f_pcalPID_C_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_C/pm_800_data_CD-Bonn_histos_shms_calPID_C.pcl'
f_pcalPID_C_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/shms_calPID_C/pm_900_data_CD-Bonn_histos_shms_calPID_C.pcl'

f_pcalPID_C = [f_pcalPID_C_pm120,f_pcalPID_C_pm580,f_pcalPID_C_pm800,
                 f_pcalPID_C_pm900]

h_pcalPID_C = load_pickle(f_pcalPID_C)

#Em Variation A
f_Em_A_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/Em_A/pm_120_data_CD-Bonn_histos_Em_A.pcl'
f_Em_A_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/Em_A/pm_580_data_CD-Bonn_histos_Em_A.pcl'
f_Em_A_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/Em_A/pm_800_data_CD-Bonn_histos_Em_A.pcl'
f_Em_A_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/Em_A/pm_900_data_CD-Bonn_histos_Em_A.pcl'

f_Em_A = [f_Em_A_pm120,f_Em_A_pm580,f_Em_A_pm800,f_Em_A_pm900]

h_Em_A = load_pickle(f_Em_A)

#Em Variation B
f_Em_B_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/Em_B/pm_120_data_CD-Bonn_histos_Em_B.pcl'
f_Em_B_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/Em_B/pm_580_data_CD-Bonn_histos_Em_B.pcl'
f_Em_B_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/Em_B/pm_800_data_CD-Bonn_histos_Em_B.pcl'
f_Em_B_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/Em_B/pm_900_data_CD-Bonn_histos_Em_B.pcl'

f_Em_B = [f_Em_B_pm120,f_Em_B_pm580,f_Em_B_pm800,f_Em_B_pm900]

h_Em_B = load_pickle(f_Em_B)

#Em Variation C
f_Em_C_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/Em_C/pm_120_data_CD-Bonn_histos_Em_C.pcl'
f_Em_C_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/Em_C/pm_580_data_CD-Bonn_histos_Em_C.pcl'
f_Em_C_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/Em_C/pm_800_data_CD-Bonn_histos_Em_C.pcl'
f_Em_C_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/Em_C/pm_900_data_CD-Bonn_histos_Em_C.pcl'

f_Em_C = [f_Em_C_pm120,f_Em_C_pm580,f_Em_C_pm800,f_Em_C_pm900]

h_Em_C = load_pickle(f_Em_C)

#Em Variation D
f_Em_D_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/Em_D/pm_120_data_CD-Bonn_histos_Em_D.pcl'
f_Em_D_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/Em_D/pm_580_data_CD-Bonn_histos_Em_D.pcl'
f_Em_D_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/Em_D/pm_800_data_CD-Bonn_histos_Em_D.pcl'
f_Em_D_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/Em_D/pm_900_data_CD-Bonn_histos_Em_D.pcl'

f_Em_D = [f_Em_D_pm120,f_Em_D_pm580,f_Em_D_pm800,f_Em_D_pm900]

h_Em_D = load_pickle(f_Em_D)

#Em Variation E
f_Em_E_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/Em_E/pm_120_data_CD-Bonn_histos_Em_E.pcl'
f_Em_E_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/Em_E/pm_580_data_CD-Bonn_histos_Em_E.pcl'
f_Em_E_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/Em_E/pm_800_data_CD-Bonn_histos_Em_E.pcl'
f_Em_E_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/Em_E/pm_900_data_CD-Bonn_histos_Em_E.pcl'

f_Em_E = [f_Em_E_pm120,f_Em_E_pm580,f_Em_E_pm800,f_Em_E_pm900]

h_Em_E = load_pickle(f_Em_E)

#ZT Variation A
f_ZT_A_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_A/pm_120_data_CD-Bonn_histos_ZT_A.pcl'
f_ZT_A_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_A/pm_580_data_CD-Bonn_histos_ZT_A.pcl'
f_ZT_A_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_A/pm_800_data_CD-Bonn_histos_ZT_A.pcl'
f_ZT_A_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_A/pm_900_data_CD-Bonn_histos_ZT_A.pcl'

f_ZT_A = [f_ZT_A_pm120,f_ZT_A_pm580,f_ZT_A_pm800,f_ZT_A_pm900]

h_ZT_A = load_pickle(f_ZT_A)

#ZT Variation B
f_ZT_B_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_B/pm_120_data_CD-Bonn_histos_ZT_B.pcl'
f_ZT_B_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_B/pm_580_data_CD-Bonn_histos_ZT_B.pcl'
f_ZT_B_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_B/pm_800_data_CD-Bonn_histos_ZT_B.pcl'
f_ZT_B_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_B/pm_900_data_CD-Bonn_histos_ZT_B.pcl'

f_ZT_B = [f_ZT_B_pm120,f_ZT_B_pm580,f_ZT_B_pm800,f_ZT_B_pm900]

h_ZT_B = load_pickle(f_ZT_B)

#ZT Variation C
f_ZT_C_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_C/pm_120_data_CD-Bonn_histos_ZT_C.pcl'
f_ZT_C_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_C/pm_580_data_CD-Bonn_histos_ZT_C.pcl'
f_ZT_C_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_C/pm_800_data_CD-Bonn_histos_ZT_C.pcl'
f_ZT_C_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/ZT_C/pm_900_data_CD-Bonn_histos_ZT_C.pcl'

f_ZT_C = [f_ZT_C_pm120,f_ZT_C_pm580,f_ZT_C_pm800,f_ZT_C_pm900]

h_ZT_C = load_pickle(f_ZT_C)

#epCoinTime OFF
f_CT_OFF_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_OFF/pm_120_data_CD-Bonn_histos_epCoinTime_OFF.pcl'
f_CT_OFF_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_OFF/pm_580_data_CD-Bonn_histos_epCoinTime_OFF.pcl'
f_CT_OFF_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_OFF/pm_800_data_CD-Bonn_histos_epCoinTime_OFF.pcl'
f_CT_OFF_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_OFF/pm_900_data_CD-Bonn_histos_epCoinTime_OFF.pcl'

f_CT_OFF = [f_CT_OFF_pm120,f_CT_OFF_pm580,f_CT_OFF_pm800,f_CT_OFF_pm900]

h_CT_OFF = load_pickle(f_CT_OFF)

#epCoinTime Variation A
f_CT_A_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_A/pm_120_data_CD-Bonn_histos_epCoinTime_A.pcl'
f_CT_A_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_A/pm_580_data_CD-Bonn_histos_epCoinTime_A.pcl'
f_CT_A_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_A/pm_800_data_CD-Bonn_histos_epCoinTime_A.pcl'
f_CT_A_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_A/pm_900_data_CD-Bonn_histos_epCoinTime_A.pcl'

f_CT_A = [f_CT_A_pm120,f_CT_A_pm580,f_CT_A_pm800,f_CT_A_pm900]

h_CT_A = load_pickle(f_CT_A)

#epCoinTime Variation B
f_CT_B_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_B/pm_120_data_CD-Bonn_histos_epCoinTime_B.pcl'
f_CT_B_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_B/pm_580_data_CD-Bonn_histos_epCoinTime_B.pcl'
f_CT_B_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_B/pm_800_data_CD-Bonn_histos_epCoinTime_B.pcl'
f_CT_B_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_B/pm_900_data_CD-Bonn_histos_epCoinTime_B.pcl'

f_CT_B = [f_CT_B_pm120,f_CT_B_pm580,f_CT_B_pm800,f_CT_B_pm900]

h_CT_B = load_pickle(f_CT_B)

#epCoinTime Variation C
f_CT_C_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_C/pm_120_data_CD-Bonn_histos_epCoinTime_C.pcl'
f_CT_C_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_C/pm_580_data_CD-Bonn_histos_epCoinTime_C.pcl'
f_CT_C_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_C/pm_800_data_CD-Bonn_histos_epCoinTime_C.pcl'
f_CT_C_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/epCoinTime_C/pm_900_data_CD-Bonn_histos_epCoinTime_C.pcl'

f_CT_C = [f_CT_C_pm120,f_CT_C_pm580,f_CT_C_pm800,f_CT_C_pm900]

h_CT_C = load_pickle(f_CT_C)

#noEDTM OFF
f_noEDTM_OFF_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/noEDTM_OFF/pm_120_data_CD-Bonn_histos_noEDTM_OFF.pcl'
f_noEDTM_OFF_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/noEDTM_OFF/pm_580_data_CD-Bonn_histos_noEDTM_OFF.pcl'
f_noEDTM_OFF_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/noEDTM_OFF/pm_800_data_CD-Bonn_histos_noEDTM_OFF.pcl'
f_noEDTM_OFF_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/noEDTM_OFF/pm_900_data_CD-Bonn_histos_noEDTM_OFF.pcl'

f_noEDTM_OFF = [f_noEDTM_OFF_pm120,f_noEDTM_OFF_pm580,f_noEDTM_OFF_pm800,f_noEDTM_OFF_pm900]

h_noEDTM_OFF = load_pickle(f_noEDTM_OFF)

#BCM1 OFF
f_BCM1_OFF_pm120 = 'SYSTEMATIC STUDIES/cut_systematics/BCM1_OFF/pm_120_data_CD-Bonn_histos_BCM1_OFF.pcl'
f_BCM1_OFF_pm580 = 'SYSTEMATIC STUDIES/cut_systematics/BCM1_OFF/pm_580_data_CD-Bonn_histos_BCM1_OFF.pcl'
f_BCM1_OFF_pm800 = 'SYSTEMATIC STUDIES/cut_systematics/BCM1_OFF/pm_800_data_CD-Bonn_histos_BCM1_OFF.pcl'
f_BCM1_OFF_pm900 = 'SYSTEMATIC STUDIES/cut_systematics/BCM1_OFF/pm_900_data_CD-Bonn_histos_BCM1_OFF.pcl'

f_BCM1_OFF = [f_BCM1_OFF_pm120,f_BCM1_OFF_pm580,f_BCM1_OFF_pm800,f_BCM1_OFF_pm900]

h_BCM1_OFF = load_pickle(f_BCM1_OFF)

# file locations
fnameBC_pm120_1 = 'yield_n_xsec/bin_centering/acos_phi/pm_120_jml_bc_corr_set_1.txt'
fnameBC_pm120_2 = 'yield_n_xsec/bin_centering/acos_phi/pm_120_jml_bc_corr_set_2.txt'
fnameBC_pm580_1 = 'yield_n_xsec/bin_centering/acos_phi/pm_580_jml_bc_corr_set_1.txt'
fnameBC_pm580_2 = 'yield_n_xsec/bin_centering/acos_phi/pm_580_jml_bc_corr_set_2.txt'
fnameBC_pm800_1 = 'yield_n_xsec/bin_centering/acos_phi/pm_800_jml_bc_corr_set_1.txt'
fnameBC_pm900_1 = 'yield_n_xsec/bin_centering/acos_phi/pm_900_jml_bc_corr_set_1.txt'

# load files with bin centering corrected xsec
f_BC_pm120_1 = B.get_file(fnameBC_pm120_1)
f_BC_pm120_2 = B.get_file(fnameBC_pm120_2)
f_BC_pm580_1 = B.get_file(fnameBC_pm580_1)
f_BC_pm580_2 = B.get_file(fnameBC_pm580_2)
f_BC_pm800_1 = B.get_file(fnameBC_pm800_1)
f_BC_pm900_1 = B.get_file(fnameBC_pm900_1)

s = ['pm_120','pm_580','pm_800','pm_900']

#%% use R_Barlow to calculate R for JML/CDB systematic

JML_xsec = {'pm_120':h_JML['pm_120_JML_FSI_data_xsec_proj'],
            'pm_580':h_JML['pm_580_JML_FSI_data_xsec_proj'],
            'pm_800':h_JML['pm_800_JML_FSI_data_xsec_proj'],
            'pm_900':h_JML['pm_900_JML_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar = RB.R_Barlow(JML_xsec,CDB_xsec,title='Xsec',settings=s)
RBar.calc_R()
# RBar.plot_R()
RBar.plot_R(th_nq=[35,45])

JML_RC = {'pm_120':h_JML['pm_120_JML_FSI_radcorr_factor_proj'],
          'pm_580':h_JML['pm_580_JML_FSI_radcorr_factor_proj'],
          'pm_800':h_JML['pm_800_JML_FSI_radcorr_factor_proj'],
          'pm_900':h_JML['pm_900_JML_FSI_radcorr_factor_proj']}
CDB_RC = {'pm_120':h_CDB['pm_120_CDB_FSI_radcorr_factor_proj'],
          'pm_580':h_CDB['pm_580_CDB_FSI_radcorr_factor_proj'],
          'pm_800':h_CDB['pm_800_CDB_FSI_radcorr_factor_proj'],
          'pm_900':h_CDB['pm_900_CDB_FSI_radcorr_factor_proj']}

RBar = RB.R_Barlow(JML_RC,CDB_RC,title='Radiative Corrections',settings=s,ana_type='factor')
RBar.calc_R()
# RBar.plot_R()

#%% use R_Barlow to calculate R for FSI/PWIA systematic

FSI_xsec = {'pm_120':h_JML['pm_120_JML_FSI_data_xsec_proj'],
            'pm_580':h_JML['pm_580_JML_FSI_data_xsec_proj'],
            'pm_800':h_JML['pm_800_JML_FSI_data_xsec_proj'],
            'pm_900':h_JML['pm_900_JML_FSI_data_xsec_proj']}
PWIA_xsec = {'pm_120':h_JML['pm_120_JML_PWIA_data_xsec_proj'],
            'pm_580':h_JML['pm_580_JML_PWIA_data_xsec_proj'],
            'pm_800':h_JML['pm_800_JML_PWIA_data_xsec_proj'],
            'pm_900':h_JML['pm_900_JML_PWIA_data_xsec_proj']}

RBar = RB.R_Barlow(FSI_xsec,PWIA_xsec,title='Xsec',settings=s)
RBar.calc_R()
# RBar.plot_R()
RBar.plot_R(th_nq=[35,45])

FSI_RC = {'pm_120':h_JML['pm_120_JML_FSI_radcorr_factor_proj'],
          'pm_580':h_JML['pm_580_JML_FSI_radcorr_factor_proj'],
          'pm_800':h_JML['pm_800_JML_FSI_radcorr_factor_proj'],
          'pm_900':h_JML['pm_900_JML_FSI_radcorr_factor_proj']}
PWIA_RC = {'pm_120':h_JML['pm_120_JML_PWIA_radcorr_factor_proj'],
          'pm_580':h_JML['pm_580_JML_PWIA_radcorr_factor_proj'],
          'pm_800':h_JML['pm_800_JML_PWIA_radcorr_factor_proj'],
          'pm_900':h_JML['pm_900_JML_PWIA_radcorr_factor_proj']}

RBar = RB.R_Barlow(FSI_RC,PWIA_RC,title='Radiative Corrections',settings=s,
                   ana_type='factor')
RBar.calc_R()
# RBar.plot_R()

#%% calculate R Barlow for HCOLL A, B, C

# HCOLL A: scale = 1.1
HCOLLA_xsec = {'pm_120':h_HCOLL_A['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_HCOLL_A['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_HCOLL_A['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_HCOLL_A['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_HCOLLA = RB.R_Barlow(HCOLLA_xsec,CDB_xsec,
                          title='Xsec HMS Collimator Scale = 1.1',settings=s)
RBar_HCOLLA.calc_R()
# RBar_HCOLLA.plot_R()

# HCOLL B: scale = 0.9
HCOLLB_xsec = {'pm_120':h_HCOLL_B['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_HCOLL_B['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_HCOLL_B['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_HCOLL_B['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_HCOLLB = RB.R_Barlow(HCOLLB_xsec,CDB_xsec,
                          title='Xsec HMS Collimator Scale = 0.9',settings=s)
RBar_HCOLLB.calc_R()
# RBar_HCOLLB.plot_R()

# HCOLL C: scale = 0.5
HCOLLC_xsec = {'pm_120':h_HCOLL_C['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_HCOLL_C['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_HCOLL_C['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_HCOLL_C['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_HCOLLC = RB.R_Barlow(HCOLLC_xsec,CDB_xsec,
                          title='Xsec HMS Collimator Scale = 0.5',settings=s)
RBar_HCOLLC.calc_R()
# RBar_HCOLLC.plot_R()
RBar_HCOLLC.plot_R(th_nq=[35,45])

#%% calculate R Barlow for shms calPID OFF

pcalPID_OFF_xsec = {'pm_120':h_pcalPID_OFF['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_pcalPID_OFF['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_pcalPID_OFF['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_pcalPID_OFF['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_pcalPID_OFF = RB.R_Barlow(pcalPID_OFF_xsec,CDB_xsec,
                          title='Xsec SHMS Calorimeter PID OFF',settings=s)
RBar_pcalPID_OFF.calc_R()
# RBar_pcalPID_OFF.plot_R()
RBar_pcalPID_OFF.plot_R(th_nq=[35,45])

#%% calculate R Barlow for shms calPID A, B, C

pcalPID_A_xsec = {'pm_120':h_pcalPID_A['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_pcalPID_A['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_pcalPID_A['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_pcalPID_A['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_pcalPID_A = RB.R_Barlow(pcalPID_A_xsec,CDB_xsec,
                          title='Xsec SHMS Calorimeter PID (0.8,1.2)',settings=s)
RBar_pcalPID_A.calc_R()
# RBar_pcalPID_A.plot_R()

pcalPID_B_xsec = {'pm_120':h_pcalPID_B['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_pcalPID_B['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_pcalPID_B['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_pcalPID_B['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_pcalPID_B = RB.R_Barlow(pcalPID_B_xsec,CDB_xsec,
                          title='Xsec SHMS Calorimeter PID (0.9,1.1)',settings=s)
RBar_pcalPID_B.calc_R()
# RBar_pcalPID_B.plot_R()
RBar_pcalPID_B.plot_R(th_nq=[35,45])

pcalPID_C_xsec = {'pm_120':h_pcalPID_C['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_pcalPID_C['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_pcalPID_C['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_pcalPID_C['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_pcalPID_C = RB.R_Barlow(pcalPID_C_xsec,CDB_xsec,
                          title='Xsec SHMS Calorimeter PID (0.95,1.05)',settings=s)
RBar_pcalPID_C.calc_R()
# RBar_pcalPID_C.plot_R()

#%% calculate R Barlow for Em A, B, C, D, E
# Em A: (-0.08,0.08)
EmA_xsec = {'pm_120':h_Em_A['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_Em_A['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_Em_A['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_Em_A['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_EmA = RB.R_Barlow(EmA_xsec,CDB_xsec,
                          title='Xsec Missing Energy Cut (-0.08,0.08)',settings=s)
RBar_EmA.calc_R()
# RBar_EmA.plot_R()

# Em B: (-0.07,0.07)
EmB_xsec = {'pm_120':h_Em_B['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_Em_B['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_Em_B['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_Em_B['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_EmB = RB.R_Barlow(EmB_xsec,CDB_xsec,
                          title='Xsec Missing Energy Cut (-0.07,0.07)',settings=s)
RBar_EmB.calc_R()
# RBar_EmB.plot_R()

# Em C: (-0.06,0.06)
EmC_xsec = {'pm_120':h_Em_C['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_Em_C['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_Em_C['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_Em_C['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_EmC = RB.R_Barlow(EmC_xsec,CDB_xsec,
                          title='Xsec Missing Energy Cut (-0.06,0.06)',settings=s)
RBar_EmC.calc_R()
# RBar_EmC.plot_R()

# Em D: (-0.04,0.04)
EmD_xsec = {'pm_120':h_Em_D['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_Em_D['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_Em_D['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_Em_D['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_EmD = RB.R_Barlow(EmD_xsec,CDB_xsec,
                          title='Xsec Missing Energy Cut (-0.04,0.04)',settings=s)
RBar_EmD.calc_R()
# RBar_EmD.plot_R()
RBar_EmD.plot_R(th_nq=[35,45])

# Em E: (-0.03,0.03)
EmE_xsec = {'pm_120':h_Em_E['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_Em_E['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_Em_E['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_Em_E['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_EmE = RB.R_Barlow(EmE_xsec,CDB_xsec,
                          title='Xsec Missing Energy Cut (-0.03,0.03)',settings=s)
RBar_EmE.calc_R()
# RBar_EmE.plot_R()

#%% calculate R Barlow for Ztar Difference A, B, C
# ZT A: +/- 2.5
ZTA_xsec = {'pm_120':h_ZT_A['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_ZT_A['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_ZT_A['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_ZT_A['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_ZTA = RB.R_Barlow(ZTA_xsec,CDB_xsec,
                          title='Xsec Ztar Difference Cut $\pm$2.5',settings=s)
RBar_ZTA.calc_R()
# RBar_ZTA.plot_R()

# ZT B: +/- 1.5
ZTB_xsec = {'pm_120':h_ZT_B['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_ZT_B['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_ZT_B['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_ZT_B['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_ZTB = RB.R_Barlow(ZTB_xsec,CDB_xsec,
                          title='Xsec Ztar Difference Cut $\pm$1.5',settings=s)
RBar_ZTB.calc_R()
# RBar_ZTB.plot_R()
RBar_ZTB.plot_R(th_nq=[35,45])

# ZT C: +/- 1.0
ZTC_xsec = {'pm_120':h_ZT_C['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_ZT_C['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_ZT_C['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_ZT_C['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_ZTC = RB.R_Barlow(ZTC_xsec,CDB_xsec,
                          title='Xsec Ztar Difference Cut $\pm$1.0',settings=s)
RBar_ZTC.calc_R()
# RBar_ZTC.plot_R()

#%% calculate R Barlow for epCoinTime OFF
CT_OFF_xsec = {'pm_120':h_CT_OFF['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CT_OFF['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CT_OFF['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CT_OFF['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_CT_OFF = RB.R_Barlow(CT_OFF_xsec,CDB_xsec,
                          title='Xsec e,p Coincidence Time Cut OFF',settings=s)
RBar_CT_OFF.calc_R()
# RBar_CT_OFF.plot_R()
RBar_CT_OFF.plot_R(th_nq=[35,45])

#%% calculate R Barlow for epCoinTime A, B, C
CT_A_xsec = {'pm_120':h_CT_A['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CT_A['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CT_A['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CT_A['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_CT_A = RB.R_Barlow(CT_A_xsec,CDB_xsec,
                          title='Xsec e,p Coincidence Time Cut $\pm1$',settings=s)
RBar_CT_A.calc_R()
# RBar_CT_A.plot_R()
RBar_CT_A.plot_R(th_nq=[35,45])

CT_B_xsec = {'pm_120':h_CT_B['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CT_B['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CT_B['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CT_B['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_CT_B = RB.R_Barlow(CT_B_xsec,CDB_xsec,
                          title='Xsec e,p Coincidence Time Cut $\pm0.5$',settings=s)
RBar_CT_B.calc_R()
# RBar_CT_B.plot_R()


CT_C_xsec = {'pm_120':h_CT_C['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CT_C['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CT_C['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CT_C['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_CT_C = RB.R_Barlow(CT_C_xsec,CDB_xsec,
                          title='Xsec e,p Coincidence Time Cut $\pm0.25$',settings=s)
RBar_CT_C.calc_R()
# RBar_CT_C.plot_R()

#%% calculate R Barlow for noEDTM OFF
noEDTM_OFF_xsec = {'pm_120':h_noEDTM_OFF['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_noEDTM_OFF['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_noEDTM_OFF['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_noEDTM_OFF['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_noEDTM_OFF = RB.R_Barlow(noEDTM_OFF_xsec,CDB_xsec,
                          title='Xsec No EDTM Cut OFF',settings=s)
RBar_noEDTM_OFF.calc_R()
# RBar_noEDTM_OFF.plot_R()
RBar_noEDTM_OFF.plot_R(th_nq=[35,45])

#%% calculate R Barlow for BCM1 OFF
BCM1_OFF_xsec = {'pm_120':h_BCM1_OFF['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_BCM1_OFF['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_BCM1_OFF['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_BCM1_OFF['pm_900_CDB_FSI_data_xsec_proj']}
CDB_xsec = {'pm_120':h_CDB['pm_120_CDB_FSI_data_xsec_proj'],
            'pm_580':h_CDB['pm_580_CDB_FSI_data_xsec_proj'],
            'pm_800':h_CDB['pm_800_CDB_FSI_data_xsec_proj'],
            'pm_900':h_CDB['pm_900_CDB_FSI_data_xsec_proj']}

RBar_BCM1_OFF = RB.R_Barlow(BCM1_OFF_xsec,CDB_xsec,
                          title='Xsec BCM1 Current Cut OFF',settings=s)
RBar_BCM1_OFF.calc_R()
# RBar_BCM1_OFF.plot_R()
RBar_BCM1_OFF.plot_R(th_nq=[35,45])

#%% calculate R Barlow for bin centering corrections FSI/PWIA
BC_xsec = {'pm_120':{'set_1':f_BC_pm120_1,'set_2':f_BC_pm120_2},
            'pm_580':{'set_1':f_BC_pm580_1,'set_2':f_BC_pm580_2},
            'pm_800':{'set_1':f_BC_pm800_1},
            'pm_900':{'set_1':f_BC_pm900_1}}

RBar_BC = RB.R_Barlow(BC_xsec,title='Xsec Vertex Correction FSI vs PWIA',settings=s,
                      ana_type='BC')
RBar_BC.calc_R()
# RBar_BC.plot_R()
RBar_BC.plot_R(th_nq=[35,45])
