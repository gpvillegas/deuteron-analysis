import LT.box as B
from LT.datafile import dfile
import numpy as np
import numpy.ma as ma
import sys           
import os                                                      
from sys import argv  

#This Code:
#This code assumes the systematics errors have been calculated on the Xsec and red Xsec, and reads them 
#to combine the data. When combining data, ONLY used statistical weight.  To combine the systematic errors
#from different data sets, take the normal average of the systematic errors for now (Mark Jones)
#When plotting combined reduced Xsec, plot statistical and total error bars (two error bars per data point)

# **NOTE** ONLY Reduced Cross Sections (Momentum Distributions) can be combined. Cross Sections CANNOT be combined, 
#          as they depend on spectrometer the kinematic setting, whereas in momentum distribution, that depenency has
#          been factored out by dividing by K*sig_cc1(GEp, GMp)
# ** combines (takes average) different data sets of a given kinematic setting and writes to file
# ** combines different kinematic settings

# some constants                                                                                                                                                     
dtr = np.pi/180.                                                                                                                                                                               
#GeV                                                                                                                                                                             
MP = 938.272 / 1000.                                                                                                                                                                       
MN = 939.566 / 1000.                                                                                                                                                 
MD = 1875.61 / 1000.                                                                                                                                                                           
me = 0.51099 / 1000. 

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)

    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr

MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#===========User Input (Dir. Name to store output)==================
#ipython combine_data.py Em_final40MeV
# sys_ext = sys.argv[1]   
#===================================================================

#check if directory exists, else creates it.
# if not os.path.exists(sys_ext):
#     os.makedirs(sys_ext)

output_file='yield_n_xsec/combine_data/redXsec_combined.txt'

fout = open(output_file, 'w') 

#------------------------------------------------------------                   
# header information for the output file    
header = """        
#This file contains combined Reduced Cross Sections from all data sets of the same kinematic bin. pm80 + pm580 (sets 1 and 2) +  pm750 (sets 1, 2 and 3)                                                                                                                                                                                           
#\\ xb = th_nq                                                                                                    
#\\ yb = pm                                                                        
#relative errors df/f=dsig/sig: kin_syst_tot (kinematic systematics),  norm_syst_tot (constant, norm. syst.)
# current header line: all averaged kinematics are either in GeV or deg : Ei_avg (beam energy), kf_avg(final e- momentum),  the_avg(final e- angle), pf_avg(final proton momentum),
# nu_avg(energy transfer, Ei-kf), Q2_avg(4-momentum transfer, GeV2), q_avg(|qlab|, magnitude of 3-momentum transfer), cthpq_cm_avg(c.o.mass angle between proton and qlab), cphi_avg(out-of-plane angle between proton and qlab)
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/  pm_avg[f,5]/  kin_syst_tot[f,6]/  norm_syst_tot[f,7]/  tot_syst_err[f,8]/  tot_stats_err[f,9]/  tot_err[f,10]/   red_dataXsec_avg[f,11]/   red_dataXsec_avg_err[f,12]/  red_dataXsec_avg_syst_err[f,13]/   red_dataXsec_avg_tot_err[f,14]/  red_pwiaXsec_avg[f,15]/ red_fsiXsec_avg[f,16]/ thnq_avg[f,17]/  
"""
fout.write(header)

def get_sys_fname(pm_set, data_set):
    # if(pm_set==80):
    #     fname = '../average_kinematics/%s/pm%i_fsi_norad_avgkin_systematics.txt' %(sys_ext, pm_set)
    # else:
    #     fname = '../average_kinematics/%s/pm%i_fsi_norad_avgkin_set%i_systematics.txt' %(sys_ext, pm_set, data_set)
    fname = ''
    return fname

def get_fname(pm_set, model, data_set):
    # if(pm_set==80):
    #     fname = '../bin_centering_corrections/%s/pm%i_laget_bc_corr.txt' %(sys_ext, pm_set)
    # else:
    #     fname = '../bin_centering_corrections/%s/pm%i_laget_bc_corr_set%i.txt' %(sys_ext, pm_set, data_set)
    fname = f'yield_n_xsec/bin_centering/{pm_set}_{model}_bc_corr_{data_set}.txt'

    return fname

def get_pm_avg(pm_set, model):
    #Code to get the average momentum from the avg. kin. file
    # if(pm_set==80):
    #     fname = '../average_kinematics/%s/pm%i_fsi_norad_avgkin.txt' %(sys_ext, pm_set)
    # else:
    #     fname = '../average_kinematics/%s/pm%i_fsi_norad_avgkin_set%i.txt' %(sys_ext, pm_set, data_set)
    fname = f'yield_n_xsec/average_kin/{pm_set}_{model}fsi_norad_avgkin.txt'

    kin = dfile(fname)
    pm = kin['pm']
    pm_bin = kin['yb']
    # lowEdge = pm_bin - 20.  #lower edged of Pm bin
    # upEdge = pm_bin + 20.   #upper edge of Pm bin
    return pm_bin, pm

def get_avg_kin(header_name, pm_set, model):                                                                                                                                                  
    #Code to get any header and average kinematics dara array                                                                                              
    # if(pm_set==80):
    #     fname = '../average_kinematics/%s/pm%i_fsi_norad_avgkin.txt' %(sys_ext, pm_set) 
    # else:
    #     fname = '../average_kinematics/%s/pm%i_fsi_norad_avgkin_set%i.txt' %(sys_ext, pm_set, data_set)
    fname = f'yield_n_xsec/average_kin/{pm_set}_{model}fsi_norad_avgkin.txt'
    kin = dfile(fname)
    data = kin[header_name]                                                                                                                                                          
    return data 

def get_sig_red(header_name, pm_set, model, data_set):
    #Code to get any header and data array from bin-center corrected data files

    #Get file containing red. Xsec and its stats. error
    f = dfile(get_fname(pm_set, model, data_set))
    
    #Get data array with desired header name
    data = f[header_name]
    return data

def get_kin_syst(header_name, pm_set, data_set):
    #Code to get any header and data array from the kin. systematic data files

    #Get file containing kin. systematics
    f = dfile(get_sys_fname(pm_set, data_set))
    
    #Get data array with desired header name (kin. syst are in %)
    data = f[header_name] / 100.   #coonvert to fractional relative error
    return data

def get_norm_syst(header_name):
    #Code to get any header and data array from the normalization. systematic data file

    #Get file containing kin. systematics
    #f = dfile('../average_kinematics/Em_final40MeV/normalization_systematics.txt')
    # f = dfile('../average_kinematics/%s/normalization_systematics_summary.txt' % (sys_ext))

    #Get data array with desired header name (norm. syst are in fractional rel. error)
    # data = f[header_name] 
    data = ''
    return data

#Target wall and spectrometer aceptance contributions
#relative error on cross sections will be added quadratically to normalization (October, 19, 2020)
tgt_wall_err = 0.029   #target wall contributes at most 2.9 % to the yield (integrated over all thnq,)
spec_acc_err = 0.014   #spectrometer acceptance contributes 1.4 % to the total cross section (study by M.K.Jones)


#Get bin information from each file
f1  = dfile(get_fname('pm_120', 'jml', 'set_1'))
i_b = f1['i_b']
i_x = f1['i_x']
i_y = f1['i_y']
xb  = f1['xb']   #th_nq bin
yb  = f1['yb']   #pm_bin

#Get Average Pmiss
pm1_b, pm120 = get_pm_avg('pm_120', 'jml')
pm2_b, pm580 = get_pm_avg('pm_580', 'jml')
pm3_b, pm800 = get_pm_avg('pm_800', 'jml')
pm4_b, pm900 = get_pm_avg('pm_900', 'jml')

#Get Average initial beam energy [MeV]           Get avg final e- energy [MeV]                     Get final e- angle [deg]                
Ei_120 = get_avg_kin('Ei', 'pm_120', 'jml') ;    kf_120 = get_avg_kin('kf', 'pm_120', 'jml') ;     the_120 = get_avg_kin('th_e', 'pm_120', 'jml')
Ei_580 = get_avg_kin('Ei', 'pm_580', 'jml') ;    kf_580 = get_avg_kin('kf', 'pm_580', 'jml') ;     the_580 = get_avg_kin('th_e', 'pm_580', 'jml')
Ei_800 = get_avg_kin('Ei', 'pm_800', 'jml') ;    kf_800 = get_avg_kin('kf', 'pm_800', 'jml') ;     the_800 = get_avg_kin('th_e', 'pm_800', 'jml')
Ei_900 = get_avg_kin('Ei', 'pm_900', 'jml') ;    kf_900 = get_avg_kin('kf', 'pm_900', 'jml') ;     the_900 = get_avg_kin('th_e', 'pm_900', 'jml')
#omega, Q2_calc, q_lab, theta_pq_calc, th_pq_cm, cos_phi
#Get final proton momentum [MeV]                 #Get omega average [MeV]                          #Get averaged Q2 [MeV2]
pf_120 = get_avg_kin('pf', 'pm_120', 'jml') ;    nu_120 = get_avg_kin('omega', 'pm_120', 'jml') ;  Q2_120 = get_avg_kin('Q2_calc', 'pm_120', 'jml')     ;
pf_580 = get_avg_kin('pf', 'pm_580', 'jml') ;    nu_580 = get_avg_kin('omega', 'pm_580', 'jml') ;  Q2_580 = get_avg_kin('Q2_calc', 'pm_580', 'jml')     ;
pf_800 = get_avg_kin('pf', 'pm_800', 'jml') ;    nu_800 = get_avg_kin('omega', 'pm_800', 'jml') ;  Q2_800 = get_avg_kin('Q2_calc', 'pm_800', 'jml')     ;
pf_900 = get_avg_kin('pf', 'pm_900', 'jml') ;    nu_900 = get_avg_kin('omega', 'pm_900', 'jml') ;  Q2_900 = get_avg_kin('Q2_calc', 'pm_900', 'jml')     ;

#Get final |q_lab| MeV                           #Get th_pq_cm                                                #Get cos(phi_pq)
q_120 = get_avg_kin('q_lab', 'pm_120', 'jml') ;  thpq_cm_120 = get_avg_kin('th_pq_cm', 'pm_120', 'jml') ;     cphi_120 = get_avg_kin('cos_phi', 'pm_120', 'jml')     ;
q_580 = get_avg_kin('q_lab', 'pm_580', 'jml') ;  thpq_cm_580 = get_avg_kin('th_pq_cm', 'pm_580', 'jml') ;     cphi_580 = get_avg_kin('cos_phi', 'pm_580', 'jml')     ;
q_800 = get_avg_kin('q_lab', 'pm_800', 'jml') ;  thpq_cm_800 = get_avg_kin('th_pq_cm', 'pm_800', 'jml') ;     cphi_800 = get_avg_kin('cos_phi', 'pm_800', 'jml')     ;
q_900 = get_avg_kin('q_lab', 'pm_900', 'jml') ;  thpq_cm_900 = get_avg_kin('th_pq_cm', 'pm_900', 'jml') ;     cphi_900 = get_avg_kin('cos_phi', 'pm_900', 'jml')     ;

#Get thnq [deg]
thnq_120 = get_avg_kin('th_nq_mc','pm_120','jml')
thnq_580 = get_avg_kin('th_nq_mc','pm_580','jml')
thnq_800 = get_avg_kin('th_nq_mc','pm_800','jml')
thnq_900 = get_avg_kin('th_nq_mc','pm_900','jml')

#Get reduced Xsec
sig_red_120_set1 = get_sig_red('red_dataXsec', 'pm_120', 'jml', 'set_1')
sig_red_120_set2 = get_sig_red('red_dataXsec', 'pm_120', 'jml', 'set_2')
sig_red_580_set1 = get_sig_red('red_dataXsec', 'pm_580', 'jml', 'set_1')
sig_red_580_set2 = get_sig_red('red_dataXsec', 'pm_580', 'jml', 'set_2')
sig_red_800_set1 = get_sig_red('red_dataXsec', 'pm_800', 'jml', 'set_1')
sig_red_900_set1 = get_sig_red('red_dataXsec', 'pm_900', 'jml', 'set_1')

#Get reduced Xsec stats. err
sig_red_120_set1_err = get_sig_red('red_dataXsec_err', 'pm_120', 'jml', 'set_1')
sig_red_120_set2_err = get_sig_red('red_dataXsec_err', 'pm_120', 'jml', 'set_2')
sig_red_580_set1_err = get_sig_red('red_dataXsec_err', 'pm_580', 'jml', 'set_1')
sig_red_580_set2_err = get_sig_red('red_dataXsec_err', 'pm_580', 'jml', 'set_2')
sig_red_800_set1_err = get_sig_red('red_dataXsec_err', 'pm_800', 'jml', 'set_1')
sig_red_900_set1_err = get_sig_red('red_dataXsec_err', 'pm_900', 'jml', 'set_1')

#Get theoretical red. Xsec
red_pwiaXsec_120_set1 = get_sig_red('red_pwiaXsec', 'pm_120', 'jml', 'set_1')
red_pwiaXsec_120_set2 = get_sig_red('red_pwiaXsec', 'pm_120', 'jml', 'set_2')
red_pwiaXsec_580_set1 = get_sig_red('red_pwiaXsec', 'pm_580', 'jml', 'set_1')
red_pwiaXsec_580_set2 = get_sig_red('red_pwiaXsec', 'pm_580', 'jml', 'set_2')
red_pwiaXsec_800_set1 = get_sig_red('red_pwiaXsec', 'pm_800', 'jml', 'set_1')
red_pwiaXsec_900_set1 = get_sig_red('red_pwiaXsec', 'pm_800', 'jml', 'set_1')

red_fsiXsec_120_set1 = get_sig_red('red_fsiXsec', 'pm_120', 'jml', 'set_1')
red_fsiXsec_120_set2 = get_sig_red('red_fsiXsec', 'pm_120', 'jml', 'set_2')
red_fsiXsec_580_set1 = get_sig_red('red_fsiXsec', 'pm_580', 'jml', 'set_1')
red_fsiXsec_580_set2 = get_sig_red('red_fsiXsec', 'pm_580', 'jml', 'set_2')
red_fsiXsec_800_set1 = get_sig_red('red_fsiXsec', 'pm_800', 'jml', 'set_1')
red_fsiXsec_900_set1 = get_sig_red('red_fsiXsec', 'pm_900', 'jml', 'set_1')

#Get Kinematic Systematics Relative error
# kin_syst_120 = get_kin_syst('sig_kin_tot', 80, 1)
# kin_syst_580_set1 = get_kin_syst('sig_kin_tot', 580, 1)
# kin_syst_580_set2 = get_kin_syst('sig_kin_tot', 580, 2)
# kin_syst_750_set1 = get_kin_syst('sig_kin_tot', 750, 1)
# kin_syst_750_set2 = get_kin_syst('sig_kin_tot', 750, 2)
# kin_syst_750_set3 = get_kin_syst('sig_kin_tot', 750, 3)


# pm_set = get_norm_syst('pm')
# data_set = get_norm_syst('set')

#Get norm. syst. errors that are the same for all data sets (these should be added as an overall norm. factor, and not in quadrature, as they affect all data sets the same)
# pAbs_syst = get_norm_syst('pAbs_syst')[0]   #proton absorption systematic relative error (dsig / sig), (0.4951 %)
# tLT_syst = get_norm_syst('tLT_syst')[0]   #total live time systematics (3 %)
# Qtot_syst = get_norm_syst('Qtot_syst')[0]   #total charge systematics (2 %)
# const_norm_syst = np.sqrt(pAbs_syst**2 + tLT_syst**2 + Qtot_syst**2 + tgt_wall_err**2 + spec_acc_err**2)   #constant norm. systematics to be added as a single value later

# print('pAbs_syst = ',pAbs_syst)
# print('tLT_syst = ', tLT_syst)
# print('Qtot_syst = ', Qtot_syst)
# print('const_norm_syst = ', const_norm_syst)

#Get norm. syst. errors that change per data set 
# htrk_eff_syst = get_norm_syst('htrk_eff_syst')
# etrk_eff_syst = get_norm_syst('etrk_eff_syst')
# tgtBoil_syst = get_norm_syst('tgtBoil_syst')
# print('htrk_eff_syst = ',htrk_eff_syst)
# print('etrk_eff_syst = ',etrk_eff_syst)
# print('tgtBoil_syst = ',tgtBoil_syst)

#norm_syst_tot = get_norm_syst('total_norm_syst')  #new syst. array (systematcis for each separate set, added in quadrature),  norm_syst_tot[i] -> i=0: 80 MeV setting, i=1: 580 (set1) setting, . . . 

# norm_syst_80 = np.sqrt(htrk_eff_syst[0]**2 + etrk_eff_syst[0]**2 + tgtBoil_syst[0]**2)
# norm_syst_580_set1 = np.sqrt(htrk_eff_syst[1]**2 + etrk_eff_syst[1]**2 + tgtBoil_syst[1]**2)
# norm_syst_580_set2 = np.sqrt(htrk_eff_syst[2]**2 + etrk_eff_syst[2]**2 + tgtBoil_syst[2]**2)
# norm_syst_750_set1 = np.sqrt(htrk_eff_syst[3]**2 + etrk_eff_syst[3]**2 + tgtBoil_syst[3]**2)
# norm_syst_750_set2 = np.sqrt(htrk_eff_syst[4]**2 + etrk_eff_syst[4]**2 + tgtBoil_syst[4]**2)
# norm_syst_750_set3 = np.sqrt(htrk_eff_syst[5]**2 + etrk_eff_syst[5]**2 + tgtBoil_syst[5]**2)
# print('norm_syst_80 = ', norm_syst_80 )
# print('norm_syst_580_set1 = ', norm_syst_580_set1 ) 
# print('norm_syst_580_set2 = ', norm_syst_580_set2 )       
# print('norm_syst_750_set1 = ', norm_syst_750_set1 )                                                                                                                                       
# print('norm_syst_750_set2 = ', norm_syst_750_set2 ) 
# print('norm_syst_750_set3 = ', norm_syst_750_set3 )                                                                                                                                                     

#Loop over all 2d bins
for ib in range(len(i_b)):

    
    #----------Calculate the Data Reduced Xsec Weighted Average for a given (Pm, th_nq) bin of a given kin. set------------
    red_dataXsec_arr = np.array([sig_red_120_set1[ib], sig_red_120_set2[ib], 
                                 sig_red_580_set1[ib], sig_red_580_set2[ib], 
                                 sig_red_800_set1[ib], sig_red_900_set1[ib]])
    red_dataXsec_m = ma.masked_values(red_dataXsec_arr, -1.)  #maske invalid values (-1)
    
    #Define the weights (Use ONLY statistical)
    red_dataXsec_weights = np.array([1./sig_red_120_set1_err[ib]**2, 
                                     1./sig_red_120_set2_err[ib]**2, 
                                     1./sig_red_580_set1_err[ib]**2,
                                     1./sig_red_580_set2_err[ib]**2, 
                                     1./sig_red_800_set1_err[ib]**2, 
                                     1./sig_red_900_set1_err[ib]**2])
    red_dataXsec_weights_m = ma.masked_values(red_dataXsec_weights, 1.)


    #=======DATA REDUCED CROSS SECTION WEIGHTED AVERAGE==========
    red_dataXsec_avg = ma.average(red_dataXsec_m, weights=red_dataXsec_weights_m)
    red_dataXsec_avg_err = 1. / np.sqrt(np.sum(red_dataXsec_weights_m))    #combined data sets statistical uncertainty per (Pm, thnq) bin

    #Construct array of kinematic uncertainty for each bin
    # kin_syst_arr = np.array([kin_syst_80[ib], kin_syst_580_set1[ib], kin_syst_580_set2[ib], kin_syst_750_set1[ib], kin_syst_750_set2[ib], kin_syst_750_set3[ib]])
    kin_syst_arr = np.zeros(6)
    
    #mask elements corresponding to masked dataXsec (We do not want to add in quadrature elements for which there is no Xsec)
    kin_syst_arr_m = ma.masked_array(kin_syst_arr, ma.getmask(red_dataXsec_m))

    #Add masked kinematic systematics array in quadrature
    kin_syst_tot2 = np.sum(kin_syst_arr_m**2)
    kin_syst_tot = np.sqrt( kin_syst_tot2 )


    #March 6, 2020
    #--Construct array of normalization uncertainties--
    # norm_syst_arr = np.array([norm_syst_80, norm_syst_580_set1, norm_syst_580_set2, norm_syst_750_set1, norm_syst_750_set2, norm_syst_750_set3])
    norm_syst_arr = np.zeros(6)

    #print('norm_syst_arr = ', norm_syst_arr)
    #print('const_norm_syst = ', const_norm_syst)
    #mask elements corresponding to masked dataXsec (We do not want to add in quadrature elements for which there is no Xsec)
    norm_syst_arr_m = ma.masked_array(norm_syst_arr, ma.getmask(red_dataXsec_m))
    #Add masked norm systematics array in quadrature
    # norm_syst_final2 = np.sum(norm_syst_arr_m**2 + const_norm_syst**2)
    # norm_syst_final2 = np.sum(norm_syst_arr_m**2)
    # norm_syst_final = np.sqrt( norm_syst_final2 )

    # norm_syst_final2_corr = np.sum(norm_syst_arr_m**2) + const_norm_syst**2
    norm_syst_final2_corr = np.sum(norm_syst_arr_m**2) 
    norm_syst_final_corr = np.sqrt( norm_syst_final2_corr )        

# =============================================================================
#     if norm_syst_final > 0.08:
#         print('------WRONG WAY-------------')
#         print('norm_syst_arr_m = ',norm_syst_arr_m)
#         print('const_norm_syst = ',const_norm_syst)
#         print('norm_syst_arr_m**2 + const_norm_syst**2 = ', norm_syst_arr_m**2 + const_norm_syst**2)
#         print('norm_syst_final2 = ', norm_syst_final2 )
#         print('norm_syst_final = ', norm_syst_final )     
# 
#         print('------RIGHT WAY-------------')                                                                                                                                                       
#         print('norm_syst_arr_m = ',norm_syst_arr_m)                                                                                                                                                                      
#         print('const_norm_syst = ',const_norm_syst)                                                                                                                                                                      
#         print('norm_syst_arr_m**2 = ', norm_syst_arr_m**2)                                                                                                                    
#         print('np.sum(norm_syst_arr_m**2) + const_norm_syst**2 = ', norm_syst_final2_corr )                                                                                                                          
#         print('norm_syst_final_correct = ', norm_syst_final_corr )  
# =============================================================================

    #print('norm_syst_final=',norm_syst_final)
    #Add total systematic error in quadrature
    tot_syst_err = np.sqrt(kin_syst_tot**2 + norm_syst_final_corr**2)

    #Define relative statistical relative error
    tot_stats_err = red_dataXsec_avg_err/red_dataXsec_avg 
    
    #Add statistical and systematics relative error in quadrature
    tot_err = np.sqrt(tot_stats_err**2 + tot_syst_err**2)
    
    #===Calculate Absolute Error on the Reduced Cross Section===
    red_dataXsec_avg_syst_err = red_dataXsec_avg*tot_syst_err
    red_dataXsec_avg_tot_err = red_dataXsec_avg*tot_err


    #================
    # THEORY AVERAGE
    #================
    
    red_pwiaXsec_arr = np.array([red_pwiaXsec_120_set1[ib], 
                                 red_pwiaXsec_120_set2[ib], 
                                 red_pwiaXsec_580_set1[ib], 
                                 red_pwiaXsec_580_set2[ib], 
                                 red_pwiaXsec_800_set1[ib], 
                                 red_pwiaXsec_900_set1[ib]])
    red_pwiaXsec_arr_m = np.ma.masked_values(red_pwiaXsec_arr, -1.)
    red_pwiaXsec_avg = np.ma.average(red_pwiaXsec_arr_m)
   
    red_fsiXsec_arr = np.array([red_fsiXsec_120_set1[ib], 
                                 red_fsiXsec_120_set2[ib], 
                                 red_fsiXsec_580_set1[ib], 
                                 red_fsiXsec_580_set2[ib], 
                                 red_fsiXsec_800_set1[ib], 
                                 red_fsiXsec_900_set1[ib]])
    red_fsiXsec_arr_m = np.ma.masked_values(red_fsiXsec_arr, -1.)
    red_fsiXsec_avg = np.ma.average(red_fsiXsec_arr_m)

    #==================================
    # Get the Average Missing Momentum
    #==================================
    # make mask for simc arrays
    pm120_m = (ma.getmask(red_dataXsec_m)[0] == False) or (ma.getmask(red_dataXsec_m)[1] == False)
    pm580_m = (ma.getmask(red_dataXsec_m)[2] == False) or (ma.getmask(red_dataXsec_m)[3] == False)
    pm800_m = ma.getmask(red_dataXsec_m)[4]
    pm900_m = ma.getmask(red_dataXsec_m)[5]
    simc_mask = np.array([~pm120_m,~pm580_m,pm800_m,pm900_m])
    
    pm_arr = np.array([pm120[ib], pm580[ib], pm800[ib], pm900[ib]])

    #pm_arr_m = np.ma.masked_values(pm_arr, 0.)
    #pm_avg = np.ma.average(pm_arr_m) / 1000.
     
    pm_arr_m = ma.masked_array(pm_arr, simc_mask)     #mask elements for which there is no cross section
    pm_avg = np.ma.average(pm_arr_m[pm_arr!=0.]) / 1000.  # average :: Convert to GeV  this is the true averaged recoil momentum
    
    # if xb[ib] == 35:
    #     print()
    #     print(ib)
    #     print('pm_bin=',yb[ib])
    #     print('pm_arr=',pm_arr)
    #     print('data_mask=',ma.getmask(red_dataXsec_m))
    #     print('simc_mask=',simc_mask)
    #     print('pm_mask=',pm_arr_m)
    #     print('pm_avg=',pm_avg)
    
    #Get average final proton momentum
    pf_arr = np.array([pf_120[ib], pf_580[ib], pf_800[ib], pf_900[ib]])
    pf_arr_m = np.ma.masked_values(pf_arr, simc_mask)
    pf_avg = np.ma.average(pf_arr_m[pf_arr_m!=0]) / 1000.   #convert to GeV
    #if(xb[ib]==25):
    #    print('pm_bin=',yb[ib])                                                                                                                                     
    #    print('pf_arr=',pf_arr)                                                                                                                                           
    #    print('pf_mask=',pf_arr_m)                                                                                                                                          
    #    print('pf_avg>>>=',pf_avg)  


    #=======Get Other Averaged Kinematics========
    #Get average neutron recoil angle
    thnq_arr = np.array([thnq_120[ib], thnq_580[ib], thnq_800[ib], thnq_900[ib]])
    thnq_arr_m = ma.masked_array(thnq_arr, simc_mask)
    thnq_avg = np.ma.average(thnq_arr_m[thnq_arr_m!=0])
    if xb[ib] == 35:
        print()
        print(ib)
        print('pm_bin=',yb[ib])
        print('thnq_arr=',thnq_arr)
        print('data_mask=',ma.getmask(red_dataXsec_m))
        print('simc_mask=',simc_mask)
        print('thnq_mask=',thnq_arr_m)
        print('thnq_avg=',thnq_avg)
    
    Ei_arr = np.array([Ei_120[ib], Ei_580[ib], Ei_800[ib], Ei_900[ib]])
    Ei_arr_m = ma.masked_array(Ei_arr, simc_mask)
    Ei_avg = np.average(Ei_arr_m[Ei_arr_m!=0]) / 1000. #convert to GeV
    #if(xb[ib]==25):                                                                                                                                         
    #    print('pm_bin=',yb[ib])                                                                                                                                        
    #    print('Ei_arr=',Ei_arr)                                                                                                                                       
    #    print('Ei_mask=',Ei_arr_m)                                                                                                                                      
    #    print('Ei_avg>>>=',Ei_avg)  
    
    kf_arr = np.array([kf_120[ib], kf_580[ib], kf_800[ib], kf_900[ib]])  
    kf_arr_m = ma.masked_array(kf_arr, simc_mask)                                                                
    kf_avg = np.average(kf_arr_m[kf_arr_m!=0]) / 1000. #convert to GeV                                                                                             
    #if(xb[ib]==25):                                                                                                                     
    #    print('pm_bin=',yb[ib])                                                                              
    #    print('kf_arr=',kf_arr)                                                                                      
    #    print('kf_mask=',kf_arr_m)                                                                                                                 
    #    print('kf_avg>>>=',kf_avg)

    the_arr = np.array([the_120[ib], the_580[ib], the_800[ib], the_900[ib]])
    the_arr_m = ma.masked_array(the_arr, simc_mask) 
    the_avg = np.average(the_arr_m[the_arr_m!=0])   #deg
    #if(xb[ib]==25):                                                                                                                                            
    #    print('pm_bin=',yb[ib])                                                                                                        
    #    print('the_arr=',the_arr)                                                                                                      
    #    print('the_mask=',the_arr_m)                                                                                                                          
    #    print('the_avg>>>=',the_avg) 
 
    nu_arr = np.array([nu_120[ib], nu_580[ib], nu_800[ib], nu_900[ib]])
    nu_arr_m = ma.masked_array(nu_arr, simc_mask) 
    nu_avg = np.average(nu_arr_m) / 1000.   
    #if(xb[ib]==15):                                                                                                                              
    #    print('pm_bin=',yb[ib])                                                                              
    #    print('nu_arr=',nu_arr)                                                                                                      
    #    print('nu_mask=',nu_arr_m)                                                                                                         
    #    print('nu_avg>>>=',nu_avg)   

    Q2_arr = np.array([Q2_120[ib], Q2_580[ib], Q2_800[ib], Q2_900[ib]])                                                          
    Q2_arr_m = ma.masked_array(Q2_arr, simc_mask)                                                                                                      
    Q2_avg = np.average(Q2_arr_m) / 1E6                                                                                                               
    #if(xb[ib]==15):                                                                                                                                               
    #    print('pm_bin=',yb[ib])                                                                                                                         
    #    print('Q2_arr=',Q2_arr)                                                                                                                                    
    #    print('Q2_mask=',Q2_arr_m)                                                                                                                               
    #    print('Q2_avg>>>=',Q2_avg)

    q_arr = np.array([q_120[ib], q_580[ib], q_800[ib], q_900[ib]])                                               
    q_arr_m = ma.masked_array(q_arr, simc_mask)                                                                                                                        
    q_avg = np.average(q_arr_m) / 1000.                                                                                                                                             
    #if(xb[ib]==35):                                                                                                                                                  
    #    print('pm_bin=',yb[ib])                                                                                                                                             
    #    print('q_arr=',q_arr)                                                                                                                        
    #    print('q_mask=',q_arr_m)                                                                                                                                
    #    print('q_avg>>>=',q_avg) 

    thpq_cm_arr = np.array([thpq_cm_120[ib], thpq_cm_580[ib], thpq_cm_800[ib], thpq_cm_900[ib]])                                                    
    thpq_cm_arr_m = ma.masked_array(thpq_cm_arr, simc_mask)                                                                                         
    thpq_cm_avg = np.average(thpq_cm_arr_m[thpq_cm_arr_m!=180.])                                                                                                                              
    cthpq_cm_avg = np.cos(thpq_cm_avg*dtr)  #calculate cos(thpq_cm)
    #if(xb[ib]==25):                                                                                                                                                           
    #    print('pm_bin=',yb[ib])                                                                                                                                                     
    #    print('thpq_cm_arr=',thpq_cm_arr)                                                                                                                                                        
    #    print('thpq_cm_mask=',thpq_cm_arr_m)                                                                                                                                                   
    #    print('thpq_cm_avg>>>=',thpq_cm_avg)  

    cphi_arr = np.array([cphi_120[ib], cphi_580[ib], cphi_800[ib], cphi_900[ib]])                
    cphi_arr_m = ma.masked_array(cphi_arr, simc_mask)                                                                              
    cphi_avg = np.average(cphi_arr_m[cphi_arr_m!=0.])                                                                                                         
    phi_avg = np.arccos(cphi_avg)*180./np.pi    #calculate phi = arccos(cos(phi)) covert to deg                                                                                                               
    #if(xb[ib]==25):                                                                                                                                                                        
    #    print('pm_bin=',yb[ib])                                                                                                                                                            
    #    print('cphi_arr=',cphi_arr)                                                                                                                                              
    #    print('cphi_mask=',cphi_arr_m)                                                                                                                                                    
    #    print('cphi_avg>>>=',cphi_avg)  
    #    print('phi_avg>>>=',phi_avg)                                                                                                                                                            
    
    l="%i  %i  %i  %f   %f   %f   %f  %f   %f   %f   %f   %.12e  %.12e  %.12e  %.12e  %.12e  %.12e %f\n" % (i_b[ib], i_x[ib], i_y[ib], xb[ib], yb[ib], pm_avg, kin_syst_tot, norm_syst_final_corr, tot_syst_err, tot_stats_err, tot_err, red_dataXsec_avg, red_dataXsec_avg_err, red_dataXsec_avg_syst_err, red_dataXsec_avg_tot_err, red_pwiaXsec_avg, red_fsiXsec_avg, thnq_avg)
    fout.write(l)
    
    

fout.close()


