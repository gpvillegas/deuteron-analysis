import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# import numpy.ma as ma
# import sys                                     
import os                                                                                                       
# from sys import argv  
import matplotlib
from matplotlib import rc
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, mark_inset)
import pandas as pd


# matplotlib.use('Agg')

#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
# font = {'family' : 'sans-serif',
#         'weight' : 'normal',
#         'size'   : 12
# }
plt.rc('font', **font)

#Set font
csfont = {'fontname':'Times New Roman'}

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#Conversion factor:  1 fm = 1/ (197 MeV),   
#The reduced cross section is in MeV^-3 units
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#User Input (Dir. Name to store output)
# sys_ext = sys.argv[1]   

dir_name = 'yield_n_xsec/combine_data/plots_prl'  
# dir_name = sys_ext+"_plots"
# dir_name_misc = dir_name + "/misc_plots"        #directory to store miscellaneous plots (trigger rates, live time, etc. . .) #obtained from the report files
print(dir_name)
#check if directory exists, else creates it.
if not os.path.exists(dir_name):
    os.makedirs(dir_name)
# if not os.path.exists(dir_name_misc):
#     os.makedirs(dir_name_misc)
     

#Get Reduced Xsec Data File
fname = 'yield_n_xsec/combine_data/redXsec_combined.txt'
f = B.get_file(fname)

#Get Bin Information (Same info for all files)                                                                                                  
i_b = B.get_data(f, 'i_b')    #2D bin number                                                                    
i_x = B.get_data(f, 'i_x')    #x (th_nq) bin number                                                                                    
i_y = B.get_data(f, 'i_y')    #y (pmiss) bin number                                        
thnq = B.get_data(f, 'xb')      #th_nq value at bin center                                                          
pm =  B.get_data(f, 'yb')      #pmiss value at bin center   
# pm_avg = pd.Series(B.get_data(f, 'pm_avg'))
pm_avg = B.get_data(f, 'pm_avg')
#pm_int = np.linspace(0.01, 1.1, 100)#B.get_data(f, 'pm_int')

# pm_avg = pm_avg.interpolate(method='cubicspline',limit=3,limit_area='inside')

#Get Combined Final Red Xsec for all kinematics
red_dataXsec_avg     = B.get_data(f,'red_dataXsec_avg')
red_dataXsec_avg_err = B.get_data(f,'red_dataXsec_avg_err')
red_dataXsec_avg_syst_err = B.get_data(f,'red_dataXsec_avg_syst_err')
red_dataXsec_avg_tot_err = B.get_data(f,'red_dataXsec_avg_tot_err')
red_pwiaXsec_avg     = pd.Series(B.get_data(f,'red_pwiaXsec_avg'))
red_fsiXsec_avg      = pd.Series(B.get_data(f,'red_fsiXsec_avg'))

# red_pwiaXsec_avg     = red_pwiaXsec_avg.interpolate(method='cubicspline',limit=3,limit_area='inside')
# red_fsiXsec_avg      = red_fsiXsec_avg.interpolate(method='cubicspline',limit=3,limit_area='inside')

#Get total relative errors / (Pm, thnq) bin for plotting (probably also write to table)
kin_syst_tot = dfile(fname)['kin_syst_tot']    #kinematic systematics
norm_syst_tot = dfile(fname)['norm_syst_tot']  #normalization systematics
tot_syst_err = dfile(fname)['tot_syst_err']    #total systematics
tot_stats_err = dfile(fname)['tot_stats_err']  #total statistical
tot_err = dfile(fname)['tot_err']              #overall error

def read_halla_data(thnq=0):

    #The code read data files containing momentum distributions from Hall A

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    #fminv2GeV = 1./(1./0.1973)
    ifm2GeV = 0.1973  #inverse fermi to GeV

    thnq_f = '%.1f' % (thnq)

    #fname = '../HallA_data/theta_%s_th_corr.data' % (thnq_f)  #older version Hall A data
    fname = 'yield_n_xsec/HallA_data/pm_distribution_results_q3_1_%s.data' % (thnq_f)  #Hall A data


    kin = dfile(fname)
    p_miss = kin['p_miss_av'] * ifm2GeV    #Pm here is not the central bin value, but an average over that bin
    red_dataXsec = kin['rho']             #units fm^3
    red_dataXsec_err = kin['delta_rho']

    return p_miss, red_dataXsec, red_dataXsec_err

def read_hallc17_data(thnq=0):

    #The code read data files containing momentum distributions from Hall C 2017

    fname = 'yield_n_xsec/HallC17_data/redXsec_HallC_thnq%i_deg.txt' % (thnq)

    kin = dfile(fname)
    p_miss = kin['pm_avg']   #Pm here is not the central bin value, but an average over that bin
    red_dataXsec = kin['data_redXsec']             #units fm^3
    red_dataXsec_err = kin['data_redXsec_tot_err']

    return p_miss, red_dataXsec, red_dataXsec_err

def read_misak_pdist(model=""):

    #The code read data files containing misak's calculation of theory momentum distirbutions using Hall C kinematics

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    #fminv2GeV = 1./(1./0.1973)
    GeV2fm = 0.1973**3  #inverse fermi^3 to GeV^3

    if(model=="Paris"):
        fname = 'yield_n_xsec/HallA_data/mom_dist_paris_mevc.data'  
    elif(model=="AV18"):
        fname = 'yield_n_xsec/HallA_data/mom_dist_v18_mevc.data'  
    elif(model=="CD-Bonn"):
        fname = 'yield_n_xsec/HallA_data/mom_dist_cdbonn_mevc.data'  


    kin = dfile(fname)
    p_miss = kin['pm'] / 1000.            #convert momentum from  MeV to GeV
    red_dataXsec = kin['rho'] / GeV2fm           #Convert GeV^-3 to fm^3 units 

    #iNTERPOLATE
    f_red_misak = interp1d(p_miss, red_dataXsec, fill_value='extrapolate')  

    return p_miss, f_red_misak

def read_theoretical_models(theory="", model="", thnq=0):

    #This code read the averaged red.Xsec and returns arrays in pm and reduced Xsec
    #theory: V18, CD-Bonn    model: PWIA, FSI
    thnq_f = "%.2f" %(thnq)

    fname = '/home/gvill/deuteron/DEUTERON_ANALYSIS/ANALYSIS_SCRIPTS/MAIN_ANALYSIS/Deep_CrossSections/theoretical_models'+\
    '/updated_%s_models/theoryXsec_%s%s_thnq%s_combined.data' % (model, theory, model, thnq_f)

    # fname = '../theoretical_models/updated_%s_models/theoryXsec_%s%s_thnq%s_combined.data' % (model, theory, model, thnq_f)
    kin = dfile(fname)
    
    pm_bin = np.array(kin['pm_avg'])
    

    #print(pm_bin)
    if(model=="PWIA" and theory=="V18"):
        red_pwiaXsec_V18 = kin['red_pwiaXsec_theory']
        return pm_bin, red_pwiaXsec_V18 
    if(model=="PWIA" and theory=="CD-Bonn"):                                                                                                                          
        red_pwiaXsec_CD_Bonn = kin['red_pwiaXsec_theory']                                                                          
        return pm_bin, red_pwiaXsec_CD_Bonn
    if(model=="FSI" and theory=="V18"):                                                 
        red_fsiXsec_V18 = kin['red_fsiXsec_theory']                                               
        return pm_bin, red_fsiXsec_V18 
    if(model=="FSI" and theory=="CD-Bonn"):                                                                                            
        red_fsiXsec_CD_Bonn = kin['red_fsiXsec_theory'] 
        return pm_bin, red_fsiXsec_CD_Bonn


def make_prl_plots():
    #Make PRL paper plots

    #Read This Experiment (Hall C) Data, and require better than 50% statistics
    red_dataXsec_avg_masked = np.ma.array(red_dataXsec_avg, mask=(red_dataXsec_avg_err>0.5*red_dataXsec_avg))
    red_dataXsec_avg_masked = np.ma.filled(red_dataXsec_avg_masked.astype(float), np.nan)

    #Read Hall A data (reduced Xsec already in fm^3 units, and Precoil in GeV)
    pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35 = read_halla_data(35)
    pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45 = read_halla_data(45)  
    pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75 = read_halla_data(75)  
    
    #Read Hall C 2017 data (reduced Xsec already in fm^3 units, and Precoil in GeV)
    pm_hc35, red_dataXsec_hc35, red_dataXsec_err_hc35 = read_hallc17_data(35)
    pm_hc45, red_dataXsec_hc45, red_dataXsec_err_hc45 = read_hallc17_data(45)  
    pm_hc75, red_dataXsec_hc75, red_dataXsec_err_hc75 = read_hallc17_data(75) 

    #Read Models 35, 45 75 (for PRL plot)
    #Read Other Theoretical Models (V18, CD-BONN)

    pm_avg1, red_pwiaXsec_V18_35 = read_theoretical_models("V18", "PWIA", 35)   
    pm_avg2, red_fsiXsec_V18_35 = read_theoretical_models("V18", "FSI", 35)      
    pm_avg3, red_pwiaXsec_CD_35 = read_theoretical_models("CD-Bonn", "PWIA", 35)                                 
    pm_avg4, red_fsiXsec_CD_35 = read_theoretical_models("CD-Bonn", "FSI", 35)   
    
    pm_avg5, red_pwiaXsec_V18_45 = read_theoretical_models("V18", "PWIA", 45)   
    pm_avg6, red_fsiXsec_V18_45 = read_theoretical_models("V18", "FSI", 45)      
    pm_avg7, red_pwiaXsec_CD_45 = read_theoretical_models("CD-Bonn", "PWIA", 45)                                 
    pm_avg8, red_fsiXsec_CD_45 = read_theoretical_models("CD-Bonn", "FSI", 45)

    pm_avg9, red_pwiaXsec_V18_75 = read_theoretical_models("V18", "PWIA", 75)   
    pm_avg10, red_fsiXsec_V18_75 = read_theoretical_models("V18", "FSI", 75)      
    pm_avg11, red_pwiaXsec_CD_75 = read_theoretical_models("CD-Bonn", "PWIA", 75)                                 
    pm_avg12, red_fsiXsec_CD_75 = read_theoretical_models("CD-Bonn", "FSI", 75)

    #Convert to sig_reduced from MeV^-3 to fm^3
    red_pwiaXsec_avg_35 = red_pwiaXsec_avg[thnq==35]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_35 = red_fsiXsec_avg[thnq==35]*MeV2fm      #Laget FSI
    # red_pwiaXsec_avg_35 = red_pwiaXsec_avg_35.interpolate(method='cubicspline',limit_area='inside')
    # red_fsiXsec_avg_35 = red_fsiXsec_avg_35.interpolate(method='cubicspline',limit_area='inside')
    
    red_pwiaXsec_V18_35 = red_pwiaXsec_V18_35*MeV2fm
    red_fsiXsec_V18_35 = red_fsiXsec_V18_35*MeV2fm
    red_pwiaXsec_CD_35 = red_pwiaXsec_CD_35*MeV2fm
    red_fsiXsec_CD_35 = red_fsiXsec_CD_35*MeV2fm

    red_pwiaXsec_avg_45 = red_pwiaXsec_avg[thnq==45]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_45 = red_fsiXsec_avg[thnq==45]*MeV2fm      #Laget FSI
    red_pwiaXsec_avg_45 = red_pwiaXsec_avg_45.interpolate(method='cubicspline',limit_area='inside')
    red_fsiXsec_avg_45 = red_fsiXsec_avg_45.interpolate(method='cubicspline',limit_area='inside')

    
    red_pwiaXsec_V18_45 = red_pwiaXsec_V18_45*MeV2fm
    red_fsiXsec_V18_45 = red_fsiXsec_V18_45*MeV2fm
    red_pwiaXsec_CD_45 = red_pwiaXsec_CD_45*MeV2fm
    red_fsiXsec_CD_45 = red_fsiXsec_CD_45*MeV2fm

    red_pwiaXsec_avg_75 = red_pwiaXsec_avg[thnq==75]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_75 = red_fsiXsec_avg[thnq==75]*MeV2fm      #Laget FSI
    red_pwiaXsec_avg_75 = red_pwiaXsec_avg_75.interpolate(method='cubicspline',limit_area='inside')
    red_fsiXsec_avg_75 = red_fsiXsec_avg_75.interpolate(method='cubicspline',limit_area='inside')

    
    red_pwiaXsec_V18_75 = red_pwiaXsec_V18_75*MeV2fm
    red_fsiXsec_V18_75 = red_fsiXsec_V18_75*MeV2fm
    red_pwiaXsec_CD_75 = red_pwiaXsec_CD_75*MeV2fm
    red_fsiXsec_CD_75 = red_fsiXsec_CD_75*MeV2fm

    #Do Interpolation to make theory curves smooth
    #35 deg
    f_red_pwiaXsec_avg_35 = interp1d(pm_avg[thnq==35], red_pwiaXsec_avg_35,fill_value='extrapolate')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_35 = interp1d(pm_avg[thnq==35], red_fsiXsec_avg_35,fill_value='extrapolate')  
    # print(f_red_pwiaXsec_avg_35(pm_avg[thnq==35]))
 
    # f_red_pwiaXsec_avg_35 = interp1d(pm[thnq==35], red_pwiaXsec_avg_35,fill_value='extrapolate', kind='cubic')   #Paris (J.M. Laget Calculation)
    # f_red_fsiXsec_avg_35 = interp1d(pm[thnq==35], red_fsiXsec_avg_35,fill_value='extrapolate', kind='cubic')    
    
    f_red_pwiaXsec_V18_35 = interp1d(pm_avg1, red_pwiaXsec_V18_35,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_35 = interp1d(pm_avg2, red_fsiXsec_V18_35,fill_value='extrapolate', kind='cubic')
    
    f_red_pwiaXsec_CD_35 = interp1d(pm_avg3, red_pwiaXsec_CD_35,fill_value='extrapolate', kind='cubic')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_35 = interp1d(pm_avg4, red_fsiXsec_CD_35,fill_value='extrapolate', kind='cubic')
    #45 eg
    f_red_pwiaXsec_avg_45 = interp1d(pm_avg[thnq==45], red_pwiaXsec_avg_45,fill_value='extrapolate')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_45 = interp1d(pm_avg[thnq==45], red_fsiXsec_avg_45,fill_value='extrapolate')    
    
    f_red_pwiaXsec_V18_45 = interp1d(pm_avg5, red_pwiaXsec_V18_45,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_45 = interp1d(pm_avg6, red_fsiXsec_V18_45,fill_value='extrapolate', kind='cubic')
    
    f_red_pwiaXsec_CD_45 = interp1d(pm_avg7, red_pwiaXsec_CD_45,fill_value='extrapolate', kind='cubic')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_45 = interp1d(pm_avg8, red_fsiXsec_CD_45,fill_value='extrapolate', kind='cubic')
    #75 eg
    f_red_pwiaXsec_avg_75 = interp1d(pm_avg[thnq==75], red_pwiaXsec_avg_75,fill_value='extrapolate', kind='linear')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_75 = interp1d(pm_avg[thnq==75], red_fsiXsec_avg_75,fill_value='extrapolate', kind='linear')    
    
    f_red_pwiaXsec_V18_75 = interp1d(pm_avg9, red_pwiaXsec_V18_75,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_75 = interp1d(pm_avg10, red_fsiXsec_V18_75,fill_value='extrapolate', kind='cubic')
    
    f_red_pwiaXsec_CD_75 = interp1d(pm_avg11, red_pwiaXsec_CD_75,fill_value='extrapolate', kind='quadratic')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_75 = interp1d(pm_avg12, red_fsiXsec_CD_75,fill_value='extrapolate', kind='quadratic')

    

    #------Get Averaged Missing Momentum-----
    # pmiss_avg_35 = pm_avg[thnq==35].interpolate(method='cubicspline',limit_area='inside')
    # pmiss_avg_45 = pm_avg[thnq==45].interpolate(method='cubicspline',limit_area='inside')
    # pmiss_avg_75 = pm_avg[thnq==75].interpolate(method='cubicspline',limit_area='inside')
    
    pmiss_avg_35 = pm_avg[thnq==35]
    pmiss_avg_45 = pm_avg[thnq==45]
    pmiss_avg_75 = pm_avg[thnq==75]
      

    '''
    #====================================================PRL PLOT 1=============================================================
    
    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(13.4, 6)) 
    gs = gridspec.GridSpec(1, 3) 
    plt.xticks([])
    plt.yticks([])

    
    #THETA_NQ = 35 DEG
    ax0 = B.pl.subplot(gs[0])

    B.pl.text(0.25, 1e-7, r'$\theta_{nq}=35\pm5^{o}$', fontsize=19)
    B.pl.text(0.1, 1e-7, r'(a)', fontsize=19)
    
    
    l2 = B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data')
    l22 = B.plot_exp(pm_hc35, red_dataXsec_hc35, red_dataXsec_err_hc35, marker='s',  markersize=5, color='c', markerfacecolor='white', capsize=0, logy=True,  label='Hall C 2017 Data')
    # l1 = B.plot_exp(pm_avg[thnq==35], red_dataXsec_avg[thnq==35]*MeV2fm, red_dataXsec_avg_tot_err[thnq==35]*MeV2fm, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    l1 = B.plot_exp(pm_avg[thnq==35], red_dataXsec_avg_masked[thnq==35]*MeV2fm, red_dataXsec_avg_tot_err[thnq==35]*MeV2fm, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.pl.hlines(3.56356E-06, 0.920, 0.960, colors='k', linestyles='-')
    #B.pl.hlines(3.62178722, 0.920, 0.960, colors='k', linestyles='-')
    #Plot theoretical curves
    # print('f_red_pwiaXsec_avg_35(pm_avg[thnq==35])=',f_red_pwiaXsec_avg_35(pm_avg[thnq==35]))
    # print('f_red_fsiXsec_avg_35(pm_avg[thnq==35])=',f_red_fsiXsec_avg_35(pm_avg[thnq==35]))

    l3, = B.plot_exp(pm_avg[thnq==35], f_red_pwiaXsec_avg_35(pm_avg[thnq==35]), linestyle='--', marker='None', color='#0000ff', logy=True, label='JML Paris PWIA')
    l4, = B.plot_exp(pm_avg[thnq==35], f_red_fsiXsec_avg_35(pm_avg[thnq==35]), linestyle='-', marker='None', color='#0000ff', logy=True, label='JML Paris FSI')
    
    l5, = B.plot_exp(pm_avg1, f_red_pwiaXsec_V18_35(pm_avg1), linestyle='--', marker='None', color='#009000', logy=True, label='MS AV18 PWIA')   
    l6, = B.plot_exp(pm_avg2, f_red_fsiXsec_V18_35(pm_avg2), linestyle='-', marker='None', color='#009000', logy=True, label='MS AV18 FSI') 
    
    l7, = B.plot_exp(pm_avg3, f_red_pwiaXsec_CD_35(pm_avg3), linestyle='--', marker='None', color='#ff00ff', logy=True, label='MS CD-Bonn PWIA')     
    l8, = B.plot_exp(pm_avg4, f_red_fsiXsec_CD_35(pm_avg4), linestyle='-', marker='None', color='#ff00ff', logy=True, label='MS CD-Bonn FSI') 


    #Set axis limits
    B.pl.xlim(0.0, 1.2)

    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$\sigma_{\mathrm{red}}$ ($\mathrm{fm}^{3}$) ', fontsize=22,  labelpad=10)                                                                                                                                                                        
    B.pl.title('') 

    #THETA_NQ = 45 DEG
    ax1 = plt.subplot(gs[1], sharey=ax0)
    
    B.pl.text(0.25, 1e-7, r'$\theta_{nq}=45\pm5^{o}$', fontsize=19)
    B.pl.text(0.1, 1e-7, r'(b)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax1.get_yticklabels(), visible=False)

    B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s', markersize=5, capsize=0, color='#ff0000', markerfacecolor='white',  logy=True, label='Hall A Data')
    B.plot_exp(pm_hc45, red_dataXsec_hc45, red_dataXsec_err_hc45, marker='s',  markersize=5, color='c', markerfacecolor='white', capsize=0, logy=True,  label='Hall C 2017 Data')
    B.plot_exp(pm_avg[thnq==45], red_dataXsec_avg_masked[thnq==45]*MeV2fm, red_dataXsec_avg_tot_err[thnq==45]*MeV2fm, marker='o', markersize=5, capsize=0, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    #B.pl.hlines(6.96669E-06, 0.980, 1.02, colors='k', linestyles='-') 

    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==45], f_red_pwiaXsec_avg_45(pm_avg[thnq==45]), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==45], f_red_fsiXsec_avg_45(pm_avg[thnq==45]), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg5, f_red_pwiaXsec_V18_45(pm_avg5), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg6, f_red_fsiXsec_V18_45(pm_avg6), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg7, f_red_pwiaXsec_CD_45(pm_avg7), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg8, f_red_fsiXsec_CD_45(pm_avg8), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI') 

    #Set axis limits
    B.pl.xlim(0.0, 1.2)

    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 1
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('')

    #THETA_NQ = 75
    ax2 = plt.subplot(gs[2], sharey=ax0)

    B.pl.text(0.25, 1e-7, r'$\theta_{nq}=75\pm5^{o}$', fontsize=19)
    B.pl.text(0.1, 1e-7, r'(c)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax2.get_yticklabels(), visible=False)

    B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s', markersize=5, capsize=0, color='#ff0000',  markerfacecolor='white', logy=True, label='Hall A Data')
    B.plot_exp(pm_hc75, red_dataXsec_hc75, red_dataXsec_err_hc75, marker='s',  markersize=5, color='c', markerfacecolor='white', capsize=0, logy=True,  label='Hall C 2017 Data')
    B.plot_exp(pm_avg[thnq==75], red_dataXsec_avg_masked[thnq==75]*MeV2fm, red_dataXsec_avg_tot_err[thnq==75]*MeV2fm, marker='o', markersize=5, capsize=0, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    #B.pl.hlines(1.41619E-03, 0.360, 0.400, colors='k', linestyles='-')
    #Plot theoretical curves
    B.plot_exp(pm[thnq==75], f_red_pwiaXsec_avg_75(pm[thnq==75]), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA')
    B.plot_exp(pm[thnq==75], f_red_fsiXsec_avg_75(pm[thnq==75]), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg9, f_red_pwiaXsec_V18_75(pm_avg9), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg10, f_red_fsiXsec_V18_75(pm_avg10), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg11, f_red_pwiaXsec_CD_75(pm_avg11), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg12, f_red_fsiXsec_CD_75(pm_avg12), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI') 

    #Set axis limits
    B.pl.xlim(0.0, 1.2)

    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 2
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('')
  
    #Remove spacing between subplots
    plt.subplots_adjust(wspace = 0.0000000001, bottom=0.13, top=0.98, left=0.085, right=0.98) #, hspace = 0.001, wspace = 0.001)
    
    #Create labels for specified plots in given order
    line_labels=['This Experiment (Hall C)', 'Hall A Data', 'Hall C (2017) Data', 'Paris PWIA', 'Paris FSI', 'AV18 PWIA', 'AV18 FSI', 'CD-Bonn PWIA', 'CD-Bonn FSI']   #define legend labels
    ax2.legend([l1, l2, l22, l3, l4, l5, l6, l7, l8], line_labels, loc='upper right', frameon=False, fontsize=12)      #subplot to use for common legend
      

    #B.pl.show()
    B.pl.savefig(dir_name+'/PRL_plot1.pdf')
    '''
    '''
    B.pl.figure()
    #Take Ratio of red. Xsec to Misak's redXsec to determine K_sig_cc1 ratio
    pm_paris, redXsec_Paris_MS = read_misak_pdist(model="Paris")
    pm_av18, redXsec_AV18_MS = read_misak_pdist(model="AV18")
    pm_cd, redXsec_CD_MS = read_misak_pdist(model="CD-Bonn")

    R_cd_35 =  f_red_pwiaXsec_CD_35(pmiss_avg_35) / redXsec_CD_MS(pmiss_avg_35) #f_red_pwiaXsec_CD_35(pm_cd) / redXsec_CD_MS(pm_cd)
    R_cd_45 =  f_red_pwiaXsec_CD_45(pmiss_avg_45) / redXsec_CD_MS(pmiss_avg_45)
    R_cd_75 =  f_red_pwiaXsec_CD_75(pmiss_avg_75) / redXsec_CD_MS(pmiss_avg_75)
    
    R_paris_35 =  f_red_pwiaXsec_avg_35(pmiss_avg_35) / redXsec_Paris_MS(pmiss_avg_35) #f_red_pwiaXsec_CD_35(pm_cd) / redXsec_CD_MS(pm_cd)
    R_paris_45 =  f_red_pwiaXsec_avg_45(pmiss_avg_45) / redXsec_Paris_MS(pmiss_avg_45)
    R_paris_75 =  f_red_pwiaXsec_avg_75(pmiss_avg_75) / redXsec_Paris_MS(pmiss_avg_75)
    
    R_AV_35 =  f_red_pwiaXsec_V18_35(pmiss_avg_35) / redXsec_AV18_MS(pmiss_avg_35) #f_red_pwiaXsec_CD_35(pm_cd) / redXsec_CD_MS(pm_cd)
    R_AV_45 =  f_red_pwiaXsec_V18_45(pmiss_avg_45) / redXsec_AV18_MS(pmiss_avg_45)
    R_AV_75 =  f_red_pwiaXsec_V18_75(pmiss_avg_75) / redXsec_AV18_MS(pmiss_avg_75)

    R_ha =  red_dataXsec_ha45 / redXsec_CD_MS(pm_ha45)
    R_ha_err =  red_dataXsec_err_ha45 / redXsec_CD_MS(pm_ha45)

    R_hc = red_dataXsec_avg_masked[thnq==45]*MeV2fm / redXsec_CD_MS(pm_avg[thnq==45])
    R_hc_err = red_dataXsec_avg_tot_err[thnq==75]*MeV2fm / redXsec_CD_MS(pm_avg[thnq==45])

    R_ref =  redXsec_CD_MS(pm_ha45) /  redXsec_CD_MS(pm_ha45)

    #B.plot_exp(pm_ha45, R_ha, R_ha_err, marker='o', markersize=5, capsize=0, color='#ff0000',   label='Hall A Data / CD-Bonn (Misak)')
    #B.plot_exp(pm_avg[thnq==45], R_hc, R_hc_err, marker='o', markersize=5, capsize=0, color='k', label='Hall C Data / CD-Bonn (Misak)')
    #B.plot_exp(pm_ha45, R_ref,  marker='None', linestyle='--', color='#ff00ff',  label='Ref.')

    B.pl.yscale('log')
    
    
    B.plot_exp(pmiss_avg_35, R_cd_35, marker='o', linestyle='--', color='r', label=r'$CD-Bonn: 35 deg$')
    B.plot_exp(pmiss_avg_45, R_cd_45, marker='^', linestyle='--', color='r', label=r'$CD-Bonn: 45 deg$')
    
    B.plot_exp(pmiss_avg_35, R_paris_35, marker='o', linestyle='--', color='b', label=r'$Paris: 35 deg$')
    B.plot_exp(pmiss_avg_45, R_paris_45, marker='^', linestyle='--', color='b', label=r'$Paris: 45 deg$')
        
    B.plot_exp(pmiss_avg_35, R_AV_35, marker='o', linestyle='--', color='g', label=r'$AV18: 35 deg$')
    B.plot_exp(pmiss_avg_45, R_AV_45, marker='^', linestyle='--', color='g', label=r'$AV18: 45 deg$')
    
    #B.plot_exp(pm_cd, redXsec_CD_MS(pm_cd), linestyle='--', marker='None', color='#ff00ff', logy=False, label='CD-Bonn PWIA MS')     
    #B.plot_exp(pm_paris, redXsec_Paris_MS(pm_paris), linestyle='--', marker='None', color='#0000ff', logy=False, label='Paris PWIA MS')     
    #B.plot_exp(pm_av18, redXsec_AV18_MS(pm_av18), linestyle='--', marker='None', color='g', logy=False, label='AV18 PWIA MS')     
    
    #B.plot_exp(pm_cd, f_red_pwiaXsec_CD_35(pm_cd), linestyle='-', marker='None', color='#ff00ff', logy=False, label='CD-Bonn PWIA')   
    #B.plot_exp(pm_paris, f_red_pwiaXsec_avg_35(pm_paris), linestyle='-', marker='None', color='#0000ff', logy=False, label='Paris PWIA')   
    #B.plot_exp(pm_av18, f_red_pwiaXsec_V18_35(pm_av18), linestyle='-', marker='None', color='g', logy=False, label='AV18 PWIA')   

    B.pl.xlabel(r'$P_{m}$ [GeV/c]')
    B.pl.ylabel(r'R = $\sigma_{red} / \sigma^{Misak}_{red}$')
    B.pl.legend()
    #B.pl.show()
    
    '''

    #Write PRL PLOT 1 data to file (avg. Pm, data Xsec and errors, and theory Xsec
    '''
    fout_name = 'redXsec_HallC_thnq35_deg.txt'
    fout = open(fout_name, 'w')
    comment1='#This datafile contains redXsec from Hall C Deuteron Experiment: E12-10-003\n'
    comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = 35 +\- 5 deg \n'
    header='#!pm_avg[f,0]/   data_redXsec[f,1]/   data_redXsec_tot_err[f,2]/   stats_rel_err[f,3]/  kin_rel_err[f,4]/   norm_rel_err[f,5]   tot_rel_syst_err[f,6]/   tot_rel_err[f,7]/ \n'
    fout.write(comment1)
    fout.write(comment2)
    fout.write(header)
    
    for i in range(len(pm_avg[thnq==35])):
        fout.write('%.5f %.5E  %.5E  %.5f  %.5f  %.5f  %.5f  %.5f\n' % ((pm_avg[thnq==35])[i], (red_dataXsec_avg_masked[thnq==35]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==35]*MeV2fm)[i], (tot_stats_err[thnq==35])[i], (kin_syst_tot[thnq==35])[i], (norm_syst_tot[thnq==35])[i], (tot_syst_err[thnq==35])[i], (tot_err[thnq==35])[i]))

    fout.close()

    fout_name = 'redXsec_HallC_thnq45_deg.txt'
    fout = open(fout_name, 'w')
    comment1='#This datafile contains redXsec from Hall C Deuteron Experiment: E12-10-003\n'
    comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = 45 +\- 5 deg \n'
    header='#!pm_avg[f,0]/   data_redXsec[f,1]/   data_redXsec_tot_err[f,2]/   stats_rel_err[f,3]/  kin_rel_err[f,4]/   norm_rel_err[f,5]   tot_rel_syst_err[f,6]/   tot_rel_err[f,7]/ \n'
    fout.write(comment1)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==45])):
        fout.write('%.5f %.5E  %.5E  %.5f  %.5f  %.5f  %.5f  %.5f\n' % ((pm_avg[thnq==45])[i], (red_dataXsec_avg_masked[thnq==45]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==45]*MeV2fm)[i], (tot_stats_err[thnq==45])[i], (kin_syst_tot[thnq==45])[i], (norm_syst_tot[thnq==45])[i], (tot_syst_err[thnq==45])[i], (tot_err[thnq==45])[i]))

    fout.close()


    fout_name = 'redXsec_HallC_thnq75_deg.txt'
    fout = open(fout_name, 'w')
    comment1='#This datafile contains redXsec from Hall C Deuteron Experiment: E12-10-003\n'
    comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = 75 +\- 5 deg \n'
    header='#!pm_avg[f,0]/   data_redXsec[f,1]/   data_redXsec_tot_err[f,2]/   stats_rel_err[f,3]/  kin_rel_err[f,4]/   norm_rel_err[f,5]   tot_rel_syst_err[f,6]/   tot_rel_err[f,7]/ \n'
    fout.write(comment1)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==75])):
        fout.write('%.5f %.5E  %.5E  %.5f  %.5f  %.5f  %.5f  %.5f\n' % ((pm_avg[thnq==75])[i], (red_dataXsec_avg_masked[thnq==75]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==75]*MeV2fm)[i], (tot_stats_err[thnq==75])[i], (kin_syst_tot[thnq==75])[i], (norm_syst_tot[thnq==75])[i], (tot_syst_err[thnq==75])[i], (tot_err[thnq==75])[i]))

    fout.close()
    '''
    
    
    #====================================================PRL PLOT 2=============================================================
    
    #Define the Data (and all models) to CD-Bonn PWIA Ratio
    #theta_nq = 35 deg
    R_ref35 = f_red_pwiaXsec_CD_35(pm[thnq==35]) / f_red_pwiaXsec_CD_35(pm[thnq==35])
    R_data35 = red_dataXsec_avg_masked[thnq==35]*MeV2fm / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_data_err35 = red_dataXsec_avg_tot_err[thnq==35]*MeV2fm / f_red_pwiaXsec_CD_35(pmiss_avg_35)

    R_HAdata35 = red_dataXsec_ha35 / f_red_pwiaXsec_CD_35(pm_ha35)   #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
    R_HAdata_err35 =  red_dataXsec_err_ha35 / f_red_pwiaXsec_CD_35(pm_ha35)
    
    R_HCdata35 = red_dataXsec_hc35 / f_red_pwiaXsec_CD_35(pm_hc35)   #Hall C data 2017
    R_HCdata_err35 =  red_dataXsec_err_hc35 / f_red_pwiaXsec_CD_35(pm_hc35)


    R_CDBonn_fsi35 = f_red_fsiXsec_CD_35(pm[thnq==35]) / f_red_pwiaXsec_CD_35(pm[thnq==35])
    # R_Paris_pwia35 = f_red_pwiaXsec_avg_35(pm[thnq==35]) / f_red_pwiaXsec_CD_35(pm[thnq==35])
    # R_Paris_fsi35 = f_red_fsiXsec_avg_35(pm[thnq==35]) / f_red_pwiaXsec_CD_35(pm[thnq==35])
    # R_Paris_pwia35 = red_pwiaXsec_avg_35 / f_red_pwiaXsec_CD_35(pm[thnq==35])
    # R_Paris_fsi35 = red_fsiXsec_avg_35 / f_red_pwiaXsec_CD_35(pm[thnq==35])
    R_Paris_pwia35 = red_pwiaXsec_avg_35 / red_pwiaXsec_CD_35
    R_Paris_fsi35 = red_fsiXsec_avg_35 / red_pwiaXsec_CD_35
    R_Paris_pwia35[98] = 1.13
    R_Paris_fsi35[98] = 1.10
    R_Paris_pwia35[99] = 1.18
    R_Paris_fsi35[99] = 1.15
    print('R_Paris_fsi35 = \n',R_Paris_fsi35)
    print('R_Paris_pwia35 = \n',R_Paris_pwia35)
    R_AV18_pwia35 = f_red_pwiaXsec_V18_35(pm[thnq==35]) / f_red_pwiaXsec_CD_35(pm[thnq==35])
    R_AV18_fsi35 = f_red_fsiXsec_V18_35(pm[thnq==35]) / f_red_pwiaXsec_CD_35(pm[thnq==35])

    #theta_nq = 45 deg    
    R_ref45 = f_red_pwiaXsec_CD_45(pm[thnq==45]) / f_red_pwiaXsec_CD_45(pm[thnq==45])
    R_data45 = red_dataXsec_avg_masked[thnq==45]*MeV2fm / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_data_err45 = red_dataXsec_avg_tot_err[thnq==45]*MeV2fm / f_red_pwiaXsec_CD_45(pmiss_avg_45)

    R_HAdata45 = red_dataXsec_ha45 / f_red_pwiaXsec_CD_45(pm_ha45)   #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
    R_HAdata_err45 =  red_dataXsec_err_ha45 / f_red_pwiaXsec_CD_45(pm_ha45)

    R_HCdata45 = red_dataXsec_hc45 / f_red_pwiaXsec_CD_45(pm_hc35)   #Hall C data 2017
    R_HCdata_err45 =  red_dataXsec_err_hc45 / f_red_pwiaXsec_CD_35(pm_hc45)

    R_CDBonn_fsi45 = f_red_fsiXsec_CD_45(pm[thnq==45]) / f_red_pwiaXsec_CD_45(pm[thnq==45])
    R_Paris_pwia45 = f_red_pwiaXsec_avg_45(pm[thnq==45]) / f_red_pwiaXsec_CD_45(pm[thnq==45])
    R_Paris_fsi45 = f_red_fsiXsec_avg_45(pm[thnq==45]) / f_red_pwiaXsec_CD_45(pm[thnq==45])
    # R_Paris_pwia45 = red_pwiaXsec_avg_45 / f_red_pwiaXsec_CD_45(pm[thnq==45])
    # R_Paris_fsi45 = red_fsiXsec_avg_45 / f_red_pwiaXsec_CD_45(pm[thnq==45])
    R_AV18_pwia45 = f_red_pwiaXsec_V18_45(pm[thnq==45]) / f_red_pwiaXsec_CD_45(pm[thnq==45])
    R_AV18_fsi45 = f_red_fsiXsec_V18_45(pm[thnq==45]) / f_red_pwiaXsec_CD_45(pm[thnq==45])

    #theta_nq = 75 deg    
    R_ref75 = f_red_pwiaXsec_CD_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_data75 = red_dataXsec_avg_masked[thnq==75]*MeV2fm / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_data_err75 = red_dataXsec_avg_tot_err[thnq==75]*MeV2fm / f_red_pwiaXsec_CD_75(pmiss_avg_75)

    R_HAdata75 = red_dataXsec_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)   #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
    R_HAdata_err75 =  red_dataXsec_err_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)

    R_HCdata75 = red_dataXsec_hc75 / f_red_pwiaXsec_CD_75(pm_hc75)   #Hall C data 2017
    R_HCdata_err75 =  red_dataXsec_err_hc75 / f_red_pwiaXsec_CD_35(pm_hc75)

    R_CDBonn_fsi75 = f_red_fsiXsec_CD_75(pm_avg[thnq==75]) / f_red_pwiaXsec_CD_75(pm_avg[thnq==75])
    R_Paris_pwia75 = f_red_pwiaXsec_avg_75(pm_avg[thnq==75]) / f_red_pwiaXsec_CD_75(pm_avg[thnq==75])
    R_Paris_fsi75 = f_red_fsiXsec_avg_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    # R_Paris_pwia75 = red_pwiaXsec_avg_75 / f_red_pwiaXsec_CD_75(pm_avg[thnq==75])
    # R_Paris_fsi75 = red_fsiXsec_avg_75 / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_pwia75 = f_red_pwiaXsec_V18_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_fsi75 = f_red_fsiXsec_V18_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)


    #B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='o', color='r', label='Hall A Data, 35 deg', capsize=0)
    #B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='^', color='g', label='Hall A Data, 45 deg', capsize=0)
    #B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='b', label='Hall A Data, 75 deg', capsize=0)
    #B.pl.ylabel('Ratio')
    #B.pl.xlabel('GeV/c')
    #B.pl.title('Hall A / CD-Bonn PWIA red. Xsec (Updated Hall A data)')
    #B.pl.legend()

    #-----Create Subplots-----
    # fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(6.7, 10)) 
    # gs = gridspec.GridSpec(3, 1)
    B.pl.figure(figsize=(8,6.5))
    # plt.xticks([])
    # plt.yticks([])
    
    # B.pl.text(0.1, 4, 'PRELIMINARY', color='gray',alpha=0.4,
    #           fontsize=48)

    #THETA_NQ = 35 DEG
    # ax0 = B.pl.subplot(gs[0])
    ax0 = B.pl.gca()

    #B.pl.text(0.1, 4, r'(a)', fontsize=19)  #original txt label
    # B.pl.text(0.1, 2, r'(a)', fontsize=19)  #ZOOM plot txt label
    #B.pl.text(1.1, 8, r'(a)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    # B.pl.setp(ax0.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.plot_exp(pm[thnq==35], R_Paris_pwia35, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA')
    B.plot_exp(pm[thnq==35], R_Paris_fsi35, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI')
    
    B.plot_exp(pm[thnq==35], R_AV18_pwia35, marker='None', linestyle='--', color='#009000', label='MS AV18 PWIA')
    B.plot_exp(pm[thnq==35], R_AV18_fsi35, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI')
    
    B.plot_exp(pm[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff', label='MS CD-Bonn PWIA')
       
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0)
    B.plot_exp(pm_hc35, R_HCdata35, R_HCdata_err35, marker='s', color='c', markerfacecolor='white',  label='Hall C (2017) Data', capsize=0)
    B.plot_exp(pm_avg[thnq==35], R_data35, R_data_err35, marker='o', color='k', label='New Data (Preliminary)')
    #print('pm, Rdata_35 = ',pm_avg[thnq==35],'::',R_data35)
    #B.pl.hlines(3.62178722, 0.920, 0.960, colors='k', linestyles='-') 
    B.plot_exp(pm[thnq==35], R_CDBonn_fsi35, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI')
    

    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    #B.pl.ylim(0.5, 2.3) #ZOOM 
    B.pl.ylim(-0.5, 10)  #ORIGINAL

    #log settings
    #B.pl.ylim(0.5, 2.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

    #Set Axes Labels for subplot 0
    # B.pl.xlabel('')                                                                                                                                                                                                   
    # B.pl.ylabel('')                                                                                                                                                                        
    # B.pl.title('')
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)
    B.pl.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=22,  labelpad=10)                                                                                                                                                                        
    B.pl.title('')

    B.pl.legend(loc='best', fontsize=12, frameon=False)
    
    #----CRTE INSET PLOT(35 DEG)
    '''
    axins35 = inset_axes(ax0, width='50%', height='50%', loc='upper left')
    B.pl.xlim(0,0.8)
    B.pl.ylim(0.3,2)
    B.plot_exp(pm_avg[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff')
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='o', color='#ff0000',   capsize=0)
    B.plot_exp(pm_avg[thnq==35], R_data35, R_data_err35, marker='o', color='k', capsize=0)

    B.plot_exp(pm_avg[thnq==35], R_CDBonn_fsi35, marker='None', linestyle='-', color='#ff00ff')
    B.plot_exp(pm_avg[thnq==35], R_Paris_pwia35, marker='None', linestyle='--', color='#0000ff')
    B.plot_exp(pm_avg[thnq==35], R_Paris_fsi35, marker='None', linestyle='-', color='#0000ff')
    
    B.plot_exp(pm_avg[thnq==35], R_AV18_pwia35, marker='None', linestyle='--', color='#009000')
    B.plot_exp(pm_avg[thnq==35], R_AV18_fsi35, marker='None', linestyle='-', color='#009000')
    #axins35.tick_params(labelleft=False, labelbottom=False)
    axins35.yaxis.tick_right()
    
    
    B.pl.xlabel('')
    B.pl.ylabel('')
    B.pl.title('')
    '''
    
    B.pl.savefig(dir_name+'/Ratio_35.png')

    #THETA_NQ = 45 DEG
    B.pl.figure(figsize=(8,6.5))
    # ax1 = B.pl.subplot(gs[1])
    ax1 = B.pl.gca()
    # B.pl.text(0.1, 5, 'PRELIMINARY', color='gray',alpha=0.4,
    #           fontsize=48)

    # B.pl.text(0.1, 16, r'(b)', fontsize=19)  #original  txt
    #B.pl.text(0.1, 1.8, r'(b)', fontsize=19)   #new txt
    #B.pl.text(1.1, 16, r'(b)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    # B.pl.setp(ax1.get_xticklabels(), visible=False)
    
    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.plot_exp(pm[thnq==45], R_Paris_pwia45, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA')
    B.plot_exp(pm[thnq==45], R_Paris_fsi45, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI')
    
    B.plot_exp(pm[thnq==45], R_AV18_pwia45, marker='None', linestyle='--', color='#009000', label='MS AV18 PWIA')
    B.plot_exp(pm[thnq==45], R_AV18_fsi45, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI')
 
    B.plot_exp(pm[thnq==45], R_ref45,  marker='None', linestyle='--', color='#ff00ff', label='MS CD-Bonn PWIA')
    
    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0)
    B.plot_exp(pm_hc45, R_HCdata45, R_HCdata_err45, marker='s', color='c', markerfacecolor='white',  label='Hall C (2017) Data', capsize=0)
    B.plot_exp(pm_avg[thnq==45], R_data45, R_data_err45, marker='o', color='k', label='New Data (Preliminary)')
    #print('pm, Rdata_45 = ',pm_avg[thnq==45],'::',R_data45)                                                                   
    #B.pl.hlines(16.65557, 0.980, 1.02, colors='k', linestyles='-')  

    B.plot_exp(pm[thnq==45], R_CDBonn_fsi45, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI')



    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-0.1, 20) #ZOOM 
    #B.pl.ylim(-0.5, 22.9) #ORIGINAL

    #log settings
    #B.pl.ylim(0.5, 2.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    # B.pl.xlabel('')                                                                                                                                                                                                   
    # B.pl.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=22,  labelpad=10)                                                                                                                                                                        
    # B.pl.title('') 
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)
    B.pl.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=22,  labelpad=10)                                                                                                                                                                        
    B.pl.title('')
    
    B.pl.legend(loc='best', fontsize=12, frameon=False)

    '''
    #----CREATE INSET PLOT(45 DEG)
    axins45 = inset_axes(ax1, width='50%', height='50%', loc='upper left')
    B.pl.xlim(0,0.8)
    B.pl.ylim(0.3,2)
    B.plot_exp(pm_avg[thnq==45], R_ref45,  marker='None', linestyle='--', color='#ff00ff')
    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='o', color='#ff0000',   capsize=0)
    B.plot_exp(pm_avg[thnq==45], R_data45, R_data_err45, marker='o', color='k', capsize=0)

    B.plot_exp(pm_avg[thnq==45], R_CDBonn_fsi45, marker='None', linestyle='-', color='#ff00ff')
    B.plot_exp(pm_avg[thnq==45], R_Paris_pwia45, marker='None', linestyle='--', color='#0000ff')
    B.plot_exp(pm_avg[thnq==45], R_Paris_fsi45, marker='None', linestyle='-', color='#0000ff')
    
    B.plot_exp(pm_avg[thnq==45], R_AV18_pwia45, marker='None', linestyle='--', color='#009000')
    B.plot_exp(pm_avg[thnq==45], R_AV18_fsi45, marker='None', linestyle='-', color='#009000')
    #axins45.tick_params(labelleft=False, labelbottom=False)
    axins45.yaxis.tick_right()
    

    B.pl.xlabel('')
    B.pl.ylabel('')
    B.pl.title('')
    '''
    B.pl.savefig(dir_name+'/Ratio_45.png')

    #THETA_NQ = 75 DEG
    B.pl.figure(figsize=(8,6.5))
    # ax2 = B.pl.subplot(gs[2])
    ax2 = B.pl.gca()
    
    # B.pl.text(0.1, 6, 'PRELIMINARY', color='gray',alpha=0.4,
    #           fontsize=48)
    
    # B.pl.text(0.1, 4, r'(c)', fontsize=19)   # original 
    #B.pl.text(0.1, 4, r'(c)', fontsize=19)   # new txt

    #B.pl.text(1.1, 15, r'(c)', fontsize=19)  #with inset txt label


    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.plot_exp(pm_avg[thnq==75], R_Paris_pwia75, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA')
    B.plot_exp(pm_avg[thnq==75], R_Paris_fsi75, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI')
    
    B.plot_exp(pm[thnq==75], R_AV18_pwia75, marker='None', linestyle='--', color='#009000', label='MS AV18 PWIA')
    B.plot_exp(pm[thnq==75], R_AV18_fsi75, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI')
 
    B.plot_exp(pm[thnq==75], R_ref75,  marker='None', linestyle='--', color='#ff00ff', label='MS CD-Bonn PWIA')
    
    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0)
    B.plot_exp(pm_hc75, R_HCdata75, R_HCdata_err75, marker='s', color='c', markerfacecolor='white',  label='Hall C (2017) Data', capsize=0)
    B.plot_exp(pm_avg[thnq==75], R_data75, R_data_err75, marker='o', color='k', label='New Data (Preliminary)')

    #print('pm, Rdata_75 = ',pm_avg[thnq==75],'::',R_data75)                                                                     
    #B.pl.hlines(3.86100677, 0.360, 0.400, colors='k', linestyles='-')   

    B.plot_exp(pm[thnq==75], R_CDBonn_fsi75, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI')


    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-0.5, 12)

    #log settings
    #B.pl.ylim(0.1, 5.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    # B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                               
    # B.pl.ylabel('')                                                                                 
    # B.pl.title('') 
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)
    B.pl.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=22,  labelpad=10)                                                                                                                                                                        
    B.pl.title('')
    
    #Remove spacing between subplots
    # plt.subplots_adjust(hspace = 0.00000001, bottom=0.09, top=0.98, right=0.95, left=0.18) #, hspace = 0.001, wspace = 0.001)
    


    B.pl.legend(loc='best', fontsize=12, frameon=False)
    '''
    #----CREATE INSET PLOT(75 DEG)
    axins75 = inset_axes(ax2, width='50%', height='50%', loc='upper left')
    B.pl.xlim(0,0.6)
    B.pl.ylim(0.3,4)
    B.plot_exp(pm_avg[thnq==75], R_ref75,  marker='None', linestyle='--', color='#ff00ff')
    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='o', color='#ff0000',   capsize=0)
    B.plot_exp(pm_avg[thnq==75], R_data75, R_data_err75, marker='o', color='k', capsize=0)

    B.plot_exp(pm_avg[thnq==75], R_CDBonn_fsi75, marker='None', linestyle='-', color='#ff00ff')
    B.plot_exp(pm_avg[thnq==75], R_Paris_pwia75, marker='None', linestyle='--', color='#0000ff')
    B.plot_exp(pm_avg[thnq==75], R_Paris_fsi75, marker='None', linestyle='-', color='#0000ff')
    
    B.plot_exp(pm_avg[thnq==75], R_AV18_pwia75, marker='None', linestyle='--', color='#009000')
    B.plot_exp(pm_avg[thnq==75], R_AV18_fsi75, marker='None', linestyle='-', color='#009000')
    #axins75.tick_params(labelleft=False, labelbottom=False)
    axins75.yaxis.tick_right()
    '''

    # #B.pl.xlabel('')
    # #B.pl.ylabel('')
    # B.pl.title('')

    # B.pl.show()
    B.pl.savefig(dir_name+'/Ratio_75.png')
    # B.pl.savefig(dir_name+'/PRL_plot2.pdf')
    
    '''
    #Write to file
    fout_name = './yield_n_xsec/combine_data/redXsec_file_35.txt'
    fout = open(fout_name, 'w')
    comment='pm_avg in GeV/c :: redXsec in fm^3 :: R = red_dataXsec / redXsec_CD-Bonn_PWIA\n'
    comment2='Conv. factor: 1 fm^-1 = 0.1973 GeV\n'
    header='#!pm_avg[f,0]/ red_dataXsec_masked[f,1]/ redXsec_CD-Bonn_PWIA[f,2]/ R_data[f,3]/ R_err[f,4]/\n'
    fout.write(comment)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==35.])):
        # print('ith=',i,' :: pmiss_avg=',pm_avg[i])
        fout.write('%.5f  %.10f   %.10f   %.10f   %.10f\n' % (pm_avg[thnq==35.][i], red_dataXsec_avg_masked[thnq==35.][i], f_red_pwiaXsec_CD_35(pm_avg[thnq==35.])[i], R_data35[i], R_data_err35[i]))

    fout.close()
    
    fout_name = './yield_n_xsec/combine_data/redXsec_file_45.txt'
    fout = open(fout_name, 'w')
    comment='pm_avg in GeV/c :: redXsec in fm^3 :: R = red_dataXsec / redXsec_CD-Bonn_PWIA\n'
    comment2='Conv. factor: 1 fm^-1 = 0.1973 GeV\n'
    header='#!pm_avg[f,0]/ red_dataXsec_masked[f,1]/ redXsec_CD-Bonn_PWIA[f,2]/ R_data[f,3]/ R_err[f,4]/\n'
    fout.write(comment)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==45.])):
        # print('ith=',i,' :: pmiss_avg=',pm_avg[i])
        fout.write('%.5f  %.10f   %.10f   %.10f   %.10f\n' % (pm_avg[thnq==45.][i], red_dataXsec_avg_masked[thnq==45.][i], f_red_pwiaXsec_CD_35(pm_avg[thnq==45.])[i], R_data45[i], R_data_err45[i]))

    fout.close()
    
    fout_name = './yield_n_xsec/combine_data/redXsec_file_75.txt'
    fout = open(fout_name, 'w')
    comment='pm_avg in GeV/c :: redXsec in fm^3 :: R = red_dataXsec / redXsec_CD-Bonn_PWIA\n'
    comment2='Conv. factor: 1 fm^-1 = 0.1973 GeV\n'
    header='#!pm_avg[f,0]/ red_dataXsec_masked[f,1]/ redXsec_CD-Bonn_PWIA[f,2]/ R_data[f,3]/ R_err[f,4]/\n'
    fout.write(comment)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==75.])):
        print('ith=',i,' :: pmiss_avg=',pm_avg[i])
        fout.write('%.5f  %.10f   %.10f   %.10f   %.10f\n' % (pm_avg[thnq==75.][i], red_dataXsec_avg_masked[thnq==75.][i], f_red_pwiaXsec_CD_35(pm_avg[thnq==75.])[i], R_data75[i], R_data_err75[i]))

    fout.close()
    '''

    #===========================================================================================================================

    
def main():
    print('Entering Main . . .')

    #plot_data_sets()

    #plot_theory_sets('pwia')
    #plot_theory_sets('fsi')
    #plot_Xsec_vs_thnq()
    #plot_final()
    make_prl_plots()

if __name__=="__main__":
    main()


