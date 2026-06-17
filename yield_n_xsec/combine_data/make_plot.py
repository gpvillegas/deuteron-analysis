import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline as cspl
import pandas as pd
# import numpy.ma as ma
# import sys                                     
import os                                                                                                       
# from sys import argv  
import matplotlib
from matplotlib import rc
# from matplotlib import *
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)


# matplotlib.use('Agg')

#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 19
}
plt.rc('font', **font)

#Set font
csfont = {'fontname':'Times New Roman'}

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr

def isfinite(arr):
    # returns an array with only finite values
    is_finite = np.isfinite(arr)
    return arr[is_finite]


#Conversion factor:  1 fm = 1/ (197 MeV),   
#The reduced cross section is in MeV^-3 units
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#User Input (Dir. Name to store output)
# sys_ext = sys.argv[1]   

dir_name = 'yield_n_xsec/combine_data/plots'
# dir_name = sys_ext+"_plots"
# dir_name_misc = dir_name + "/misc_plots"        #directory to store miscellaneous plots (trigger rates, live time, etc. . .) #obtained from the report files
# print(dir_name)
# #check if directory exists, else creates it.
if not os.path.exists(dir_name):
    os.makedirs(dir_name)
# if not os.path.exists(dir_name_misc):
#     os.makedirs(dir_name_misc)
     

#Get Reduced Xsec Data File (for Q2: (4,5) GeV2)
fname = 'yield_n_xsec/combine_data/redXsec_combined.txt'
f = B.get_file(fname)

#Get Bin Information (Same info for all files)                                                                                                  
i_b = B.get_data(f, 'i_b')    #2D bin number                                                                    
i_x = B.get_data(f, 'i_x')    #x (th_nq) bin number                                                                                    
i_y = B.get_data(f, 'i_y')    #y (pmiss) bin number                                        
thnq = B.get_data(f, 'xb')      #th_nq value at bin center                                                          
pm =  B.get_data(f, 'yb')      #pmiss value at bin center   
pm_avg = pd.Series(B.get_data(f, 'pm_avg'))
# pm_avg = B.get_data(f, 'pm_avg')

# interpolate nan values to make curves smooth
# pm_avg = pm_avg.interpolate(method='cubicspline',limit=3)

#Get Combined Final Red Xsec for all kinematics
red_dataXsec_avg     = B.get_data(f,'red_dataXsec_avg')
red_dataXsec_avg_err = B.get_data(f,'red_dataXsec_avg_err')
red_dataXsec_avg_syst_err = B.get_data(f,'red_dataXsec_avg_syst_err')
red_dataXsec_avg_tot_err = B.get_data(f,'red_dataXsec_avg_tot_err')
red_pwiaXsec_avg     = pd.Series(B.get_data(f,'red_pwiaXsec_avg'))
red_fsiXsec_avg      = pd.Series(B.get_data(f,'red_fsiXsec_avg'))
# red_pwiaXsec_avg     = B.get_data(f,'red_pwiaXsec_avg')
# red_fsiXsec_avg      = B.get_data(f,'red_fsiXsec_avg')


#interpolate nan values to make theoretical curves smooth
# red_pwiaXsec_avg     = red_pwiaXsec_avg.interpolate(method='cubicspline',limit=3)
# red_fsiXsec_avg      = red_fsiXsec_avg.interpolate(method='cubicspline',limit=3)


#Get total relative errors / (Pm, thnq) bin for plotting (probably also write to table)
kin_syst_tot = dfile(fname)['kin_syst_tot']    #kinematic systematics
norm_syst_tot = dfile(fname)['norm_syst_tot']  #normalization systematics
tot_syst_err = dfile(fname)['tot_syst_err']    #total systematics
tot_stats_err = dfile(fname)['tot_stats_err']  #total statistical
tot_err = dfile(fname)['tot_err']              #overall error

'''
#----READ Q2 = (3,4) GeV2-----
fname2 = './Q2_3to4GeV/redXsec_combined.txt'
f2 = B.get_file(fname2)

#Get Bin Information (Same info for all files)                                                                                                  
i_b0 = B.get_data(f2, 'i_b')    #2D bin number                                                                    
i_x0 = B.get_data(f2, 'i_x')    #x (th_nq) bin number                                                                                    
i_y0 = B.get_data(f2, 'i_y')    #y (pmiss) bin number                                        
thnq0 = B.get_data(f2, 'xb')      #th_nq value at bin center                                                          
pm0 =  B.get_data(f2, 'yb')      #pmiss value at bin center   
pm_avg0 = B.get_data(f2, 'pm_avg')

#Get Combined Final Red Xsec for all kinematics
red_dataXsec_avg0     = B.get_data(f2,'red_dataXsec_avg')
red_dataXsec_avg_err0 = B.get_data(f2,'red_dataXsec_avg_err')
red_dataXsec_avg_syst_err0 = B.get_data(f2,'red_dataXsec_avg_syst_err')
red_dataXsec_avg_tot_err0 = B.get_data(f2,'red_dataXsec_avg_tot_err')
red_pwiaXsec_avg0     = B.get_data(f2,'red_pwiaXsec_avg')
red_fsiXsec_avg0      = B.get_data(f2,'red_fsiXsec_avg')

#Get total relative errors / (Pm, thnq) bin for plotting (probably also write to table)
kin_syst_tot0 = dfile(fname2)['kin_syst_tot']    #kinematic systematics
norm_syst_tot0 = dfile(fname2)['norm_syst_tot']  #normalization systematics
tot_syst_err0 = dfile(fname2)['tot_syst_err']    #total systematics
tot_stats_err0 = dfile(fname2)['tot_stats_err']  #total statistical
tot_err0 = dfile(fname2)['tot_err']              #overall error
'''

def plot_final():

    #----MASKING ARRAYS TO CHOOSE ERRORS BELOW 50 %-------
    #The Limit on the red. Xsec error should be placed at the very end, when the combined data is plotted (all data sets and overlapping bins have been combined)
    red_dataXsec_avg_masked = np.ma.array(red_dataXsec_avg, mask=(red_dataXsec_avg_err>0.*red_dataXsec_avg))
    red_dataXsec_avg_masked = np.ma.filled(red_dataXsec_avg_masked.astype(float), np.nan)
    
    #DATA Q2: (3,4) GeV2
    #red_dataXsec_avg_masked0 = np.ma.array(red_dataXsec_avg0, mask=(red_dataXsec_avg_err0>0.5*red_dataXsec_avg0))
    #red_dataXsec_avg_masked0 = np.ma.filled(red_dataXsec_avg_masked0.astype(float), np.nan)
    
    #Read Hall A data (reduced Xsec already in fm^3 units, and Precoil in GeV)
    pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35 = read_halla_data(35)
    pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45 = read_halla_data(45)  
    pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75 = read_halla_data(75)  

    #Read HallC 2017 data (red Xsec in fm^3 and Pm in GeV)
    pm_hc35, red_dataXsec_hc35, red_dataXsec_err_hc35 = read_hallc17_data(35)
    pm_hc45, red_dataXsec_hc45, red_dataXsec_err_hc45 = read_hallc17_data(45)
    pm_hc75, red_dataXsec_hc75, red_dataXsec_err_hc75 = read_hallc17_data(75)

    #Plot momentum distribution vs. Pmiss for different theta_nq
    thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

    for i, ithnq in enumerate(thnq_arr):

        th_nq_min = ithnq - 5
        th_nq_max = ithnq + 5

        #Create output files to write relative uncertainties for (Pm) bins for each th_nq setting
        fout_name = 'relative_errors_thnq%i.txt' % (ithnq)
        fout = open(fout_name, 'w')
        comment='#theta_nq = %i +/- 5 deg :: All errors are relative, dsig / sig [%%] \n' %(ithnq)
        header='#!pm_avg[f,0]/   kin_syst[f,1]/   norm_syst[f,2]/   tot_syst[f,3]/   tot_stats[f,4]/   tot_err[f,5]/\n'
        fout.write(comment)
        fout.write(header)

        #Read Other Theoretical Models (V18, CD-BONN)
        pm_avg1, red_pwiaXsec_V18 = read_theoretical_models("V18", "PWIA", ithnq)
        pm_avg2, red_fsiXsec_V18 = read_theoretical_models("V18", "FSI", ithnq) 
        pm_avg3, red_pwiaXsec_CD = read_theoretical_models("CD-Bonn", "PWIA", ithnq)                                   
        pm_avg4, red_fsiXsec_CD = read_theoretical_models("CD-Bonn", "FSI", ithnq)    

    
        # if(ithnq==35 or ithnq==45 or ithnq==75):
        #     #Read Hall A Data
        #     pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)


        #------Convert MeV^-3 to fm^3 for theory--------

        #THEORY
        red_pwiaXsec_avg[thnq==ithnq] = red_pwiaXsec_avg[thnq==ithnq]*MeV2fm    #Laget PWIA
        red_fsiXsec_avg[thnq==ithnq] = red_fsiXsec_avg[thnq==ithnq]*MeV2fm    #Laget FSI
        
        #Q2: (3,4) GeV2
        #red_pwiaXsec_avg0[thnq==ithnq] = red_pwiaXsec_avg0[thnq==ithnq]*MeV2fm    #Laget PWIA
        #red_fsiXsec_avg0[thnq==ithnq] = red_fsiXsec_avg0[thnq==ithnq]*MeV2fm    #Laget FSI

        red_pwiaXsec_V18 = red_pwiaXsec_V18*MeV2fm
        red_fsiXsec_V18 = red_fsiXsec_V18*MeV2fm
        red_pwiaXsec_CD = red_pwiaXsec_CD*MeV2fm
        red_fsiXsec_CD = red_fsiXsec_CD*MeV2fm


        #Do Interpolation to make theory curves smooth
        pmiss_avg = pm_avg[thnq==ithnq].interpolate(method='cubicspline',limit_area='inside')
        red_pwiaXsec_avg_ithnq = red_pwiaXsec_avg[thnq==ithnq].interpolate(method='cubicspline',limit_area='inside')
        red_fsiXsec_avg_ithnq = red_fsiXsec_avg[thnq==ithnq].interpolate(method='cubicspline',limit_area='inside')

        # f_red_pwiaXsec_avg = interp1d(pm_avg[sel1], red_pwiaXsec_avg[sel1],fill_value='extrapolate',bounds_error=False)   #Paris (J.M. Laget Calculation)
        # f_red_fsiXsec_avg = interp1d(pm_avg[sel1], red_fsiXsec_avg[sel1],fill_value='extrapolate',bounds_error=False)    
        f_red_pwiaXsec_avg = interp1d(pmiss_avg, red_pwiaXsec_avg_ithnq,fill_value='extrapolate',bounds_error=False,kind='cubic')   #Paris (J.M. Laget Calculation)
        f_red_fsiXsec_avg = interp1d(pmiss_avg, red_fsiXsec_avg_ithnq,fill_value='extrapolate',bounds_error=False,kind='cubic')    

        #Paris Potential for Q2 = (3,4) GeV2
        #f_red_pwiaXsec_avg0 = interp1d(pm_avg0[thnq==ithnq], red_pwiaXsec_avg0[thnq==ithnq],fill_value='extrapolate')   #Paris (J.M. Laget Calculation)
        #f_red_fsiXsec_avg0 = interp1d(pm_avg0[thnq==ithnq], red_fsiXsec_avg0[thnq==ithnq],fill_value='extrapolate')    

        f_red_pwiaXsec_V18 = interp1d(pm_avg1, red_pwiaXsec_V18,fill_value='extrapolate')                        #AV18 (M. Sargsian calculation)
        f_red_fsiXsec_V18 = interp1d(pm_avg2, red_fsiXsec_V18,fill_value='extrapolate')
       
        f_red_pwiaXsec_CD = interp1d(pm_avg3, red_pwiaXsec_CD,fill_value='extrapolate')                          #CD-Bonn (M. Sargsian calculation)
        f_red_fsiXsec_CD = interp1d(pm_avg4, red_fsiXsec_CD,fill_value='extrapolate')

        print(ithnq)


        B.pl.clf()
        B.pl.figure(figsize=(7.5,6.5))

        #--------------HALL A DATA---------------
        if(ithnq==35):
            #Plot Hall A Data for thnq==35 or 45 deg
            B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, label='Hall A Data')
        if(ithnq==45):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, label='Hall A Data')
        if(ithnq==75):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, label='Hall A Data')
        #----------------------------------------
        
        #--------------HALL C 2017 DATA---------------
        if(ithnq==35):
            #Plot Hall A Data for thnq==35 or 45 deg
            B.plot_exp(pm_hc35, red_dataXsec_hc35, red_dataXsec_err_hc35, marker='s',  markersize=5, color='c', markerfacecolor='white', capsize=0, label='Hall C 2017 Data')
        if(ithnq==45):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_hc45, red_dataXsec_hc45, red_dataXsec_err_hc45, marker='s',  markersize=5, color='c', markerfacecolor='white', capsize=0, label='Hall C 2017 Data')
        if(ithnq==75):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_hc75, red_dataXsec_hc75, red_dataXsec_err_hc75, marker='s',  markersize=5, color='c', markerfacecolor='white', capsize=0, label='Hall C 2017 Data')
        #----------------------------------------

        #Plot theoretical curves
        
        # B.plot_exp(pmiss_avg, f_red_pwiaXsec_avg(pmiss_avg), linestyle='--', marker='None', color='blue', logy=True, label='JML Paris PWIA')
        # B.plot_exp(pmiss_avg, f_red_fsiXsec_avg(pmiss_avg), linestyle='-', marker='None', color='blue', logy=True, label='JML Paris FSI')

        # B.plot_exp(pm_avg[thnq==ithnq], f_red_pwiaXsec_avg(pm_avg[thnq==ithnq]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
        # B.plot_exp(pm_avg[thnq==ithnq], f_red_fsiXsec_avg(pm_avg[thnq==ithnq]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
        
        # B.plot_exp(pm_avg[thnq==ithnq], red_pwiaXsec_avg[thnq==ithnq], linestyle='--', marker='None', color='blue', logy=True, label='JML Paris PWIA')
        # B.plot_exp(pm_avg[thnq==ithnq], red_fsiXsec_avg[thnq==ithnq], linestyle='-', marker='None', color='blue', logy=True, label='JML Paris FSI')
        B.plot_exp(pmiss_avg, red_pwiaXsec_avg_ithnq, linestyle='--', marker='None', color='blue', logy=True, label='JML Paris PWIA')
        B.plot_exp(pmiss_avg, red_fsiXsec_avg_ithnq, linestyle='-', marker='None', color='blue', logy=True, label='JML Paris FSI')


        #Q2 = (3,4) GeV2
        #B.plot_exp(pm_avg0[thnq==ithnq], f_red_pwiaXsec_avg0(pm_avg0[thnq==ithnq]), linestyle='--', marker='None', color='red', logy=True, label=r'Paris PWIA ($Q^{2}=3.5\pm0.5$ GeV$^{2}$)')
        #B.plot_exp(pm_avg0[thnq==ithnq], f_red_fsiXsec_avg0(pm_avg0[thnq==ithnq]), linestyle='-', marker='None', color='red', logy=True, label=r'Paris FSI ($Q^{2}=3.5\pm0.5$ GeV$^{2}$)')


        B.plot_exp(pm_avg1, f_red_pwiaXsec_V18(pm_avg1), linestyle='--', marker='None', color='green', logy=True, label='MS AV18 PWIA')   
        B.plot_exp(pm_avg2, f_red_fsiXsec_V18(pm_avg2), linestyle='-', marker='None', color='green', logy=True, label='MS AV18 FSI') 

        B.plot_exp(pm_avg3, f_red_pwiaXsec_CD(pm_avg3), linestyle='--', marker='None', color='magenta', logy=True, label='MS CD-Bonn PWIA')     
        B.plot_exp(pm_avg4, f_red_fsiXsec_CD(pm_avg4), linestyle='-', marker='None', color='magenta', logy=True, label='MS CD-Bonn FSI') 
        
        #----------PLOT FINAL REDUCED XSEC VS. THETA_nq------------
        
        B.plot_exp(pm[thnq==ithnq], red_dataXsec_avg[thnq==ithnq]*MeV2fm, red_dataXsec_avg_err[thnq==ithnq]*MeV2fm, marker='o', color='black', logy=True, label='New Data' )
        # B.plot_exp(pm[thnq==ithnq], red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm, red_dataXsec_avg_err[thnq==ithnq]*MeV2fm, marker='o', color='r', markerfacecolor='none', label='Statistical Error' )
        # B.plot_exp(pm_avg[thnq==ithnq], red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm, red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm, marker='o', markersize=6, color='k', markerfacecolor='k', label='This Experiment (Hall C)' )
                
        #B.plot_exp(pm_avg0[thnq==ithnq], red_dataXsec_avg_masked0[thnq==ithnq]*MeV2fm, red_dataXsec_avg_tot_err0[thnq==ithnq]*MeV2fm, marker='o', markersize=6, color='r', markerfacecolor='r', label=r'Hall C, $Q^{2}=3.5\pm0.5$ GeV$^{2}$' )
        if ithnq == 35:
            B.pl.text(0.1, 1e-4, 'PRELIMINARY', color='gray',alpha=0.4,
                      fontsize=54)
        if ithnq == 45:
            B.pl.text(0.1, 1e-4, 'PRELIMINARY', color='gray',alpha=0.4,
                      fontsize=54)
        if ithnq == 75:
            B.pl.text(0.1, 1e-3, 'PRELIMINARY', color='gray',alpha=0.4,
                      fontsize=54)

        B.pl.xlabel('$p_{r}$ (GeV/c)', fontsize=22)
        B.pl.ylabel(r'$\sigma_{red}\, (\textrm{fm}^{3})$', fontsize=22)
        #B.pl.ylim(0, 1.e-4)
        B.pl.xlim(0., 1.2)
        B.pl.title(r'Reduced Cross Section, $\theta_{nq} = %i \pm 5$ deg '%(ithnq))
        #B.pl.yscale('linear')
        B.pl.legend(frameon=False, fontsize='x-small')            
        B.pl.savefig(dir_name+'/final_redXsec_thnq%i.pdf'%(ithnq))        

        #------Plot the Exp. Data relative errors dsig / sig  to know how large are the contributions from systematics and statistical------ 
        
        B.pl.clf()
        B.pl.figure(i+1)

        y = np.array([0. for ii in range(len(pm_avg[thnq==ithnq]))])

        dsig_sig_stats = red_dataXsec_avg_err[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)
        dsig_sig_syst = red_dataXsec_avg_syst_err[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)
        dsig_sig_tot = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)
        
        #Q2:(3,4)
        #dsig_sig_stats = red_dataXsec_avg_err0[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked0[thnq==ithnq]*MeV2fm)
        #dsig_sig_syst = red_dataXsec_avg_syst_err0[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked0[thnq==ithnq]*MeV2fm)
        #dsig_sig_tot = red_dataXsec_avg_tot_err0[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked0[thnq==ithnq]*MeV2fm)


        for iii in range(len(y)):
            if (np.isnan(dsig_sig_tot[iii])):
                y[iii] = np.nan

        B.plot_exp(pm_avg[thnq==ithnq],  y, dsig_sig_stats*100., marker='o', markersize=6, color='b', label='Statistical Error')
        B.plot_exp(pm_avg[thnq==ithnq],  y, dsig_sig_syst*100., marker='o', markersize=6, color='r', label='Systematics Error')
        B.plot_exp(pm_avg[thnq==ithnq],  y, dsig_sig_tot*100., marker='o', markersize=6, color='k', label='Total Error')

        B.pl.xlabel(r' $p_{r}$ (GeV/c)')
        B.pl.ylabel(r'Relative Error (\%)')
        B.pl.ylim(-100, 100)
        B.pl.title(r'Data Relative Error, $\theta_{nq} = %i \pm 5$ deg'%(ithnq))
        B.pl.legend()            
        B.pl.savefig(dir_name+'/final_redXsec_relativeError_thnq%i.pdf'%(ithnq))
        
        #------Plot the Systematic Relative Errors dsig / sig  to know how large are the contributions from (kin, norm)systematics------ 
        
        B.pl.clf()
        B.pl.figure(i+2)

        y = np.array([0. for i in range(len(pm_avg[thnq==ithnq]))])

        dsig_sig_kin_syst = kin_syst_tot[thnq==ithnq]*100.
        dsig_sig_norm_syst = norm_syst_tot[thnq==ithnq]*100.
        dsig_sig_tot_syst = tot_syst_err[thnq==ithnq]*100.
        rel_stats_err = tot_stats_err[thnq==ithnq]*100.
        rel_tot_err = tot_err[thnq==ithnq]*100.
        pmiss_avg = pm_avg[thnq==ithnq].reset_index(drop=True)
        # pmiss_avg = pm_avg[thnq==ithnq]
        
        #Q2;(3,4) GeV2
        #dsig_sig_kin_syst = kin_syst_tot0[thnq==ithnq]*100.
        #dsig_sig_norm_syst = norm_syst_tot0[thnq==ithnq]*100.
        #dsig_sig_tot_syst = tot_syst_err0[thnq==ithnq]*100.
        #rel_stats_err = tot_stats_err0[thnq==ithnq]*100.
        #rel_tot_err = tot_err0[thnq==ithnq]*100.
        #pmiss_avg = pm_avg0[thnq==ithnq]

        for i in range(len(y)):
            if (np.isnan(dsig_sig_tot[i])):
                y[i] = np.nan

        B.plot_exp(pm_avg[thnq==ithnq],  y, dsig_sig_kin_syst, marker='o', markersize=6, color='r', label='Kinematics Systematic Error')
        B.plot_exp(pm_avg[thnq==ithnq],  y, dsig_sig_norm_syst, marker='o', markersize=6, color='b', label='Normalization Systematic Error')
        B.plot_exp(pm_avg[thnq==ithnq],  y, dsig_sig_tot_syst, marker='o', markersize=6, color='k', label='Total Systematic Error')

        B.pl.xlabel(r' $p_{r}$ (GeV/c)')
        B.pl.ylabel(r'Relative Error (\%)')
        B.pl.ylim(-25, 25)
        B.pl.title(r'Systematic Relative Error, $\theta_{nq} = %i \pm 5$ deg'%(ithnq))
        B.pl.legend()            
        B.pl.savefig(dir_name+'/final_redXsec_SystrelativeError_thnq%i.png'%(ithnq))
        
        #Write Relative Errors to file
        #print('====LENGTHS====')
        # print('pm_avg >>>>>> = ',pmiss_avg)

        #np.savetxt(fout_name,np.c_[pmiss_avg, dsig_sig_kin_syst, dsig_sig_norm_syst, dsig_sig_tot_syst,  dsig_sig_stats, dsig_sig_tot], delimiter='  ', fmt='%.4f')
        for ith in range(len(pmiss_avg)):
            print('ith=',ith,' :: pmiss_avg=',pmiss_avg[ith])
            fout.write('%.5f  %.5f   %.5f   %.5f  %.5f   %.5f\n' % (pmiss_avg[ith], dsig_sig_kin_syst[ith], dsig_sig_norm_syst[ith], dsig_sig_tot_syst[ith], rel_stats_err[ith], rel_tot_err[ith]))

        fout.close()

        '''
        #--------------Plot the % Deviation of data and all models relative to the CD-Bonn PWIA model------------------------
        
        B.pl.clf()
        B.pl.figure(i+3)
        
        #Calculate percent deviation from data and all models relative

        #print('red_dataXsec=',(red_dataXsec_avg_masked[thnq==ithnq]))

        f_red_pwiaXsec_CD_arr = f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        ref_line = (f_red_pwiaXsec_CD_arr - f_red_pwiaXsec_CD_arr)/f_red_fsiXsec_CD_arr * 100.

        #print('f_red_CD=',(f_red_fsiXsec_CD_arr))
        
        data_dev = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_pwiaXsec_CD_arr) / (f_red_pwiaXsec_CD_arr) * 100.
        data_dev_err = np.sqrt((red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm/f_red_pwiaXsec_CD_arr )**2) * 100.
        #print('red_dataXsec=',red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)
        #print('f_red_pwiaXsec=',f_red_pwiaXsec_CD_arr)
        #print('data_dev=',data_dev)
        #print('data_dev_err=',data_dev_err)
        #print('len of pm=',len(pm[thnq==ithnq]))
        #print('len of pm_avg=',len(pm_avg[thnq==ithnq]))
        CD_fsi_dev = (f_red_fsiXsec_CD(pm_avg[thnq==ithnq]) - f_red_pwiaXsec_CD_arr)/f_red_pwiaXsec_CD_arr * 100.
        
        V18_fsi_dev = (f_red_fsiXsec_V18(pm_avg[thnq==ithnq]) - f_red_pwiaXsec_CD_arr)/f_red_pwiaXsec_CD_arr * 100.
        V18_pwia_dev = (f_red_pwiaXsec_V18(pm_avg[thnq==ithnq]) - f_red_pwiaXsec_CD_arr)/f_red_pwiaXsec_CD_arr * 100.
        
        paris_fsi_dev = (f_red_fsiXsec_avg(pm_avg[thnq==ithnq]) - f_red_pwiaXsec_CD_arr)/f_red_pwiaXsec_CD_arr * 100.
        paris_pwia_dev = (f_red_pwiaXsec_avg(pm_avg[thnq==ithnq]) - f_red_pwiaXsec_CD_arr)/f_red_pwiaXsec_CD_arr * 100.

        B.plot_exp(pm_avg[thnq==ithnq],  ref_line, marker='None', linestyle='-', color='magenta', label='CD-Bonn FSI (reference)')
        B.plot_exp(pm_avg[thnq==ithnq],  data_dev, data_dev_err, marker='o', color='k', label='Data')
        B.plot_exp(pm_avg[thnq==ithnq],  CD_pwia_dev, marker='None', linestyle='--', color='magenta', label='CD-Bonn PWIA')
        B.plot_exp(pm_avg[thnq==ithnq],  V18_fsi_dev, marker='None', linestyle='-', color='g', label='V18 FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  V18_pwia_dev, marker='None', linestyle='--', color='g', label='V18 PWIA')
        B.plot_exp(pm_avg[thnq==ithnq],  paris_fsi_dev, marker='None', linestyle='-', color='b', label='Paris FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  paris_pwia_dev, marker='None', linestyle='--', color='b', label='Paris PWIA')


        B.pl.xlabel('Neutron Recoil Momenta [GeV/c]')
        B.pl.ylabel(r'Relative Error [\%]')
        B.pl.ylim(-400, 300)
        B.pl.title(r'Relative Error, $\theta_{nq} = %i \pm 5$ deg'%(ithnq))
        B.pl.legend()            
        B.pl.savefig(dir_name+'/final_redXsec_relativeErrorModel_thnq%i.pdf'%(ithnq))
        
        #---------Plot Ratio  sig_data (or other models) / sig_model_CDBonn_FSI, 
        #-----OR, plot % deviation of data to all models
        B.pl.clf()
        B.pl.figure(i+4)
        
        #Calculate percent deviation from data and all models relative

      
        #Take (Data - Model) / Data  ratio  (a - b) / a,  sig = b/a**2 *sig_a
        f_red_fsiXsec_CD_arr = f_red_fsiXsec_CD(pm_avg[thnq==ithnq])
        ref_line = (red_dataXsec_avg_masked[thnq==ithnq] - red_dataXsec_avg_masked[thnq==ithnq]) / red_dataXsec_avg_masked[thnq==ithnq] 
        
        fsiCD_ratio = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_fsiXsec_CD_arr) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)  * 100.
        fsiCD_ratio_err = f_red_fsiXsec_CD_arr / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)**2  * (red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm) * 100.

        pwiaCD_ratio = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_pwiaXsec_CD(pm_avg[thnq==ithnq]) ) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm) * 100.                  
        pwiaCD_ratio_err = f_red_pwiaXsec_CD(pm_avg[thnq==ithnq]) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)**2  * (red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm) * 100.
        
        V18_fsi_ratio = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_fsiXsec_V18(pm_avg[thnq==ithnq])) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm) * 100.
        V18_fsi_ratio_err = f_red_fsiXsec_V18(pm_avg[thnq==ithnq]) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)**2  * (red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm) * 100.

        V18_pwia_ratio = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_pwiaXsec_V18(pm_avg[thnq==ithnq])) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm) * 100.
        V18_pwia_ratio_err = f_red_pwiaXsec_V18(pm_avg[thnq==ithnq]) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)**2  * (red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm) * 100.

        paris_fsi_ratio = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_fsiXsec_avg(pm_avg[thnq==ithnq])) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm) * 100.
        paris_fsi_ratio_err = f_red_fsiXsec_avg(pm_avg[thnq==ithnq]) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)**2  * (red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm) * 100.

        paris_pwia_ratio = (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm - f_red_pwiaXsec_avg(pm_avg[thnq==ithnq])) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm) * 100.
        paris_pwia_ratio_err = f_red_pwiaXsec_avg(pm_avg[thnq==ithnq]) / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm)**2  * (red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm)* 100. 
        
        B.plot_exp(pm_avg[thnq==ithnq],  ref_line, marker='None', linestyle='-', color='k', label='Data (reference)')
        B.plot_exp(pm_avg[thnq==ithnq],  fsiCD_ratio, fsiCD_ratio_err, marker='o', ms=6, color='magenta', label='CD-Bonn FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  pwiaCD_ratio,pwiaCD_ratio_err, marker='o', ms=6, ecolor='magenta', mec='magenta', mfc='white', label='CD-Bonn PWIA')
        
        B.plot_exp(pm_avg[thnq==ithnq],  V18_fsi_ratio, V18_fsi_ratio_err,  marker='s', ms=6, color='g', label='V18 FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  V18_pwia_ratio, V18_pwia_ratio_err, marker='s', ms=6, ecolor='g', mec='g', mfc='white', label='V18 PWIA') 

        B.plot_exp(pm_avg[thnq==ithnq],  paris_fsi_ratio, paris_fsi_ratio_err, marker='D', ms=6, color='b', label='Paris FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  paris_pwia_ratio,paris_pwia_ratio_err, marker='D', ms=6, ecolor='b', mec='b', mfc='white', label='Paris PWIA')

        B.pl.xlabel('Neutron Recoil Momenta [GeV/c]')
        B.pl.ylabel(r'$(\sigma^{Data}_{red} - \sigma^{Model}_{red}) /\sigma^{Data}_{red} [\%]$')
        B.pl.ylim(-200, 200)
        B.pl.title(r'Percent Deviation of Model from Data, $\theta_{nq} = %i \pm 5$ deg'%(ithnq))
        B.pl.legend(loc='upper right')            
        B.pl.savefig(dir_name+'/final_redXsec_Data2ModelDev_thnq%i.pdf'%(ithnq))
        
        
        #Take Data / Model ratio
        f_red_fsiXsec_CD_arr = f_red_fsiXsec_CD(pm_avg[thnq==ithnq])
        ref_line = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / (red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm) 
        
        fsiCD_ratio = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_fsiXsec_CD_arr
        fsiCD_ratio_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / f_red_fsiXsec_CD_arr 

        pwiaCD_ratio = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])                   
        pwiaCD_ratio_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm/f_red_pwiaXsec_CD(pm_avg[thnq==ithnq]) 
        
        V18_fsi_ratio = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_fsiXsec_V18(pm_avg[thnq==ithnq])
        V18_fsi_ratio_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / f_red_fsiXsec_V18(pm_avg[thnq==ithnq]) 

        V18_pwia_ratio = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_V18(pm_avg[thnq==ithnq]) 
        V18_pwia_ratio_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_V18(pm_avg[thnq==ithnq])

        paris_fsi_ratio = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_fsiXsec_avg(pm_avg[thnq==ithnq])
        paris_fsi_ratio_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / f_red_fsiXsec_avg(pm_avg[thnq==ithnq]) 

        paris_pwia_ratio = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_avg(pm_avg[thnq==ithnq])
        paris_pwia_ratio_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_avg(pm_avg[thnq==ithnq])  

        B.plot_exp(pm_avg[thnq==ithnq],  ref_line, marker='None', linestyle='-', color='k', label='Data (reference)')
        B.plot_exp(pm_avg[thnq==ithnq],  fsiCD_ratio, fsiCD_ratio_err, marker='o', ms=4, color='magenta', label='CD-Bonn FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  pwiaCD_ratio,pwiaCD_ratio_err, marker='o', ms=4, ecolor='magenta', mec='magenta', mfc='white', label='CD-Bonn PWIA')
        
        B.plot_exp(pm_avg[thnq==ithnq],  V18_fsi_ratio, V18_fsi_ratio_err,  marker='s', ms=4, color='g', label='V18 FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  V18_pwia_ratio, V18_pwia_ratio_err, marker='s', ms=4, ecolor='g', mec='g', mfc='white', label='V18 PWIA') 

        B.plot_exp(pm_avg[thnq==ithnq],  paris_fsi_ratio, paris_fsi_ratio_err, marker='^', ms=4, color='b', label='Paris FSI')
        B.plot_exp(pm_avg[thnq==ithnq],  paris_pwia_ratio,paris_pwia_ratio_err, marker='^', ms=4, ecolor='b', mec='b', mfc='white', label='Paris PWIA')

        B.pl.xlabel('Neutron Recoil Momenta [GeV/c]')
        B.pl.ylabel(r'Ratio $\sigma^{Data}_{red}/\sigma^{Model}_{red}$')
        B.pl.ylim(-1.0, 6)
        B.pl.title(r'Ratio $\sigma^{Data}_{red}/\sigma^{Model}_{red}$, $\theta_{nq} = %i \pm 5$ deg'%(ithnq))
        B.pl.legend()            
        B.pl.savefig(dir_name+'/final_redXsec_Data2Model_thnq%i.pdf'%(ithnq))
       

        #-------Plot the Ratio  sig_red_exp(pm) / sig_red_exp(p0=0.5 GeV/c) for pm >=0.5 GeV/c (same for models), to compare shapes
        
                
        #Require ONLY thnq = 35, 45 deg (FOCUS POINT: pm= 0.5 GeV)

        if (ithnq==35 or ithnq==45):
                            
            fp = np.array([0.5, 0.82])
            for j in range(len(fp)):

                B.pl.clf()                                                                                                          
                B.pl.figure(j+4) 

                print('fp=',fp[j])
                #DATA
                print('pm=',pm[thnq==ithnq])
                sig_exp =  red_dataXsec_avg_masked[thnq==ithnq] 
                sig_exp_p0 = red_dataXsec_avg_masked[(thnq==ithnq)][np.where(pm[thnq==ithnq]==fp[j])]
                sig_exp_err =  red_dataXsec_avg_tot_err[thnq==ithnq] 
                sig_exp_p0_err = red_dataXsec_avg_tot_err[thnq==ithnq][np.where(pm[thnq==ithnq]==fp[j])]

                print('sig_exp=',sig_exp)
                print('sig_exp_p0=', sig_exp_p0)
                #Paris
                sig_paris_pwia = f_red_pwiaXsec_avg(pm[thnq==ithnq])
                sig_paris_pwia_p0 = f_red_pwiaXsec_avg(pm[thnq==ithnq][np.where(pm[thnq==ithnq]==fp[j])]) 
                
                sig_paris_fsi = f_red_fsiXsec_avg(pm[thnq==ithnq])                                                                                                            
                sig_paris_fsi_p0 = f_red_fsiXsec_avg(pm[thnq==ithnq][np.where(pm[thnq==ithnq]==fp[j])]) 
                
                
                #AV18                                                                                                                   
                sig_V18_pwia = f_red_pwiaXsec_V18(pm[thnq==ithnq])        #V18 Xsec                                                                             
                sig_V18_pwia_p0 = f_red_pwiaXsec_V18(pm[np.where(pm[thnq==ithnq]==fp[j])])  #V18 Xsec given that Pm = 0.5, 0.8 GeV 
                
                sig_V18_fsi = f_red_fsiXsec_V18(pm[thnq==ithnq])        #V18 Xsec                                                                      
                sig_V18_fsi_p0 = f_red_fsiXsec_V18(pm[np.where(pm[thnq==ithnq]==fp[j])])  #V18 Xsec given that Pm = 0.5, 0.8 GeV   
                
                #CD-Bonn
                sig_CD_pwia = f_red_pwiaXsec_CD(pm[thnq==ithnq])        #CD-Bonn Xsec
                sig_CD_pwia_p0 = f_red_pwiaXsec_CD(pm[np.where(pm[thnq==ithnq]==fp[j])])  #CD-Bonn Xsec given that Pm = 0.5, 0.8 GeV 
                
                sig_CD_fsi = f_red_fsiXsec_CD(pm[thnq==ithnq])        #CD-Bonn Xsec                                                                     
                sig_CD_fsi_p0 = f_red_fsiXsec_CD(pm[np.where(pm[thnq==ithnq]==fp[j])])  #CD-Bonn Xsec given that Pm = 0.5, 0.8 GeV  

                #print('pm=',pm_bin3)
                #print('sig_CD=',sig_CD_pwia)
                #print('sig_CD_p0=',sig_CD_pwia_p0)

                #DEFINE THE RATIOS
                #DATA
                R_exp = sig_exp / float(sig_exp_p0)
                R_exp_err = np.sqrt(((1. / sig_exp_p0)**2 *  sig_exp_err**2 )  + (sig_exp / sig_exp_p0**2)**2 * sig_exp_p0_err**2)
                #Paris
                R_paris_pwia = sig_paris_pwia / sig_paris_pwia_p0
                R_paris_fsi = sig_paris_fsi / sig_paris_fsi_p0 
                #AV18
                R_V18_pwia = sig_V18_pwia / sig_V18_pwia_p0                                                                
                R_V18_fsi = sig_V18_fsi / sig_V18_fsi_p0 
                #CD-Bonn                                                                                                                          
                R_CD_pwia = sig_CD_pwia / sig_CD_pwia_p0                                                                                            
                R_CD_fsi = sig_CD_fsi / sig_CD_fsi_p0 
                
                B.plot_exp(pm[thnq==ithnq], R_exp, R_exp_err, marker='o', color='k', label='DATA', logy=True)
                
                B.plot_exp(pm[thnq==ithnq], R_paris_pwia, linestyle='--',  marker='None', color='blue', label='Paris PWIA', logy=True)  
                B.plot_exp(pm[thnq==ithnq], R_paris_fsi, linestyle='-',  marker='None', color='blue', label='Paris FSI', logy=True) 
                
                B.plot_exp(pm[thnq==ithnq], R_V18_pwia, linestyle='--',  marker='None', color='green', label='V18 PWIA', logy=True)                             
                B.plot_exp(pm[thnq==ithnq], R_V18_fsi, linestyle='-',  marker='None', color='green', label='V18 FSI', logy=True)
                
                B.plot_exp(pm[thnq==ithnq], R_CD_pwia, linestyle='--',  marker='None', color='magenta', label='CD-Bonn PWIA', logy=True)                                         
                B.plot_exp(pm[thnq==ithnq], R_CD_fsi, linestyle='-',  marker='None', color='magenta', label='CD-Bonn FSI', logy=True) 
                
                #xdata = pm[thnq==ithnq][np.where(pm[thnq==ithnq]>=fp[j])]
                #ydata =  R_exp[np.where(pm[thnq==ithnq]>=fp[j])]   
                #ydata_err =  R_exp_err[np.where(pm[thnq==ithnq]>=fp[j])]
                if(ithnq==35):
                    xdata = pm[thnq==ithnq]
                    ydata = R_exp
                    ydata_err = R_exp_err

                    xd = xdata[(~np.isnan(ydata)) & (xdata>=0.5)]
                    yd = ydata[(~np.isnan(ydata)) & (xdata>=0.5)]
                    yd_err = ydata_err[(~np.isnan(ydata)) & (xdata>=0.5)]


                if(ithnq==45):
                    xdata=pm[thnq==ithnq]
                    ydata=R_exp
                    ydata_err = R_exp_err
               
                    xd = xdata[(~np.isnan(ydata)) & (xdata>=0.5)]
                    yd = ydata[(~np.isnan(ydata)) & (xdata>=0.5)]
                    yd_err = ydata_err[(~np.isnan(ydata)) & (xdata>=0.5)]
                
                #Define Fit Function
                def f(x):
                    y = b() * np.exp(m()*x)    #lny = m*x + ln b
                    return y

                m = B.Parameter(1., 'm')
                b = B.Parameter(1., 'b')

                F = B.genfit(f,[m,b], xd, yd, yd_err)

                mp = m.get()[1] ; mp_err=m.get()[2]
                bp = np.log(b.get()[1]) ; bp_err=np.log(b.get()[2])

                B.plot_line(F.xpl, F.ypl, color='r', lw=2, label='FIT \n slope: %.3f $\pm$ %.3f \n y-int.= %.3f $\pm$ %.3f'%(mp, mp_err, bp, bp_err))
                B.pl.yscale('log')
                B.pl.xlim(0., 1.19)
                B.pl.ylim(0, 1e6)
                
                B.pl.xlabel(r'Neutron Recoil Momenta, $p_{r}$ [GeV/c]')                                                                                                                                                            
                B.pl.ylabel(r'Cross Section Ratio')                                                                                                                                                                       
                B.pl.title(r'Cross Section Ratio, $\theta_{nq} = %i \pm 5$ deg'%(ithnq))                                                                                                                             
                B.pl.legend(loc='upper right', fontsize='small')                                                                                                                                                                                            
                B.pl.savefig(dir_name+'/ratio_test_fp%f_thnq%i.pdf'%(fp[j], ithnq))   
        '''

    B.pl.figure()
    R_560 = red_dataXsec_avg[pm==0.56]/
        
        
        


def read_halla_data(thnq=0):

    #The code read data files containing momentum distributions from Hall A

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    fminv2GeV = 1./(1./0.1973)
    ifm2GeV = 0.1973  #inverse fermi to GeV

    thnq_f = '%.1f' % (thnq)
    fname = 'yield_n_xsec/HallA_data/old_hallA/theta_%s_th_corr.data' % (thnq_f)

    kin = dfile(fname)
    p_miss = kin['p_miss'] * ifm2GeV    #Pm here is not the central bin value, but an average over that bin
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
    
def read_theoretical_models(theory="", model="", thnq=0):

    #This code read the averaged red.Xsec and returns arrays in pm and reduced Xsec
    #theory: V18, CD-Bonn    model: PWIA, FSI
    thnq_f = "%.2f" %(thnq)

    fname = '/home/gvill/deuteron/DEUTERON_ANALYSIS/ANALYSIS_SCRIPTS/MAIN_ANALYSIS/Deep_CrossSections/theoretical_models'+\
    '/updated_%s_models/theoryXsec_%s%s_thnq%s_combined.data' % (model, theory, model, thnq_f)
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
    red_pwiaXsec_V18_35 = red_pwiaXsec_V18_35*MeV2fm
    red_fsiXsec_V18_35 = red_fsiXsec_V18_35*MeV2fm
    red_pwiaXsec_CD_35 = red_pwiaXsec_CD_35*MeV2fm
    red_fsiXsec_CD_35 = red_fsiXsec_CD_35*MeV2fm

    red_pwiaXsec_avg_45 = red_pwiaXsec_avg[thnq==45]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_45 = red_fsiXsec_avg[thnq==45]*MeV2fm      #Laget FSI
    red_pwiaXsec_V18_45 = red_pwiaXsec_V18_45*MeV2fm
    red_fsiXsec_V18_45 = red_fsiXsec_V18_45*MeV2fm
    red_pwiaXsec_CD_45 = red_pwiaXsec_CD_45*MeV2fm
    red_fsiXsec_CD_45 = red_fsiXsec_CD_45*MeV2fm

    red_pwiaXsec_avg_75 = red_pwiaXsec_avg[thnq==75]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_75 = red_fsiXsec_avg[thnq==75]*MeV2fm      #Laget FSI
    red_pwiaXsec_V18_75 = red_pwiaXsec_V18_75*MeV2fm
    red_fsiXsec_V18_75 = red_fsiXsec_V18_75*MeV2fm
    red_pwiaXsec_CD_75 = red_pwiaXsec_CD_75*MeV2fm
    red_fsiXsec_CD_75 = red_fsiXsec_CD_75*MeV2fm

    #Do Interpolation to make theory curves smooth
    #35 deg
    f_red_pwiaXsec_avg_35 = interp1d(pm_avg[thnq==35], red_pwiaXsec_avg_35,fill_value='extrapolate')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_35 = interp1d(pm_avg[thnq==35], red_fsiXsec_avg_35,fill_value='extrapolate')    
    
    f_red_pwiaXsec_V18_35 = interp1d(pm_avg1, red_pwiaXsec_V18_35,fill_value='extrapolate')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_35 = interp1d(pm_avg2, red_fsiXsec_V18_35,fill_value='extrapolate')
    
    f_red_pwiaXsec_CD_35 = interp1d(pm_avg3, red_pwiaXsec_CD_35,fill_value='extrapolate')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_35 = interp1d(pm_avg4, red_fsiXsec_CD_35,fill_value='extrapolate')
    #45 eg
    f_red_pwiaXsec_avg_45 = interp1d(pm_avg[thnq==45], red_pwiaXsec_avg_45,fill_value='extrapolate')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_45 = interp1d(pm_avg[thnq==45], red_fsiXsec_avg_45,fill_value='extrapolate')    
    
    f_red_pwiaXsec_V18_45 = interp1d(pm_avg5, red_pwiaXsec_V18_45,fill_value='extrapolate')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_45 = interp1d(pm_avg6, red_fsiXsec_V18_45,fill_value='extrapolate')
    
    f_red_pwiaXsec_CD_45 = interp1d(pm_avg7, red_pwiaXsec_CD_45,fill_value='extrapolate')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_45 = interp1d(pm_avg8, red_fsiXsec_CD_45,fill_value='extrapolate')
    #75 eg
    f_red_pwiaXsec_avg_75 = interp1d(pm_avg[thnq==75], red_pwiaXsec_avg_75,fill_value='extrapolate')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_75 = interp1d(pm_avg[thnq==75], red_fsiXsec_avg_75,fill_value='extrapolate')    
    
    f_red_pwiaXsec_V18_75 = interp1d(pm_avg9, red_pwiaXsec_V18_75,fill_value='extrapolate')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_75 = interp1d(pm_avg10, red_fsiXsec_V18_75,fill_value='extrapolate')
    
    f_red_pwiaXsec_CD_75 = interp1d(pm_avg11, red_pwiaXsec_CD_75,fill_value='extrapolate')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_75 = interp1d(pm_avg12, red_fsiXsec_CD_75,fill_value='extrapolate')
    
    
    #====THETA_NQ = 35 DEG====
    ax1 = B.pl.subplot(131)
    B.plot_exp(pm_avg[thnq==35], red_dataXsec_avg_masked[thnq==35]*MeV2fm, red_dataXsec_avg_tot_err[thnq==35]*MeV2fm, marker='o', markersize=5, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='o', markersize=5, color='c', logy=True,  label='Hall A Data')
    
    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==35], f_red_pwiaXsec_avg_35(pm_avg[thnq==35]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==35], f_red_fsiXsec_avg_35(pm_avg[thnq==35]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg1, f_red_pwiaXsec_V18_35(pm_avg1), linestyle='--', marker='None', color='green', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg2, f_red_fsiXsec_V18_35(pm_avg2), linestyle='-', marker='None', color='green', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg3, f_red_pwiaXsec_CD_35(pm_avg3), linestyle='--', marker='None', color='magenta', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg4, f_red_fsiXsec_CD_35(pm_avg4), linestyle='-', marker='None', color='magenta', logy=True, label='CD-Bonn FSI') 

    #B.pl.xlabel('$p_{r}$ (GeV/c)')
    B.pl.title('')
    B.pl.xlabel('')
    B.pl.ylabel(r'$\sigma_{red} (fm^{3})$', fontsize=15)
    B.pl.subplots_adjust(wspace=.001)
    B.pl.xlim(0., 1.2)
    B.pl.xticks(np.arange(0.1, 1.3, 0.2))
    B.pl.ylim(0., 50)
    

    #====THETA_NQ = 45 DEG====
    B.pl.subplot(132)
    B.plot_exp(pm_avg[thnq==45], red_dataXsec_avg_masked[thnq==45]*MeV2fm, red_dataXsec_avg_tot_err[thnq==45]*MeV2fm, marker='o', markersize=5, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='o', markersize=5, color='c', logy=True, label='Hall A Data')
    
    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==45], f_red_pwiaXsec_avg_45(pm_avg[thnq==45]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==45], f_red_fsiXsec_avg_45(pm_avg[thnq==45]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg5, f_red_pwiaXsec_V18_45(pm_avg5), linestyle='--', marker='None', color='green', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg6, f_red_fsiXsec_V18_45(pm_avg6), linestyle='-', marker='None', color='green', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg7, f_red_pwiaXsec_CD_45(pm_avg7), linestyle='--', marker='None', color='magenta', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg8, f_red_fsiXsec_CD_45(pm_avg8), linestyle='-', marker='None', color='magenta', logy=True, label='CD-Bonn FSI') 

    B.pl.gca().set_yticklabels([''])
    B.pl.title('')
    B.pl.ylabel('')
    B.pl.xlabel('$p_{r}$ (GeV/c)', fontsize=15)
    B.pl.subplots_adjust(wspace=.001)
    B.pl.xlim(0., 1.2)
    B.pl.xticks(np.arange(0.1, 1.3, 0.2))
    B.pl.ylim(0., 50)
    
    #====THETA_NQ = 75 DEG====
    B.pl.subplot(133)
    B.plot_exp(pm_avg[thnq==75], red_dataXsec_avg_masked[thnq==75]*MeV2fm, red_dataXsec_avg_tot_err[thnq==75]*MeV2fm, marker='o', markersize=5, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.plot_exp(0.,0.,0., marker='o', markersize=5, color='c',  label='Hall A Data')
    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==75], f_red_pwiaXsec_avg_75(pm_avg[thnq==75]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==75], f_red_fsiXsec_avg_75(pm_avg[thnq==75]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg9, f_red_pwiaXsec_V18_75(pm_avg9), linestyle='--', marker='None', color='green', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg10, f_red_fsiXsec_V18_75(pm_avg10), linestyle='-', marker='None', color='green', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg11, f_red_pwiaXsec_CD_75(pm_avg11), linestyle='--', marker='None', color='magenta', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg12, f_red_fsiXsec_CD_75(pm_avg12), linestyle='-', marker='None', color='magenta', logy=True, label='CD-Bonn FSI') 

    B.pl.gca().set_yticklabels([''])
    B.pl.subplots_adjust(wspace=.001)
    B.pl.title('')
    B.pl.ylabel('')
    B.pl.xlabel('')
    B.pl.xlim(0., 1.2)
    B.pl.xticks(np.arange(0.1, 1.3, 0.2))
    B.pl.ylim(0., 50)
    B.pl.legend(frameon=False, loc='upper right')   
    B.pl.show()

def plot_report():

    #-------Plot the report file quantities: trigger rates, efficiencies (detector, live time, etc.), 


    def get_fname(pm_set, data_set):
        fname='../../root_files/pm%i_fsiXsec_set%i_Em_final40MeV/report_deep_pm%i_set%i.dat'%(pm_set, data_set, pm_set, data_set)
        return fname

    def get_var(pm_set=0, data_set=0, var=''):
        fvar = np.array(dfile(get_fname(pm_set,data_set))[var])
        return fvar
         
    def get_heep_var(var=''):
        fname='../../root_files/HEEP_ELASTICS/report_heep.dat'
        fvar = np.array(dfile(fname)[var])
        return fvar

    
    #---------------Plot Accepted COin. Triggers Counts per charge---------------
    CPQ80_set1 = get_var(80,1,'ptrig6_accp')/get_var(80,1,'charge') ; CPQ80_set1_err = get_var(80,1,'ptrig6_accp')/get_var(80,1,'charge')**2 * (0.02*get_var(80,1,'charge'))
    CPQ580_set1 = get_var(580,1,'ptrig6_accp')/get_var(580,1,'charge') ; CPQ580_set1_err = get_var(580,1,'ptrig6_accp')/get_var(580,1,'charge')**2 * (0.02*get_var(580,1,'charge'))
    CPQ580_set2 = get_var(580,2,'ptrig6_accp')/get_var(580,2,'charge') ; CPQ580_set2_err = get_var(580,2,'ptrig6_accp')/get_var(580,2,'charge')**2 * (0.02*get_var(580,2,'charge'))
    CPQ750_set1 = get_var(750,1,'ptrig6_accp')/get_var(750,1,'charge') ; CPQ750_set1_err = get_var(750,1,'ptrig6_accp')/get_var(750,1,'charge')**2 * (0.02*get_var(750,1,'charge'))
    CPQ750_set2 = get_var(750,2,'ptrig6_accp')/get_var(750,2,'charge') ; CPQ750_set2_err = get_var(750,2,'ptrig6_accp')/get_var(750,2,'charge')**2 * (0.02*get_var(750,2,'charge'))
    CPQ750_set3 = get_var(750,3,'ptrig6_accp')/get_var(750,3,'charge') ; CPQ750_set3_err = get_var(750,3,'ptrig6_accp')/get_var(750,3,'charge')**2 * (0.02*get_var(750,3,'charge'))

    B.plot_exp(get_var(80,1,'Run'),   CPQ80_set1*0.06,   CPQ80_set1_err*0.06,  marker='^', color='k', label='80 (set1), scaled x0.06' )
    B.plot_exp(get_var(580,1,'Run'),  CPQ580_set1,  CPQ580_set1_err, marker='^', color='b', label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  CPQ580_set2,  CPQ580_set2_err, marker='^', color='g', label='580 (set2)' )
    B.plot_exp(get_var(750,1,'Run'),  CPQ750_set1,  CPQ750_set1_err, marker='^', color='r', label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  CPQ750_set2,  CPQ750_set2_err, marker='^', color='m', label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  CPQ750_set3,  CPQ750_set3_err, marker='^', color='c', label='750 (set3)' )

    B.pl.xlabel('Run Number')
    B.pl.ylabel('Accepted Coin. Triggers / mC')
    B.pl.title('Accepted Coincidence Triggers / Charge vs. Run Number')

    B.pl.legend(loc='upper right')
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/counts_per_charge.pdf')
    #-------------------------------------------------------
    
    #---------------Plot Run vs Total Live Time---------------
    #Get the average
    #w80 = 1. / (get_var(80,1,'tLT')*0.05)**2
    #w580_1 = 1. / (get_var(580,1,'tLT')*0.05)**2
    #w580_2 = 1. / (get_var(580,2,'tLT')*0.05)**2
    #w750_1 = 1. / (get_var(750,1,'tLT')*0.05)**2
    #w750_2 = 1. / (get_var(750,2,'tLT')*0.05)**2
    #w750_3 = 1. / (get_var(750,3,'tLT')*0.05)**2
    #x80 = np.average(get_var(80,1,'tLT')) + np.average(get_var(580,1,'tLT')) +  np.average(get_var(580,2,'tLT')) +  np.average(get_var(750,1,'tLT')) +  np.average(get_var(750,2,'tLT')) + np.average(get_var(750,3,'tLT'))) / 6.

    avgLT = ( np.average(get_var(80,1,'tLT')) + np.average(get_var(580,1,'tLT')) +  np.average(get_var(580,2,'tLT')) +  np.average(get_var(750,1,'tLT')) +  np.average(get_var(750,2,'tLT')) + np.average(get_var(750,3,'tLT'))) / 6.

    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'tLT'), get_var(80,1,'tLT')*0.05,  marker='s', color='k', label='80 (set1)', zorder=2 )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'tLT'), get_var(580,1,'tLT')*0.05,  marker='s', color='b', label='580 (set1)', zorder=2 )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'tLT'), get_var(580,2,'tLT')*0.05,  marker='s', color='g', label='580 (set2)', zorder=2 )  
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'tLT'), get_var(750,1,'tLT')*0.05,  marker='s', color='r', label='750 (set1)', zorder=2 )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'tLT'), get_var(750,2,'tLT')*0.05,  marker='s', color='m', label='750 (set2)', zorder=2 )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'tLT'), get_var(750,3,'tLT')*0.05,  marker='s', color='c', label='750 (set3)', zorder=2 )
    B.pl.xlim(3288, 3410)
    B.pl.ylim(0.80, 1.05)
    B.pl.axhline(y=avgLT, color='k', linestyle='--', label='avg=%f'%(avgLT), zorder=10)

    B.pl.xlabel('Run Number')
    B.pl.ylabel('Total Live Time')
    B.pl.title('Total EDMT Live Time vs. Run Number')

    B.pl.legend(loc='upper right')
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/total_livetime.pdf')

    #-------------------------------------------------------
    
    
    #---------------Plot Run vs Tracking Efficiencies---------------

    avg_htrk = ( np.average(get_var(80,1,'hTrkEff')) + np.average(get_var(580,1,'hTrkEff')) + np.average(get_var(580,2,'hTrkEff')) + np.average(get_var(750,1,'hTrkEff')) + np.average(get_var(750,2,'hTrkEff')) + np.average(get_var(750,3,'hTrkEff')) ) / 6.
    avg_etrk = ( np.average(get_var(80,1,'eTrkEff')) + np.average(get_var(580,1,'eTrkEff')) + np.average(get_var(580,2,'eTrkEff')) + np.average(get_var(750,1,'eTrkEff')) + np.average(get_var(750,2,'eTrkEff')) + np.average(get_var(750,3,'eTrkEff')) ) / 6.

    print('avg_htrk=',avg_htrk)
    print('avg_etrk=',avg_etrk)

    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'hTrkEff'), get_var(80,1,'hTrkEff_err'),  marker='s', color='k', label='HMS' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'hTrkEff'), get_var(580,1,'hTrkEff_err'),  marker='s', color='b', label='HMS ' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'hTrkEff'), get_var(580,2,'hTrkEff_err'),  marker='s', color='g', label='HMS  ' )
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'hTrkEff'), get_var(750,1,'hTrkEff_err'),  marker='s', color='r', label='HMS   ' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'hTrkEff'), get_var(750,2,'hTrkEff_err'),  marker='s', color='m', label='HMS     ' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'hTrkEff'), get_var(750,3,'hTrkEff_err'),  marker='s', color='c', label='HMS      ' )

    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'eTrkEff'), get_var(80,1,'eTrkEff_err'),  marker='s', mec='k', mfc='white', ecolor='k', label='SHMS ---- 80(set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'eTrkEff'), get_var(580,1,'eTrkEff_err'),  marker='s', mec='b', mfc='white', ecolor='b', label='SHMS ---- 580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'eTrkEff'), get_var(580,2,'eTrkEff_err'),  marker='s', mec='g', mfc='white', ecolor='g', label='SHMS ---- 580 (set2)' )
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'eTrkEff'), get_var(750,1,'eTrkEff_err'),  marker='s', mec='r', mfc='white', ecolor='r', label='SHMS ---- 750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'eTrkEff'), get_var(750,2,'eTrkEff_err'),  marker='s', mec='m', mfc='white', ecolor='m', label='SHMS ---- 750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'eTrkEff'), get_var(750,3,'eTrkEff_err'),  marker='s', mec='c', mfc='white', ecolor='c', label='SHMS ---- 750 (set3)' )
        
    #B.pl.axhline(y=avg_htrk, color='k', linestyle='--', label='avg_htrk=%f'%(avg_htrk), zorder=10)
    #B.pl.axhline(y=avg_etrk, color='k', linestyle='--', label='avg_etrk=%f'%(avg_etrk), zorder=10)

    B.pl.xlabel('Run Number')
    B.pl.ylabel('Tracking Efficiency')
    B.pl.title('Tracking Efficiency vs. Run Number')
    B.pl.ylim(0.9, 1.05)
    B.pl.xlim(3285, 3400)

    B.pl.legend(ncol=2, loc='upper right')
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/tracking_eff.pdf')

    
    #---------------Plot Run vs Target Boiling Factor---------------
    #tgt_Boil = 1 - m*I  ; apply errror propagation
    m =  0.00080029      #LD2 slope
    sig_m = 0.00007037   
    dI_I = 2.0/100.      #relative error in current (assume 2%)
    sig_I_80 = dI_I * get_var(80,1,'avg_current') 
    sig_I_580set1 = dI_I * get_var(580,1,'avg_current') ;  sig_I_580set2 = dI_I * get_var(580,2,'avg_current')
    sig_I_750set1 = dI_I * get_var(750,1,'avg_current') ;  sig_I_750set2 = dI_I * get_var(750,2,'avg_current') ;  sig_I_750set3 = dI_I * get_var(750,3,'avg_current')

    tb_80_err = np.sqrt(  get_var(80,1,'avg_current')**2 * sig_m**2 + m**2 * sig_I_80**2  ) 
    tb_580set1_err = np.sqrt(  get_var(580,1,'avg_current')**2 * sig_m**2 + m**2 * sig_I_580set1**2  ) 
    tb_580set2_err = np.sqrt(  get_var(580,2,'avg_current')**2 * sig_m**2 + m**2 * sig_I_580set2**2  ) 
    tb_750set1_err = np.sqrt(  get_var(750,1,'avg_current')**2 * sig_m**2 + m**2 * sig_I_750set1**2  ) 
    tb_750set2_err = np.sqrt(  get_var(750,2,'avg_current')**2 * sig_m**2 + m**2 * sig_I_750set2**2  ) 
    tb_750set3_err = np.sqrt(  get_var(750,3,'avg_current')**2 * sig_m**2 + m**2 * sig_I_750set3**2  ) 

    avg_boil = (np.average(get_var(80,1,'tgtBoil_factor')) + np.average(get_var(580,1,'tgtBoil_factor')) + np.average(get_var(580,2,'tgtBoil_factor')) + np.average(get_var(750,1,'tgtBoil_factor')) + np.average(get_var(750,2,'tgtBoil_factor')) + np.average(get_var(750,3,'tgtBoil_factor'))) / 6.

    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'tgtBoil_factor'), tb_80_err,  marker='s', ecolor = 'k', mec='k', mfc='white', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'tgtBoil_factor'), tb_580set1_err,  marker='s',  ecolor = 'b', mec='b', mfc='white',label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'tgtBoil_factor'), tb_580set2_err,  marker='s',  ecolor = 'g', mec='g', mfc='white',label='580 (set2)' )
  
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'tgtBoil_factor'), tb_750set1_err,  marker='s',  ecolor = 'r', mec='r', mfc='white',label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'tgtBoil_factor'), tb_750set2_err,  marker='s',  ecolor = 'm', mec='m', mfc='white',label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'tgtBoil_factor'), tb_750set3_err,  marker='s',  ecolor = 'c', mec='c', mfc='white',label='750 (set3)' )
    
    B.pl.axhline(y=avg_boil, color='k', linestyle='--', label='avg_boil=%f'%(avg_boil), zorder=10)

    B.pl.xlabel('Run Number')
    B.pl.ylabel('Target Boiling Factor')
    B.pl.title('LD2 Boiling Factor vs. Run Number')
    B.pl.ylim(0.9, 1.05)
    B.pl.xlim(3285, 3400)

    B.pl.legend(loc='upper right')
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/target_boil.pdf')

    #---------------------------------------------------------------
    
    
    #---------------Plot Run vs Average Beam Current---------------
    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'avg_current'), get_var(80,1,'avg_current')*0.02,  marker='s', color='k', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'avg_current'), get_var(580,1,'avg_current')*0.02,  marker='s', color='b', label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'avg_current'), get_var(580,2,'avg_current')*0.02,  marker='s', color='g', label='580 (set2)' )
  
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'avg_current'), get_var(750,1,'avg_current')*0.02,  marker='s', color='r', label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'avg_current'), get_var(750,2,'avg_current')*0.02,  marker='s', color='m', label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'avg_current'), get_var(750,3,'avg_current')*0.02,  marker='s', color='c', label='750 (set3)' )

    B.pl.xlabel('Run Number')
    B.pl.ylabel(r'Average Beam Current [\mu A]')
    B.pl.title('Average Beam Current vs. Run Number')
    B.pl.ylim(40, 70)
    B.pl.xlim(3285, 3400)

    B.pl.legend(loc='upper right')
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/beam_current.pdf')

    #-------------------------------------------------------
    

    
    #---------------Plot Run vs Trigger Rates---------------   
    B.pl.subplot(311)
    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'ptrig1_rate'),   marker='s', color='k', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'ptrig1_rate'),   marker='s', color='b', label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'ptrig1_rate'),  marker='s', color='g', label='580 (set2)' )
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'ptrig1_rate'),   marker='s', color='r', label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'ptrig1_rate'),  marker='s', color='m', label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'ptrig1_rate'),   marker='s', color='c', label='750 (set3)' )
    B.pl.title('Trigger Rates vs. Run Number')
    B.pl.ylabel(r'SHMS Trigger Rate [kHz]')
    B.pl.xlim(3285, 3410)
    B.pl.legend(loc='upper right')

    #---------------------------
    B.pl.subplot(312)
    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'ptrig4_rate')*0.2,   marker='s', color='k', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'ptrig4_rate'),   marker='s', color='b', label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'ptrig4_rate'),  marker='s', color='g', label='580 (set2)' )
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'ptrig4_rate'),   marker='s', color='r', label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'ptrig4_rate'),  marker='s', color='m', label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'ptrig4_rate'),   marker='s', color='c', label='750 (set3)' )
    B.pl.ylabel(r'HMS Trigger Rate [kHz]')
    B.pl.text(3400, 0.165, '$P_{m}=80$ MeV/c \n Scaled $x0.2$')
    B.pl.title('')
    B.pl.xlim(3285, 3410)
    #---------------------------
    B.pl.subplot(313)
    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'ptrig6_rate')*1000*0.02,   marker='s', color='k', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'ptrig6_rate')*1000,   marker='s', color='b', label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'ptrig6_rate')*1000,  marker='s', color='g', label='580 (set2)' )
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'ptrig6_rate')*1000,   marker='s', color='r', label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'ptrig6_rate')*1000,  marker='s', color='m', label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'ptrig6_rate')*1000,   marker='s', color='c', label='750 (set3)' )
    B.pl.ylabel(r'Coincidence Trigger Rate [Hz]')
    B.pl.text(3400, 3.0, '$P_{m}=80$ MeV/c \nScaled $x0.02$')
    B.pl.title('')

    B.pl.subplots_adjust(hspace=.001)
    B.pl.xlabel('Run Number')
    #B.pl.ylabel(r'Trigger Rate [kHz]')
    B.pl.xlim(3285, 3410)
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/trigger_rates.pdf')

    #-------------------------------------------------------

    
    
    #---------------Plot Run vs BPM position---------------
    #Assume relative uncertainty of 0.1 mm (0.01 cm) (from D. Gaskell) for now
    bpm_err = 0.01 #in cm,  df/f
    B.plot_exp(get_heep_var('Run'),  get_heep_var('xBPM'), get_heep_var('xBPM')*bpm_err, marker='D', color='gray', label='H(e,e\'p)' )

    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'xBPM'),  get_var(80,1,'xBPM')*bpm_err, marker='s', color='k', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'xBPM'),get_var(580,1,'xBPM')*bpm_err,  marker='s', color='b', label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'xBPM'),get_var(580,2,'xBPM')*bpm_err,  marker='s', color='g', label='580 (set2)' )
  
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'xBPM'), get_var(750,1,'xBPM')*bpm_err, marker='s', color='r', label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'xBPM'), get_var(750,2,'xBPM')*bpm_err, marker='s', color='m', label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'xBPM'), get_var(750,3,'xBPM')*bpm_err, marker='s', color='c', label='750 (set3)' )
    #-------------
    B.plot_exp(get_heep_var('Run'),  get_heep_var('yBPM'), get_heep_var('yBPM')*bpm_err, marker='D', color='gray', mfc='none', label='H(e,e\'p)' )

    B.plot_exp(get_var(80,1,'Run'),  get_var(80,1,'yBPM'),  get_var(80,1,'yBPM')*bpm_err, marker='s', color='k', mfc='none', label='80 (set1)' )
    B.plot_exp(get_var(580,1,'Run'),  get_var(580,1,'yBPM'),get_var(580,1,'yBPM')*bpm_err,  marker='s', color='b',  mfc='none',label='580 (set1)' )
    B.plot_exp(get_var(580,2,'Run'),  get_var(580,2,'yBPM'),get_var(580,2,'yBPM')*bpm_err,  marker='s', color='g', mfc='none', label='580 (set2)' )
  
    B.plot_exp(get_var(750,1,'Run'),  get_var(750,1,'yBPM'),get_var(750,1,'yBPM')*bpm_err,  marker='s', color='r',  mfc='none',label='750 (set1)' )
    B.plot_exp(get_var(750,2,'Run'),  get_var(750,2,'yBPM'),get_var(750,2,'yBPM')*bpm_err,  marker='s', color='m',  mfc='none',label='750 (set2)' )
    B.plot_exp(get_var(750,3,'Run'),  get_var(750,3,'yBPM'),get_var(750,3,'yBPM')*bpm_err,  marker='s', color='c',  mfc='none',label='750 (set3)' )


    B.pl.xlabel('Run Number')
    B.pl.ylabel(r'BPM Position [cm]')
    B.pl.title('Beam Position Monitor vs. Run Number')
    B.pl.text(3290, 0.0325, r'X BPMs')
    B.pl.text(3290, 0.02, r'Y BPMs')
    B.pl.xlim(3285, 3410)

    B.pl.legend()
    B.pl.show()
    B.pl.savefig(dir_name_misc+'/beam_position.pdf')
    
    #------------------------------------------------------------

def main():
    print('Entering Main . . .')

    #plot_data_sets()

    #plot_theory_sets('pwia')
    #plot_theory_sets('fsi')
    #plot_Xsec_vs_thnq()
    plot_final()
    #make_prl_plots()
    # plot_report()
if __name__=="__main__":
    main()


