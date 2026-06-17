#This code checks the Bin-Centering corrections are done properly and make
#sense. 
#The code makes plots of: 
#1) BC corr. ratio for PWIA, FSI
#2) dataXsec before/after FSI BC. Corrected
#3) fsiXsec_avg_simc, fsiXsec_theory  (the two cross sections ratio give bc corr. factor)

from LT.datafile import dfile
import LT.box as B
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
import sys
import os

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#User Input (Usage: python check_bc_corr.py 80 1 )
# pm_set = int(sys.argv[1])
# data_set = int(sys.argv[2])
# sys_ext = sys.argv[3]  
pm_set, model, data_set = input('Enter pm setting, model, and data set: \n').split()

#Make relevant directories to store plots
dir_name = f"./yield_n_xsec/bin_centering/plots/{pm_set}"

#check if directory exists, else creates it.
if not os.path.exists(dir_name):
    os.makedirs(dir_name)

#Read data file
# if pm_set==80:
#     fname = sys_ext + '/pm%i_laget_bc_corr.txt' % (pm_set)
# else:
fname = f'yield_n_xsec/bin_centering/{pm_set}_{model}_bc_corr_{data_set}.txt'

f = dfile(fname)

#Get Relevant headers to plot
pm = f['yb']    #missing momentum bin
thnq = f['xb']  #theta_nq bin 

#Get Bin-Centering Factors
bc_fact_fsi = f['bc_fact_fsi']
bc_fact_fsi_err = f['bc_fact_fsi_err']

bc_fact_pwia = f['bc_fact_pwia']
bc_fact_pwia_err = f['bc_fact_pwia_err']

#Get Data Xsec
fsiRC_dataXsec = f['fsiRC_dataXsec']               #only rad. corrected  data
fsiRC_dataXsec_err = f['fsiRC_dataXsec_err']
fsiRC_dataXsec_fsibc_corr = f['fsiRC_dataXsec_fsibc_corr']          #rad + b.c. corr.  data
fsiRC_dataXsec_fsibc_corr_err = f['fsiRC_dataXsec_fsibc_corr_err']

#Get SIMC avg. Xsec
# fsiXsec_SIMC_avg = f['fsiXsec_SIMC_avg']
# fsiXsec_SIMC_avg_err = f['fsiXsec_SIMC_avg_err']

#Get Laget Theory Xsec (@ the avg. kinematics)
fsiXsec_theory = f['fsiXsec_theory']

#Convert to nan is value is -1:
convert2NaN(bc_fact_fsi, value=-1)
convert2NaN(bc_fact_fsi_err, value=-1)
convert2NaN(bc_fact_pwia, value=-1)
convert2NaN(bc_fact_pwia_err, value=-1)
convert2NaN(fsiRC_dataXsec, value=-1)
convert2NaN(fsiRC_dataXsec_err, value=-1)
convert2NaN(fsiRC_dataXsec_fsibc_corr, value=-1)
convert2NaN(fsiRC_dataXsec_fsibc_corr_err, value=-1)
# convert2NaN(fsiXsec_SIMC_avg, value=-1)
# convert2NaN(fsiXsec_SIMC_avg_err, value=-1)
convert2NaN(fsiXsec_theory, value=-1)


thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

for i, ithnq in enumerate(thnq_arr):
    
    print('thnq=',ithnq)
         
    B.pl.clf()
    B.pl.figure(i)
    
    #Plot data Xsec
    B.plot_exp(pm[thnq==ithnq], fsiRC_dataXsec[thnq==ithnq], fsiRC_dataXsec_err[thnq==ithnq], marker='o', color='k', logy=True, label='Rad. Corrected')
    B.plot_exp(pm[thnq==ithnq], fsiRC_dataXsec_fsibc_corr[thnq==ithnq], fsiRC_dataXsec_fsibc_corr_err[thnq==ithnq], marker='o', color='r', logy=True, label='Rad. + BC. Corrected')
    
    B.pl.legend() 

    B.pl.title('Data Cross Sections, $\\theta_{nq}=$'+f'{ithnq}'+'$ \pm 5$ deg', fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{r}$ (GeV/c)',  fontsize=12)                           #Set common x-label
    B.pl.ylabel(r'Data Cross Section $\mu b$ / MeV sr$^{2}$', fontsize=12)  #Set common y-label

    B.pl.savefig(dir_name+f'/{pm_set}{data_set}_dataXsecBC_{ithnq}.png')
  
    #--------------PLOT BC. Corr. Factor------------------
    B.pl.clf()
    B.pl.figure(i+1)

    #Plot BC factor
    B.plot_exp(pm[thnq==ithnq], bc_fact_fsi[thnq==ithnq], bc_fact_fsi_err[thnq==ithnq], marker='o', color='k', logy=False, label='Bin-Centering Factor, FSI')
    B.plot_exp(pm[thnq==ithnq], bc_fact_pwia[thnq==ithnq], bc_fact_pwia_err[thnq==ithnq], marker='o', color='r', logy=False, label='Bin-Centering Factor, PWIA')
    
    xmi,xma = B.pl.xlim()
    B.pl.hlines([0.9,1.1],xmin=xmi,xmax=xma,linestyles='--',colors='r')
    B.pl.legend() 

    B.pl.title('Bin Centering Corr. Factor, $\\theta_{nq}=$'+f'{ithnq}'+'$ \pm 5$ deg', fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{r}$ (GeV/c)',  fontsize=12)                           #Set common x-label
    B.pl.ylabel('BC Factor', fontsize=12)  #Set common y-label

    B.pl.savefig(dir_name+f'/{pm_set}{data_set}_BCfactor_{ithnq}.png')
  
    #------Plot Average SIMC Xsed and theory Xsec-----
    B.pl.clf()
    B.pl.figure(i+2)
    
    #Plot data Xsec
    # B.plot_exp(pm[thnq==ithnq], fsiXsec_SIMC_avg[thnq==ithnq], fsiXsec_SIMC_avg_err[thnq==ithnq], marker='o', color='k', logy=True, label='Avg. SIMC Xsec')
    B.plot_exp(pm[thnq==ithnq], fsiXsec_theory[thnq==ithnq], marker='o', color='r', logy=True, label='Laget Theory Xsec')
    
    B.pl.legend() 

    B.pl.title('Model Cross Sections, $\\theta_{nq}=$'+f'{ithnq}'+'$ \pm 5$ deg', fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{r}$ (GeV/c)',  fontsize=12)                           #Set common x-label
    B.pl.ylabel(r'Data Cross Section $\mu b$ / MeV sr$^{2}$', fontsize=12)  #Set common y-label

    B.pl.savefig(dir_name+f'/{pm_set}{data_set}_avgXsec_{ithnq}.png')
    
