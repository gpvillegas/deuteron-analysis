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
from matplotlib import colormaps as cmp

# Set default font to times new roman
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14
}
B.pl.rc('font', **font)
#%%
def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr

#%%

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
fname = f'yield_n_xsec/bin_centering/acos_phi/{pm_set}_{model}_bc_corr_{data_set}.txt'

f = dfile(fname)

#Get Relevant headers to plot
pm = f['yb']    #missing momentum bin
thnq = f['xb']  #theta_nq bin 
pm_avg = f['pm_avg']

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
    
    rel_error = fsiRC_dataXsec_err[thnq==ithnq]/fsiRC_dataXsec[thnq==ithnq]
    mask = rel_error < 0.5

    print('thnq=',ithnq)
         
    B.pl.clf()
    B.pl.figure(figsize=(8,6))
    
    #Plot data Xsec
    B.plot_exp(pm[thnq==ithnq][mask], fsiRC_dataXsec[thnq==ithnq][mask], fsiRC_dataXsec_err[thnq==ithnq][mask], marker='o', color='k', logy=True, label='Rad. Corrected')
    B.plot_exp(pm_avg[thnq==ithnq][mask], fsiRC_dataXsec_fsibc_corr[thnq==ithnq][mask], fsiRC_dataXsec_fsibc_corr_err[thnq==ithnq][mask], marker='o', color='r', logy=True, label='Rad. + VC. Corrected')
    
    B.pl.legend() 

    B.pl.title('$\\theta_{nq}=$'+f'{ithnq}'+'$ \pm 5$ deg', fontsize=18)                          #Set common title
    B.pl.xlabel(r'$P_{m}$ [GeV]',  fontsize=18)                           #Set common x-label
    B.pl.ylabel(r'Data Cross Section $\mu b$ / MeV sr$^{2}$', fontsize=18)  #Set common y-label

    B.pl.savefig(dir_name+f'/{pm_set}{data_set}_dataXsecBC_{ithnq}.pdf')
  
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
    
#%%

fname_120_1 = 'yield_n_xsec/bin_centering/acos_phi/pm_120_jml_bc_corr_set_1.txt'
fname_120_2 = 'yield_n_xsec/bin_centering/acos_phi/pm_120_jml_bc_corr_set_2.txt'
fname_580_1 = 'yield_n_xsec/bin_centering/acos_phi/pm_580_jml_bc_corr_set_1.txt'
fname_580_2 = 'yield_n_xsec/bin_centering/acos_phi/pm_580_jml_bc_corr_set_2.txt'
fname_800 = 'yield_n_xsec/bin_centering/acos_phi/pm_800_jml_bc_corr_set_1.txt'
fname_900 = 'yield_n_xsec/bin_centering/acos_phi/pm_900_jml_bc_corr_set_1.txt'

f_120_1 = dfile(fname_120_1)
f_120_2 = dfile(fname_120_2)
f_580_1 = dfile(fname_580_1)
f_580_2 = dfile(fname_580_2)
f_800 = dfile(fname_800)
f_900 = dfile(fname_900)

#Get Relevant headers to plot
pm = f_120_1['yb']    #missing momentum bin
thnq = f_120_1['xb']  #theta_nq bin 

#Get Bin-Centering Factors
# pm 120
bc_fact_fsi_120_1 = convert2NaN(f_120_1['bc_fact_fsi'],-1)
bc_fact_fsi_err_120_1 = convert2NaN(f_120_1['bc_fact_fsi_err'],-1)

bc_fact_pwia_120_1 = convert2NaN(f_120_1['bc_fact_pwia'],-1)
bc_fact_pwia_err_120_1 = convert2NaN(f_120_1['bc_fact_pwia_err'],-1)

bc_fact_fsi_120_2 = convert2NaN(f_120_2['bc_fact_fsi'],-1)
bc_fact_fsi_err_120_2 = convert2NaN(f_120_2['bc_fact_fsi_err'],-1)

bc_fact_pwia_120_2 = convert2NaN(f_120_2['bc_fact_pwia'],-1)
bc_fact_pwia_err_120_2 = convert2NaN(f_120_2['bc_fact_pwia_err'],-1)

# pm 580
bc_fact_fsi_580_1 = convert2NaN(f_580_1['bc_fact_fsi'],-1)
bc_fact_fsi_err_580_1 = convert2NaN(f_580_1['bc_fact_fsi_err'],-1)

bc_fact_pwia_580_1 = convert2NaN(f_580_1['bc_fact_pwia'],-1)
bc_fact_pwia_err_580_1 = convert2NaN(f_580_1['bc_fact_pwia_err'],-1)

bc_fact_fsi_580_2 = convert2NaN(f_580_2['bc_fact_fsi'],-1)
bc_fact_fsi_err_580_2 = convert2NaN(f_580_2['bc_fact_fsi_err'],-1)

bc_fact_pwia_580_2 = convert2NaN(f_580_2['bc_fact_pwia'],-1)
bc_fact_pwia_err_580_2 = convert2NaN(f_580_2['bc_fact_pwia_err'],-1)

# pm 800
bc_fact_fsi_800 = convert2NaN(f_800['bc_fact_fsi'],-1)
bc_fact_fsi_err_800 = convert2NaN(f_800['bc_fact_fsi_err'],-1)

bc_fact_pwia_800 = convert2NaN(f_800['bc_fact_pwia'],-1)
bc_fact_pwia_err_800 = convert2NaN(f_800['bc_fact_pwia_err'],-1)

# pm 900
bc_fact_fsi_900 = convert2NaN(f_900['bc_fact_fsi'],-1)
bc_fact_fsi_err_900 = convert2NaN(f_900['bc_fact_fsi_err'],-1)

bc_fact_pwia_900 = convert2NaN(f_900['bc_fact_pwia'],-1)
bc_fact_pwia_err_900 = convert2NaN(f_900['bc_fact_pwia_err'],-1)

thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

for i, ithnq in enumerate(thnq_arr): #enumerate(thnq_arr)
    
    rel_err_120 = bc_fact_pwia_err_120_1[thnq==ithnq]/bc_fact_pwia_120_1[thnq==ithnq]
    mask_120 = rel_err_120 < 0.5
    
    rel_err_580 = bc_fact_pwia_err_580_1[thnq==ithnq]/bc_fact_pwia_580_1[thnq==ithnq]
    mask_580 = rel_err_580 < 0.5
    
    rel_err_800 = bc_fact_pwia_err_800[thnq==ithnq]/bc_fact_pwia_800[thnq==ithnq]
    mask_800 = rel_err_800 < 0.5
    
    rel_err_900 = bc_fact_pwia_err_900[thnq==ithnq]/bc_fact_pwia_900[thnq==ithnq]
    mask_900 = rel_err_900 < 0.5

    print('thnq=',ithnq)
    #--------------PLOT BC. Corr. Factor------------------
    B.pl.figure(figsize=(8,6))

    #Plot BC factor
    # pm 120
    B.plot_exp(pm[thnq==ithnq][mask_120], bc_fact_pwia_120_1[thnq==ithnq][mask_120], bc_fact_pwia_err_120_1[thnq==ithnq][mask_120], marker='s', color='gray', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq][mask_120], bc_fact_fsi_120_1[thnq==ithnq][mask_120], bc_fact_fsi_err_120_1[thnq==ithnq][mask_120], marker='s', markerfacecolor='white', color='black', logy=False, label='FSI')
    
    B.plot_exp(pm[thnq==ithnq][mask_120], bc_fact_pwia_120_1[thnq==ithnq][mask_120], bc_fact_pwia_err_120_1[thnq==ithnq][mask_120], marker='s', color=cmp['Paired'](5), logy=False, markersize=6)
    B.plot_exp(pm[thnq==ithnq][mask_120], bc_fact_fsi_120_1[thnq==ithnq][mask_120], bc_fact_fsi_err_120_1[thnq==ithnq][mask_120], marker='s', markerfacecolor='white', color=cmp['tab20b'](13), logy=False, label='$P_m =$ 120 MeV', markersize=6)
    
    # B.plot_exp(pm[thnq==ithnq], bc_fact_pwia_120_2[thnq==ithnq], bc_fact_pwia_err_120_2[thnq==ithnq], marker='s', color='c', logy=False, label='$P_m =$ 120 set 2')
    # B.plot_exp(pm[thnq==ithnq], bc_fact_fsi_120_2[thnq==ithnq], bc_fact_fsi_err_120_2[thnq==ithnq], marker='s', markerfacecolor='white', color='c', logy=False)
    
    # pm 580
    B.plot_exp(pm[thnq==ithnq][mask_580], bc_fact_pwia_580_1[thnq==ithnq][mask_580], bc_fact_pwia_err_580_1[thnq==ithnq][mask_580], marker='s', color=cmp['tab20c'](5), logy=False, markersize=6)
    B.plot_exp(pm[thnq==ithnq][mask_580], bc_fact_fsi_580_1[thnq==ithnq][mask_580], bc_fact_fsi_err_580_1[thnq==ithnq][mask_580], marker='s', markerfacecolor='white', color=cmp['tab20c'](4), logy=False, label='$P_m =$ 580 MeV', markersize=6)

    # pm 800
    B.plot_exp(pm[thnq==ithnq][mask_800], bc_fact_pwia_800[thnq==ithnq][mask_800], bc_fact_pwia_err_800[thnq==ithnq][mask_800], marker='s', color=cmp['tab20c'](9), logy=False, markersize=6)
    B.plot_exp(pm[thnq==ithnq][mask_800], bc_fact_fsi_800[thnq==ithnq][mask_800], bc_fact_fsi_err_800[thnq==ithnq][mask_800], marker='s', markerfacecolor='white', color=cmp['tab20c'](8), logy=False, label='$P_m =$ 800 MeV', markersize=6)

    # pm 900
    B.plot_exp(pm[thnq==ithnq][mask_900], bc_fact_pwia_900[thnq==ithnq][mask_900], bc_fact_pwia_err_900[thnq==ithnq][mask_900], marker='s', color=cmp['tab20b'](17), logy=False, markersize=6)
    B.plot_exp(pm[thnq==ithnq][mask_900], bc_fact_fsi_900[thnq==ithnq][mask_900], bc_fact_fsi_err_900[thnq==ithnq][mask_900], marker='s', markerfacecolor='white', color=cmp['tab20b'](16), logy=False, label='$P_m =$ 900 MeV', markersize=6)

    xmi,xma = B.pl.xlim()
    B.pl.hlines([0.9,1.1,0.8,1.2],xmin=xmi,xmax=xma,linestyles='--',colors=['r','r','black','black'])
    B.pl.legend(fontsize='small') 

    B.pl.title(' $\\theta_{nq}=$'+f'{ithnq}'+'$ \pm 5$ deg', fontsize=18)                          #Set common title
    B.pl.xlabel('$P_{m}$ [GeV]',  fontsize=18)                           #Set common x-label
    B.pl.ylabel('VC Factor', fontsize=18)  #Set common y-label

    B.pl.savefig(f'yield_n_xsec/bin_centering/plots/BCfactor_allpm_{ithnq}.pdf')


