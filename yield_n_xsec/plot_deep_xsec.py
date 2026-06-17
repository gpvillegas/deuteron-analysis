#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 12:53:24 2025

@author: gvill

This program reads in a bin file with xsec histograms and plots them

"""
import LT.box as B

import pickle as P
from matplotlib import colormaps as cmp
# import matplotlib.lines as mlines
from matplotlib.ticker import MultipleLocator
import numpy as np

# Set default font to times new roman
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14
}
B.pl.rc('font', **font)
#%% functions
def set_nans(histo,set_nan_to=0.):
    array1 = histo.bin_content
    array2 = histo.bin_error
    for i in range(len(array1)):
        if array1[i] == array1[i]:
            array1[i] = array1[i]
        else:
            array1[i] = set_nan_to
        
        if array2[i] == array2[i]:
            array2[i] = array2[i]
        else:
            array2[i] = set_nan_to
    return

def set_inf(inpt,set_inf_to=0.):
    is_inf = np.where(inpt == np.inf)
    for i in is_inf[0]:
        inpt[i] = set_inf_to
    print(f'Set {is_inf[0].size} inf to {set_inf_to}\n')
    return

#%% load pickle file with histograms
plot_setting = input('Select setting to plot:\npm_120, pm_580, pm_800, pm_900\n')

filename = f'yield_n_xsec/{plot_setting}_histos_comm_bin.pcl'
with open(filename, 'rb') as f:
    histos = P.load(f)
    
filename2 = f'yield_n_xsec/{plot_setting}_histos_comm_bin2.pcl'
with open(filename2, 'rb') as f:
    histos.update(P.load(f))

filename3 = f'yield_n_xsec/{plot_setting}_histos_comm_bin3.pcl'
with open(filename3, 'rb') as f:
    histos.update(P.load(f))
#%% plot 2d histograms
save = input('save all figures?\n (0) no\n (1) yes\n')

nrow = 2
ncol = 5

fig, ax = B.pl.subplots(nrow,ncol,figsize=(15,8),
                        constrained_layout=True)
B.pl.suptitle(f'{plot_setting} 2D histograms',fontsize=14)
i = 1

for s in histos['data_raw']:
    # data_raw
    ax = B.pl.subplot(nrow,ncol,i)
    histos['data_raw'][s].plot(colormap=cmp['jet'])
    B.pl.title(f'Raw Counts - Data {s}')
    B.pl.xlabel('Neutron Recoil Angle, $\\theta_{nq}$ [$^\circ$]',fontsize=14)
    B.pl.ylabel('Missing Momentum, $P_m$ [GeV]',fontsize=14)
    i+=1

    # yield (Data with Cuts)
    ax = B.pl.subplot(nrow,ncol,i)
    histos[f'{plot_setting}_yield'][s].plot(colormap=cmp['jet'])
    B.pl.title(f'Raw Yield {s}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    i+=1
    
    # normalized yield
    ax = B.pl.subplot(nrow,ncol,i)
    histos[f'{plot_setting}_yield_norm'][s].plot(colormap=cmp['jet'])
    B.pl.title(f'Normalized Yield {s}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    i+=1
    
    # FSI xsec
    ax = B.pl.subplot(nrow,ncol,i)
    histos[f'{plot_setting}_FSI_xsec'][s].plot(colormap=cmp['jet'])
    B.pl.title(f'Data FSI Xsec {s}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    i+=1
    
    # PWIA xsec
    ax = B.pl.subplot(nrow,ncol,i)
    histos[f'{plot_setting}_PWIA_xsec'][s].plot(colormap=cmp['jet'])
    B.pl.title(f'Data PWIA Xsec {s}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    i+=1

nrow = 1
ncol = 5
i = 1
fig, ax = B.pl.subplots(nrow,ncol,figsize=(15,4),
                        constrained_layout=True)
B.pl.suptitle(f'{plot_setting} 2D histograms',fontsize=14)

ax = B.pl.subplot(nrow,ncol,i)
histos[f'{plot_setting}_FSI_radcorr_factor'].plot(colormap=cmp['jet'])
B.pl.title('SIMC FSI Rad. Corr. ($Y_{norad}/Y_{rad}$)')
B.pl.xlabel('Neutron Recoil Angle, $\\theta_{nq}$ [$^\circ$]',fontsize=14)
B.pl.ylabel('Missing Momentum, $P_m$ [GeV]',fontsize=14)
i+=1

ax = B.pl.subplot(nrow,ncol,i)
histos[f'{plot_setting}_yield_FSI_PS'].plot(colormap=cmp['jet'])
B.pl.title('SIMC FSI Phase Space')
B.pl.xlabel('')
B.pl.ylabel('')
i+=1

ax = B.pl.subplot(nrow,ncol,i)
histos[f'{plot_setting}_yield_PWIA_PS'].plot(colormap=cmp['jet'])
B.pl.title('SIMC FSI Phase Space')
B.pl.xlabel('')
B.pl.ylabel('')
i+=1

ax = B.pl.subplot(nrow,ncol,i)
histos[f'{plot_setting}_JML_fsi_xsec'].plot(colormap=cmp['jet'])
B.pl.title('SIMC JML FSI Cross Section')
B.pl.xlabel('')
B.pl.ylabel('')
i+=1

ax = B.pl.subplot(nrow,ncol,i)
histos[f'{plot_setting}_JML_pwia_xsec'].plot(colormap=cmp['jet'])
B.pl.title('SIMC JML PWIA Cross Section')
B.pl.xlabel('')
B.pl.ylabel('')
i+=1

if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_2D_histograms.png')
#%%
# plot raw yield projections
rows = 4
cols = 5
i=1
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Raw Yield Projections')
for h in histos[f'{plot_setting}_yield_proj']:
    ax=B.pl.subplot(rows,cols,i)
    # print(h.bin_error)
    # h.bin_error = raw_yield_err[i-1]
    # print(h.bin_error)
    h['set_1'].plot(facecolor=cmp['Paired'].colors[0])
    h['set_1'].plot_exp(markersize=1.2,elinewidth=1.2,capsize=0)
    
    h['set_1'].plot(facecolor=cmp['Paired'].colors[0])
    h['set_1'].plot_exp(markersize=1.2,elinewidth=1.2,capsize=0)
    
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    i+=1

if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Raw_Yield_Proj.png')
#%%  
# plot simc radiative correction factor projections
FSI_radcorr = histos[f'{plot_setting}_FSI_radcorr_factor_proj']
PWIA_radcorr = histos[f'{plot_setting}_PWIA_radcorr_factor_proj']

rows = 4
cols = 5
i=1
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Radiative Correction Factor Projections')
for h in range(len(FSI_radcorr)):
    ax=B.pl.subplot(rows,cols,i)
    # print(h.bin_error)
    # h.bin_error = raw_yield_err[i-1]
    # print(h.bin_error)
    # h.plot(facecolor=cmp['Paired'].colors[0])
    FSI_radcorr[h].plot_exp(markersize=4,elinewidth=1.2,capsize=0,
                            color='black',label='FSI')
    PWIA_radcorr[h].plot_exp(markersize=4,elinewidth=1.2,capsize=0,
                             color='green',label='PWIA')

    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    if i == 1:
        B.pl.legend()
    i+=1

# if save:
#     B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Rad_Corr_Fac_Proj.png')

# plot the PWIA/FSI ratio for radcorr 
rows = 4
cols = 5
i=1
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Radiative Correction Factor PWIA/FSI Ratio Projections')
for h in range(len(FSI_radcorr)):
    ax=B.pl.subplot(rows,cols,i)
    r = PWIA_radcorr[h].bin_content/FSI_radcorr[h].bin_content
    re = PWIA_radcorr[h].bin_error/FSI_radcorr[h].bin_content
    ratio = B.histo(bin_center=PWIA_radcorr[h].bin_center,
                    bin_content=r,
                    bin_error=re,
                    title=PWIA_radcorr[h].title,
                    xlabel='',
                    ylabel='')
    ratio.plot_exp(markersize=4,elinewidth=1.2,capsize=0,
                            color='black',ignore_zeros=True)
    
    B.pl.hlines(1,xmin=0,xmax=1,linestyles='--',color='r')

    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    i+=1    
#%% plot simc phase space projections
rows = 4
cols = 5
i=1
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('SIMC Phase Space Projections')
for h in histos[f'{plot_setting}_JML_FSI_PS_proj']:
    ax=B.pl.subplot(rows,cols,i)
    # print(h.bin_error)
    # h.bin_error = raw_yield_err[i-1]
    # print(h.bin_error)
    h.plot(facecolor=cmp['Paired'].colors[0])
    h.plot_exp(markersize=1.2,elinewidth=1.2,capsize=0)
    
    # B.pl.yscale('log')
    B.pl.xlim(-0.1,1.5)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    i+=1
    
rows = 4
cols = 5
i=1
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('SIMC Phase Space Projections')
for h in histos[f'{plot_setting}_JML_PWIA_PS_proj']:
    ax=B.pl.subplot(rows,cols,i)
    # print(h.bin_error)
    # h.bin_error = raw_yield_err[i-1]
    # print(h.bin_error)
    h.plot(facecolor=cmp['Paired'].colors[0])
    h.plot_exp(markersize=1.2,elinewidth=1.2,capsize=0)
    
    # B.pl.yscale('log')
    B.pl.xlim(-0.1,1.5)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    i+=1

if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Phase_Space_Proj.png')    
    
#%% plot norm yield projections
norm_yield_proj = histos[f'{plot_setting}_yield_norm_proj']
simc_norm_yield_proj = histos[f'{plot_setting}_JML_Yield_FSI_rad_proj']                       

rows = 4
cols = 5
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Normalized Yield Projections')
first_plot = True
plot_PS = True
for i in range(len(norm_yield_proj['set_1'])):
    ax=B.pl.subplot(rows,cols,i+1)
    # print(h.bin_error)
    # h.bin_error = raw_yield_err[i-1]
    # print(h.bin_error)
    h1 = norm_yield_proj['set_1'][i]
    h2 = norm_yield_proj['set_2'][i]
    hs = simc_norm_yield_proj[i]
    
    hs.plot(facecolor=cmp['Paired'].colors[4],alpha=0.3)
    hs.plot_exp(markersize=2.5,elinewidth=1.2,capsize=0,
                color=cmp['Paired'].colors[5],label='SIMC')
    
    h1.plot(facecolor=cmp['Paired'].colors[0],alpha=0.3)
    h1.plot_exp(markersize=2.5,elinewidth=1.2,capsize=0,label='DATA set_1')
    
    h2.plot(facecolor=cmp['Paired'].colors[3],alpha=0.3)
    h2.plot_exp(markersize=2.5,elinewidth=1.2,capsize=0,label='DATA set_2')
    
    if first_plot:
        B.pl.legend()
        B.pl.ylabel('Counts')
        
        first_plot = False
 
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    B.pl.xlabel('$P_m$')  

if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Norm_Yield_Proj.png')
    
#%% plot cross sections projections
xsec_proj = histos[f'{plot_setting}_FSI_xsec_proj']
JML_fsi_xsec_proj = histos[f'{plot_setting}_JML_fsi_xsec_proj']
JML_pwia_xsec_proj = histos[f'{plot_setting}_JML_pwia_xsec_proj']
# CDB_fsi_xsec_proj = histos[f'{plot_setting}_CDB_fsi_xsec_proj']
# CDB_pwia_xsec_proj = histos[f'{plot_setting}_CDB_pwia_xsec_proj']

if plot_setting == 'pm120':
    pm80_xsec_proj = histos['pm80_xsec_proj']
elif plot_setting == 'pm580':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
elif plot_setting == 'pm800' or plot_setting == 'pm900':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
    pm750_xsec_proj = histos['pm750_xsec_proj']

rows = 4
cols = 5
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Cross Section Projections')
first_plot = True
for i in range(len(xsec_proj['set_1'])):
    ax=B.pl.subplot(rows,cols,i+1)
    h0 = xsec_proj['set_1'][i]
    h00 = xsec_proj['set_2'][i]
    h1 = JML_fsi_xsec_proj[i]
    h2 = JML_pwia_xsec_proj[i]
    # h3 = CDB_fsi_xsec_proj[i]
    # h4 = CDB_pwia_xsec_proj[i]

    
    h1.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                label='JML FSI')
    # B.plot_line(h1.bin_center, h1.bin_content,label='JML FSI')
    h2.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                label='JML PWIA')
    # B.plot_line(h2.bin_center, h2.bin_content,label='JML PWIA')
    # h4.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
    #             label='MS CD-Bonn PWIA')
    # h3.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
    #             label='MS CD-Bonn FSI')    
    if plot_setting == 'pm120':
        h3 = pm80_xsec_proj[i]
        h3.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                marker='^',color=cmp['Set1'].colors[0],label='pm_80 Comm Data')
    elif plot_setting == 'pm580':
        h3 = pm580_comm_xsec_proj[i]
        h3.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                marker='^',color=cmp['Set1'].colors[0],label='pm_580 Comm Data')
    elif plot_setting == 'pm800' or plot_setting == 'pm900':
        h3 = pm580_comm_xsec_proj[i]
        h4 = pm750_xsec_proj[i]
        h3.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                marker='^',color=cmp['Set1'].colors[0],label='pm_580 Comm Data')
        h4.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                marker='^',color=cmp['Set1'].colors[2],label='pm_750 Comm Data')
    
    h0.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                color='black',label='DATA set_1')
    h00.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                color='black',label='DATA set_2')
    
    B.pl.yscale('log')
    B.pl.xlim(-0.1,1.5)
    if first_plot:
        B.pl.legend()
        B.pl.ylabel('$\\frac{d^5\sigma}{d\Omega_e\Omega_pdE\'}$')
        first_plot = False
    
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    B.pl.xlabel('$P_m$')

B.pl.show()
    
if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Xsec_Proj.png')

# plot yields with and without radiative corrections
norm_yield_radcorr_proj = histos[f'{plot_setting}_yield_radcorr_proj']

rows = 4
cols = 5
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Rad. Corr. Yield Projections')
first_plot = True
for i in range(len(norm_yield_proj)):
    ax=B.pl.subplot(rows,cols,i+1)
    h0 = norm_yield_proj[i]
    h1 = norm_yield_radcorr_proj[i]
    
    h1.plot(facecolor=cmp['Paired'].colors[6])
    h1.plot_exp(markersize=1.2,elinewidth=1.2,capsize=0,
                color=cmp['Paired'].colors[7],label='Rad. Corr. Yield')
    
    h0.plot(facecolor=cmp['Paired'].colors[0])
    h0.plot_exp(markersize=1.2,elinewidth=1.2,capsize=0,
                color=cmp['Paired'].colors[1],label='No Rad. Corr. Yield')
    
    if first_plot:
        B.pl.legend()
        first_plot = False
    
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Rad_Corr_Yield_Proj.png')


#%% plot relative errors
rel_err_proj = histos[f'{plot_setting}_Relative_Stats_Error']
raw_yield_h = histos[f'{plot_setting}_yield_proj']

rows = 4
cols = 5
j=0
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Relative Error Projections [%]')
first_plot = True
for i in rel_err_proj:
    ax=B.pl.subplot(rows,cols,j+1)
    
    dy = rel_err_proj[i]
    x = raw_yield_h[j].bin_center
    y = np.zeros(x.size)
    
    B.plot_exp(x,y,dy,markersize=1.2,elinewidth=1.2,
                color=cmp['Paired'].colors[1])
    B.pl.hlines([-10,10],xmin=-0.15,xmax=1.5,linestyles='--',color='green')
    B.pl.hlines([-20,20],xmin=-0.15,xmax=1.5,linestyles='--',color='red')
    
    if first_plot:
        first_plot=False
        B.pl.text(1.25,12,'10%',color='green')
        B.pl.text(1.25,23,'20%',color='red')
    
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    j+=1

if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Relative_Stats_Error_Proj.png')
#%%
# rel_err_proj = histos[f'{plot_setting}_Relative_Stats_Error']['set_1']
# raw_yield_h = histos[f'{plot_setting}_yield_proj']['set_1']
rel_err_proj = histos[f'{plot_setting}_Relative_Stats_Error']
raw_yield_h = histos[f'{plot_setting}_yield_proj']

j=0
first_plot = True
for i in rel_err_proj:
    if j == 4 or j==7 or j==3:
        B.pl.figure()
        
        dy = rel_err_proj[i]
        x = raw_yield_h[j].bin_center
        y = np.zeros(x.size)
        
        B.plot_exp(x,y,dy,markersize=2,elinewidth=2.4,
                    color=cmp['Paired'].colors[1])
        B.pl.hlines([-10,10],xmin=-0.15,xmax=1.5,linestyles='--',color='green')
        B.pl.hlines([-20,20],xmin=-0.15,xmax=1.5,linestyles='--',color='red')
        
        B.pl.text(1.25,12,'10%',color='green')
        B.pl.text(1.25,23,'20%',color='red')
        
        B.pl.xlim(0.,1.2)
        B.pl.ylim(-25,25)
        B.pl.xlabel('$P_m$')
        B.pl.ylabel('Relative Error [%]')
        B.pl.title(i)
        ax = B.pl.gca()
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.025))
        
    j+=1

#%% plot cross section projections for specific angle bins
FSI_xsec_proj = histos[f'{plot_setting}_FSI_xsec_proj']
# PWIA_xsec_proj = histos[f'{plot_setting}_PWIA_xsec_proj']
simc_jmlfsi_xsec_proj = histos[f'{plot_setting}_JML_fsi_xsec_proj']
simc_jmlpwia_xsec_proj = histos[f'{plot_setting}_JML_pwia_xsec_proj']
simc_parfsi_xsec_proj = histos[f'{plot_setting}_PAR_fsi_xsec_proj']
simc_parpwia_xsec_proj = histos[f'{plot_setting}_PAR_pwia_xsec_proj']
simc_v18fsi_xsec_proj = histos[f'{plot_setting}_V18_fsi_xsec_proj']
simc_v18pwia_xsec_proj = histos[f'{plot_setting}_V18_pwia_xsec_proj']
simc_cdbfsi_xsec_proj = histos[f'{plot_setting}_CDB_fsi_xsec_proj']
simc_cdbpwia_xsec_proj = histos[f'{plot_setting}_CDB_pwia_xsec_proj']

if plot_setting == 'pm_120':
    pm80_xsec_proj = histos['pm80_xsec_proj']
elif plot_setting == 'pm_580':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
    pm750_xsec_proj = histos['pm750_xsec_proj']


for s in FSI_xsec_proj:
    first_plot = True
    for i in range(len(FSI_xsec_proj[s])):
        h0 = FSI_xsec_proj[s][i]
        # h1 = PWIA_xsec_proj[s][i]
        h2 = simc_jmlfsi_xsec_proj[i]
        h3 = simc_jmlpwia_xsec_proj[i]
        h4 = simc_parfsi_xsec_proj[i]
        h5 = simc_parpwia_xsec_proj[i]
        h6 = simc_v18fsi_xsec_proj[i]
        h7 = simc_v18pwia_xsec_proj[i]
        h8 = simc_cdbfsi_xsec_proj[i]
        h9 = simc_cdbpwia_xsec_proj[i]
        
        if h0.title[16:18] == '35' or h0.title[16:18] == '45' or\
            h0.title[16:18] == '75':
            B.pl.figure(figsize=(8,6))     
            h2.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='s',color=cmp['tab20'](2),label='JML FSI')
            h3.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='D',color=cmp['tab20'](3),label='JML PWIA')
            h4.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='s',color=cmp['tab20'](4),label='MS Paris FSI')
            h5.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='D',color=cmp['tab20'](5),label='MS Paris PWIA')
            h6.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='s',color=cmp['tab20'](6),label='MS V18 FSI')
            h7.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='D',color=cmp['tab20'](7),label='MS V18 PWIA')
            h8.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='s',color=cmp['tab20'](8),label='MS CD-Bonn FSI')
            h9.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='D',color=cmp['tab20'](9),label='MS CD-Bonn PWIA')
            if plot_setting == 'pm_120':
                h3 = pm80_xsec_proj[i]
                h3.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='^',color=cmp['Set1'].colors[0],label='pm_80 Comm Data')
            elif plot_setting == 'pm_580':
                h3 = pm580_comm_xsec_proj[i]
                h3.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='^',color=cmp['Set1'].colors[0],label='pm_580 Comm Data')
            elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
                h3 = pm580_comm_xsec_proj[i]
                h4 = pm750_xsec_proj[i]
                h3.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='^',color=cmp['Set1'].colors[0],label='pm_580 Comm Data')
                h4.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        marker='^',color=cmp['Set1'].colors[2],label='pm_750 Comm Data')
            
            h0.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                        color='black',label=f'FSI DATA {s}')
            # h1.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
            #             color='black',label=f'PWIA DATA {s}')
            
            B.pl.yscale('log')
            B.pl.xlim(-0.1,2.0)
            B.pl.ylabel('$\\frac{d^5\sigma}{d\Omega_e\Omega_pdE\'}$ [$\mu b \cdot sr^2 GeV^{-1}$]')
            B.pl.xlabel('$P_m$ [Gev/c]')
            B.pl.legend()
            ax = B.pl.gca()
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))

#%%
first_plot = True
for i in range(len(FSI_xsec_proj['set_1'])):
    h01 = FSI_xsec_proj['set_1'][i]
    set_nans(h01)
    h02 = FSI_xsec_proj['set_2'][i]
    set_nans(h02)
    h0 = h01+h02/2
    h0.title = h01.title
    h2 = simc_jmlfsi_xsec_proj[i]
    h3 = simc_jmlpwia_xsec_proj[i]
    h4 = simc_parfsi_xsec_proj[i]
    h5 = simc_parpwia_xsec_proj[i]
    h6 = simc_v18fsi_xsec_proj[i]
    h7 = simc_v18pwia_xsec_proj[i]
    h8 = simc_cdbfsi_xsec_proj[i]
    h9 = simc_cdbpwia_xsec_proj[i]
    
    if h0.title[16:18] == '35' or h0.title[16:18] == '45' or\
        h0.title[16:18] == '75':
        B.pl.figure(figsize=(8,6))     
        h2.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='s',color=cmp['tab20'](2),label='JML FSI')
        h3.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='D',color=cmp['tab20'](3),label='JML PWIA')
        h4.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='s',color=cmp['tab20'](4),label='MS Paris FSI')
        h5.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='D',color=cmp['tab20'](5),label='MS Paris PWIA')
        h6.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='s',color=cmp['tab20'](6),label='MS V18 FSI')
        h7.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='D',color=cmp['tab20'](7),label='MS V18 PWIA')
        h8.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='s',color=cmp['tab20'](8),label='MS CD-Bonn FSI')
        h9.plot_exp(markersize=4,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='D',color=cmp['tab20'](9),label='MS CD-Bonn PWIA')
        if plot_setting == 'pm_120':
            h3 = pm80_xsec_proj[i]
            h3.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='^',color=cmp['Set1'].colors[0],label='pm_80 Comm Data')
        elif plot_setting == 'pm_580':
            h3 = pm580_comm_xsec_proj[i]
            h3.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='^',color=cmp['Set1'].colors[0],label='pm_580 Comm Data')
        elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
            h3 = pm580_comm_xsec_proj[i]
            h4 = pm750_xsec_proj[i]
            h3.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='^',color=cmp['Set1'].colors[0],label='pm_580 Comm Data')
            h4.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    marker='^',color=cmp['Set1'].colors[2],label='pm_750 Comm Data')
        
        h01.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
                    color='black',label='DATA')
        # h1.plot_exp(markersize=8,elinewidth=1.2,capsize=0,ignore_zeros=True,
        #             color='black',label=f'PWIA DATA {s}')
        
        B.pl.yscale('log')
        B.pl.ylim(1e-8,1e-2)
        B.pl.xlim(-0.1,2.0)
        B.pl.ylabel('$\\frac{d^5\sigma}{d\Omega_e\Omega_pdE\'}$ [$\mu b \cdot sr^2 GeV^{-1}$]')
        B.pl.xlabel('$P_m$ [Gev/c]')
        B.pl.legend()
        ax = B.pl.gca()
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        
B.pl.show()

#%% plot cross section ratios data/FSI, data/PWIA

        
xsec_proj = histos[f'{plot_setting}_xsec_proj']
JML_fsi_xsec_proj = histos[f'{plot_setting}_JML_fsi_xsec_proj']
JML_pwia_xsec_proj = histos[f'{plot_setting}_JML_pwia_xsec_proj']
# CDB_fsi_xsec_proj = histos[f'{plot_setting}_CDB_fsi_xsec_proj']
# CDB_pwia_xsec_proj = histos[f'{plot_setting}_CDB_pwia_xsec_proj']
if plot_setting == 'pm_120':
    pm80_xsec_proj = histos['pm80_xsec_proj']
elif plot_setting == 'pm_580':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
    pm750_xsec_proj = histos['pm750_xsec_proj']

rows = 4
cols = 5
fig, ax=B.pl.subplots(rows,cols,figsize=(15,15),constrained_layout=True)
B.pl.suptitle('Cross Section Ratios Projections')
first_plot = True
for i in range(len(xsec_proj)):
    ax=B.pl.subplot(rows,cols,i+1)
    
    bc = xsec_proj[i].bin_content/JML_fsi_xsec_proj[i].bin_content
    bce = xsec_proj[i].bin_error/JML_fsi_xsec_proj[i].bin_content
    bc2 = xsec_proj[i].bin_content/JML_pwia_xsec_proj[i].bin_content
    bce2 = xsec_proj[i].bin_error/JML_pwia_xsec_proj[i].bin_content
    
    if plot_setting == 'pm_120':
        bc3 = xsec_proj[i].bin_content/pm80_xsec_proj[i].bin_content
        bce3 = xsec_proj[i].bin_error/pm80_xsec_proj[i].bin_content
    elif plot_setting == 'pm_580':
        bc3 = xsec_proj[i].bin_content/pm580_comm_xsec_proj[i].bin_content
        bce3 = xsec_proj[i].bin_error/pm580_comm_xsec_proj[i].bin_content
    elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
        bc3 = xsec_proj[i].bin_content/pm580_comm_xsec_proj[i].bin_content
        bce3 = xsec_proj[i].bin_error/pm580_comm_xsec_proj[i].bin_content
        bc4 = xsec_proj[i].bin_content/pm750_xsec_proj[i].bin_content
        bce4 = xsec_proj[i].bin_error/pm750_xsec_proj[i].bin_content
    
    r1 = B.histo(bin_center=xsec_proj[i].bin_center,
                 bin_content=bc,
                 bin_error=bce) # DATA/FSI
    r2 = B.histo(bin_center=xsec_proj[i].bin_center,
                 bin_content=bc2,
                 bin_error=bce2) # DATA/PWIA
    #DATA/comm_data
    r3 = B.histo(bin_center=xsec_proj[i].bin_center,
                 bin_content=bc3,
                 bin_error=bce3)
    if plot_setting == 'pm_800' or plot_setting == 'pm900':
       r4 = B.histo(bin_center=xsec_proj[i].bin_center,
                    bin_content=bc4,
                    bin_error=bce4) 
       r4.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                   label='DATA/COMM_DATA')

    r1.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                label='DATA/JML_FSI')
    r2.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                label='DATA/JML_PWIA')
    r3.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                label='DATA/COMM_DATA')
    
    # B.pl.yscale('log')
    B.pl.title(xsec_proj[i].title)
    B.pl.ylim(0.5,2.5)
    B.pl.xlim(0,1.2)
    if first_plot:
        B.pl.legend()
        B.pl.ylabel('$\\frac{\sigma_{DATA}}{\sigma_{THEORY}}$')
        B.pl.xlabel('$P_m$')
        first_plot = False
    else:
        B.pl.ylabel('')
        B.pl.xlabel('')
    
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    B.pl.hlines([0.9,1.1,0.8,1.2],xmin=0,xmax=0.2,linestyles='--',
                color=['green','green','black','black'])

B.pl.show()
    
if save:
    B.pl.savefig(f'yield_n_xsec/{plot_setting}_plots_alleff/{plot_setting}_Xsec_Ratios_Proj.png')
    
#%%
xsec_proj = histos[f'{plot_setting}_FSI_xsec_proj']['set_2']
JML_fsi_xsec_proj = histos[f'{plot_setting}_JML_fsi_xsec_proj']
JML_pwia_xsec_proj = histos[f'{plot_setting}_JML_pwia_xsec_proj']
# CDB_fsi_xsec_proj = histos[f'{plot_setting}_CDB_fsi_xsec_proj']
# CDB_pwia_xsec_proj = histos[f'{plot_setting}_CDB_pwia_xsec_proj']
if plot_setting == 'pm_120':
    pm80_xsec_proj = histos['pm80_xsec_proj']
elif plot_setting == 'pm_580':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
    pm580_comm_xsec_proj = histos['pm580_comm_xsec_proj']
    pm750_xsec_proj = histos['pm750_xsec_proj']

first_plot = True
for i in range(len(xsec_proj)):
    if i == 4 or i == 7:
        B.pl.figure()
        
        bc = xsec_proj[i].bin_content/JML_fsi_xsec_proj[i].bin_content
        bce = xsec_proj[i].bin_error/JML_fsi_xsec_proj[i].bin_content
        bc2 = xsec_proj[i].bin_content/JML_pwia_xsec_proj[i].bin_content
        bce2 = xsec_proj[i].bin_error/JML_pwia_xsec_proj[i].bin_content
        
        if plot_setting == 'pm_120':
            bc3 = xsec_proj[i].bin_content/pm80_xsec_proj[i].bin_content
            bce3 = xsec_proj[i].bin_error/pm80_xsec_proj[i].bin_content
        elif plot_setting == 'pm_580':
            bc3 = xsec_proj[i].bin_content/pm580_comm_xsec_proj[i].bin_content
            bce3 = xsec_proj[i].bin_error/pm580_comm_xsec_proj[i].bin_content
        elif plot_setting == 'pm_800' or plot_setting == 'pm_900':
            bc3 = xsec_proj[i].bin_content/pm580_comm_xsec_proj[i].bin_content
            bce3 = xsec_proj[i].bin_error/pm580_comm_xsec_proj[i].bin_content
            bc4 = xsec_proj[i].bin_content/pm750_xsec_proj[i].bin_content
            bce4 = xsec_proj[i].bin_error/pm750_xsec_proj[i].bin_content
        
        r1 = B.histo(bin_center=xsec_proj[i].bin_center,
                     bin_content=bc,
                     bin_error=bce) # DATA/FSI
        r2 = B.histo(bin_center=xsec_proj[i].bin_center,
                     bin_content=bc2,
                     bin_error=bce2) # DATA/PWIA
        #DATA/comm_data
        r3 = B.histo(bin_center=xsec_proj[i].bin_center,
                     bin_content=bc3,
                     bin_error=bce3)
        if plot_setting == 'pm_800' or plot_setting == 'pm_900':
           r4 = B.histo(bin_center=xsec_proj[i].bin_center,
                        bin_content=bc4,
                        bin_error=bce4) 
           r4.plot_exp(markersize=5,elinewidth=1.2,capsize=0,ignore_zeros=True,
                       label='DATA/750COMM_DATA')
    
        r1.plot_exp(markersize=8,elinewidth=2.4,capsize=0,ignore_zeros=True,
                    label='DATA/JML_FSI')
        r2.plot_exp(markersize=8,elinewidth=2.4,capsize=0,ignore_zeros=True,
                    label='DATA/JML_PWIA')
        r3.plot_exp(markersize=8,elinewidth=2.4,capsize=0,ignore_zeros=True,
                    label='DATA/COMM_DATA')
        
        # B.pl.yscale('log')
        B.pl.title(xsec_proj[i].title)
        B.pl.ylim(0,2.5)
        B.pl.xlim(0.,1.2)
    
        B.pl.legend()
        B.pl.ylabel('$\\frac{\sigma_{DATA}}{\sigma_{THEORY}}$')
        B.pl.xlabel('$P_m$')
        first_plot = False
        
        ax = B.pl.gca()
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        B.pl.hlines([0.9,1.1,0.8,1.2],xmin=0,xmax=1.2,linestyles='--',
                    color=['green','green','black','black'])

B.pl.show()