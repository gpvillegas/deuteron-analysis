#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 11:00:22 2026

@author: gvill
"""

import data_init as D
import cut_handler as C
import numpy as np
import LT.box as B
import matplotlib as mpl

mpl.rcParams["font.family"] = 'sans-serif'
mpl.rcParams["font.size"] = 14

T_branches = ['P.gtr.beta','H.gtr.beta','H.dc.1u1.time','H.dc.1u1.wirenum',
              'P.dc.1u1.time','P.dc.1u1.wirenum','P.cal.etottracknorm',
              'H.cal.etottracknorm']

R = 20873

ROOTfiles_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

RUN = D.DATA_INIT(data_type='deut23_data',run=R,
                         select_branches={'T':T_branches},
                         select_trees=['T'],
                         ROOTfiles_path=ROOTfiles_DIR) 
#%%
#to get rid of nested arrays of 1 element
def get_array(arr):
    #arr = data.Branches[r][time_name]
    newarr = np.zeros(arr.size)
    for i in range(arr.size):
        value = arr[i].tolist()
        try:    
            newarr[i] = value[0]
        except IndexError:
            newarr[i] = 0.
    return newarr

#%% plot hopdoscope beta
hbeta = RUN.Branches['H.gtr.beta']
pbeta = RUN.Branches['P.gtr.beta']

h_hbeta = B.histo(hbeta,range=(0.1,1.5),bins=100,
                  title='HMS Hodoscope Beta with Tracking',
                  xlabel='HMS $\\beta$',
                  ylabel='Counts')
h_hbeta.fit(0.87,0.97,plot_fit=False)

h_pbeta = B.histo(pbeta,range=(0.1,1.5),bins=100,
                  title='SHMS Hodoscope Beta with Tracking',
                  xlabel='SHMS $\\beta$',
                  ylabel='')
h_pbeta.fit(0.9,1.1,plot_fit=False)

B.pl.figure(figsize=(16,5))
B.pl.subplot(1,2,1)
h_hbeta.plot(filled=False)
h_hbeta.plot_fit(color='r',
                 label=f'A = {h_hbeta.A.value:0.3e}$\pm${h_hbeta.A.err:0.3e}\n'+\
                         f'Mean = {h_hbeta.mean.value:0.4f}$\pm${h_hbeta.mean.err:0.4f}\n'+\
                             f'Sigma = {h_hbeta.sigma.value:0.4f}$\pm${h_hbeta.sigma.err:0.4f}')
bp = 2.2622/np.sqrt(0.938**2 + 2.2622**2)
bd = 2.2622/np.sqrt(1.876**2 + 2.2622**2)
bt = 2.2622/np.sqrt(2.809**2 + 2.2622**2)
B.pl.vlines([bp,bd,bt],0,[h_hbeta.A.value,4900,4400],linestyles='--',color='black')
B.pl.text(0.95,h_hbeta.A.value,'$\\beta_{proton}$')
B.pl.text(0.66,5000,'$\\beta_{deuteron}$')
B.pl.text(0.56,4500,'$\\beta_{triton}$')
B.pl.legend(loc='upper left',fontsize=12)

B.pl.subplot(1,2,2)
h_pbeta.plot(filled=False)
h_pbeta.plot_fit(color='r',
                 label=f'A = {h_pbeta.A.value:0.3e}$\pm${h_pbeta.A.err:0.3e}\n'+\
                         f'Mean = {h_pbeta.mean.value:0.4f}$\pm${h_pbeta.mean.err:0.4f}\n'+\
                             f'Sigma = {h_pbeta.sigma.value:0.4f}$\pm${h_pbeta.sigma.err:0.4f}')
B.pl.vlines(1.0,0,h_pbeta.A.value,linestyles='--',color='black')
B.pl.text(1.05,h_pbeta.A.value-100,'$\\beta_{electron}$')
B.pl.legend(loc='upper left',fontsize=12)

#%% plot plane 1u1 drift time vs wirenum
hwirenum = get_array(RUN.Branches['H.dc.1u1.wirenum'])
hdctime = get_array(RUN.Branches['H.dc.1u1.time'])
pwirenum = get_array(RUN.Branches['P.dc.1u1.wirenum']) 
pdctime = get_array(RUN.Branches['P.dc.1u1.time'])

h_dctimeWirenum = B.histo2d(hwirenum[hdctime>0],hdctime[hdctime>0],
                            range=[(0,100),(-50,350)],bins=100,logz=True,
                            title='HMS Drift Time vs. Wire Number: Plane 1u1',
                            xlabel='Wire Number',
                            ylabel='Drift Time [ns]')
p_dctimeWirenum = B.histo2d(pwirenum[pdctime>0],pdctime[pdctime>0],
                            range=[(0,100),(-50,350)],bins=100,logz=True,
                            title='SHMS Drift Time vs. Wire Number: Plane 1u1',
                            xlabel='Wire Number',
                            ylabel='')

B.pl.figure(figsize=(15,5),layout='constrained')
B.pl.subplot(1,2,1)
h_dctimeWirenum.plot(colormap=mpl.colormaps['viridis'])

B.pl.subplot(1,2,2)
p_dctimeWirenum.plot(colormap=mpl.colormaps['viridis'])

#%% calorimeter plots
hcaletottracknorm = RUN.Branches['H.cal.etottracknorm']
pcaletottracknorm = RUN.Branches['P.cal.etottracknorm']

h_hcal = B.histo(hcaletottracknorm, range=(0.5,1.5), bins=100,
                 title='HMS Normalized Energy Deposited',
                 xlabel='$E_{dep}/P$',
                 ylabel='Counts')
# h_hcal.fit(0.95,1.08,plot_fit=False)
h_pcal = B.histo(pcaletottracknorm, range=(0.5,1.5), bins=100,
                 title='SHMS Normalized Energy Deposited',
                 xlabel='$E_{dep}/P$',
                 ylabel='')
h_pcal.fit(0.95,1.08,plot_fit=False)

B.pl.figure(figsize=(15,5))
B.pl.subplot(1,2,1)
h_hcal.plot(filled=False)
# h_hcal.plot_fit(color='r',
#                 label=f'A = {h_hcal.A.value:0.3e}$\pm${h_hcal.A.err:0.3e}\n'+\
#                         f'Mean = {h_hcal.mean.value:0.4f}$\pm${h_hcal.mean.err:0.4f}\n'+\
#                             f'Sigma = {h_hcal.sigma.value:0.4f}$\pm${h_hcal.sigma.err:0.4f}')
# B.pl.legend(loc='upper left')

B.pl.subplot(1,2,2)
h_pcal.plot(filled=False)
h_pcal.plot_fit(color='r',
                label=f'A = {h_pcal.A.value:0.3e}$\pm${h_pcal.A.err:0.3e}\n'+\
                        f'Mean = {h_pcal.mean.value:0.4f}$\pm${h_pcal.mean.err:0.4f}\n'+\
                            f'Sigma = {h_pcal.sigma.value:0.4f}$\pm${h_pcal.sigma.err:0.4f}')
B.pl.legend(loc='upper left')

