#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 11:11:23 2025

@author: gvill

SHMS Momentum alignment from Emiss:
    This code calculates the difference between SIMC and data
    on the Emiss spctrum. This difference is used to calculate
    an offset on the electron (SHMS) momentum.
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import LT.box as B

#import deut_analysis_tools as dat
import database_operations as db
import data_init as D
import cut_handler as C

#%% helper functions
def fit_legend(fit_params):
    ax = B.pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    for par in fit_params:
        fit_info = (
               f"A = {par['A'].value:.2f} $\pm$ {par['A'].err:.2f}\n"
               f"$\mu$ = {par['mean'].value:.4e} $\pm$ {par['mean'].err:.4e}\n"
               f"$\sigma$ = {par['sigma'].value:.4e} $\pm$ {par['sigma'].err:.4e}\n"
           )
        handles.append(mlines.Line2D([], [], color='None', label=fit_info))
        labels.append(fit_info)
    
    #ax.legend(handles=handles, labels=labels)
    return (handles, labels)

#%% file locations
deut_db = '/home/gvill/deuteron/deuteron_db/deuteron_db.db' 

#%% select hcana variables to load

my_branches = ['H.kin.secondary.emiss','H.gtr.p','P.gtr.p','H.gtr.dp',
               'P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph',
               'P.cal.etottracknorm']

my_branches_sim = ['Em','e_pf','h_pf','e_delta','h_delta','h_xptar',
                   'h_yptar','e_xptar','e_yptar','Weight','Normfac']
                   
#%% load data root files
# make sure files are where they should be!
RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                  select_branches=my_branches)
# load simc root files
# make sure files are where they should be!
SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                   select_branches=my_branches_sim)

#%% Declare constants
Eb = 10.542         # GeV
Mp = 0.938272       # GeV

#%% assign variables
# run list
runs = RUN.many
settings = SIMC.many

Em = {}
kf_meas = {}
Pf_meas = {}

for r in runs:
    # reconstructed missing energy in GeV
    Em[r] = RUN.Branches[r]['H.kin.secondary.emiss']
    
    # measured electron final momentum in GeV
    kf_meas[r] = RUN.Branches[r]['P.gtr.p']
    
    # measured proton final momentum in GeV
    Pf_meas[r] = RUN.Branches[r]['H.gtr.p']

sEm = {}
skf_meas = {}
sPf_meas = {}

for s in settings:
    sEm[s] = SIMC.Branches[s]['Em']
    skf_meas[s] = SIMC.Branches[s]['e_pf']*(1/1e3)
    sPf_meas[s] = SIMC.Branches[s]['h_pf']*(1/1e3)
    
#%% calculated quantiies

Ef = {}
Em_calc = {}

for r in runs:
    pf = Pf_meas[r]
    kf = kf_meas[r]
    
    Ef[r] = np.sqrt(Mp*Mp + pf*pf)
    Em_calc[r] = Eb + Mp - Ef[r] - kf
    
sEf = {}
sEm_calc = {}

for s in settings:
    pf = sPf_meas[s]
    kf = skf_meas[s]
    
    sEf[s] = np.sqrt(Mp*Mp + pf*pf)
    sEm_calc[s] = Eb + Mp - sEf[s] - kf

#%% calculate differences

dEm = {}
for r in runs:
    dEm[r] = Em_calc[r] - Em[r]

sdEm = {}    
for s in settings:    
    sdEm[s] = sEm_calc[s] - sEm[s]

#%% get norm and weights
NORM = D.get_norm(20851)
WEIGHTS = D.calc_weights(SIMC.Branches)

#%% make histos with no cuts

Em_histo = NORM*B.histo(Em,range=(-0.1,0.5),bins=100)
dEm_histo = NORM*B.histo(dEm,range=(-0.002,0.002),bins=100)

sEm_histo = B.histo(sEm,range=(-0.1,0.5),bins=100,weights=WEIGHTS,
                     calc_w2=True)
sdEm_histo = B.histo(sdEm,range=(-0.002,0.002),bins=100,weights=WEIGHTS,
                     calc_w2=True)

h_to_plot = {'$\Delta E_m$':dEm_histo,'$E_m$':Em_histo}

sh_to_plot = {'$\Delta E_m$':sdEm_histo,'$E_m$':sEm_histo}

for p in h_to_plot:
    h = h_to_plot[p]
    hS = sh_to_plot[p]
    B.pl.figure()
    h.plot()
    hS.plot()
    
    B.pl.title(p)
    B.pl.xlabel('')
    B.pl.ylabel('')
#%% define cuts
   
cuts_list = C.acceptance_cuts + [C.shms_calPID]

cuts_to_apply = []
for cut in cuts_list:
    br = RUN.Branches[C.HCANA_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply.append(cut_array)

all_cuts = cuts_to_apply[0]
for arr in cuts_to_apply:    
    all_cuts = all_cuts & arr    

#%% define cuts for SIMC

cuts_list_sim = C.acceptance_cuts

cuts_to_apply_sim = []
for cut in cuts_list_sim:
    cut.init()
    br = SIMC.Branches[C.SIMC_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply_sim.append(cut_array)

all_cuts_sim = cuts_to_apply_sim[0]
for arr in cuts_to_apply_sim:    
    all_cuts_sim = all_cuts_sim & arr 

#%% apply cuts
Em = Em[all_cuts]
dEm = dEm[all_cuts]

sEm = sEm[all_cuts_sim]
sdEm = sdEm[all_cuts_sim]
WEIGHTS = WEIGHTS[all_cuts_sim]

#%% make histos

Em_histo = NORM*B.histo(Em,range=(-0.05,0.05),bins=100)
dEm_histo = NORM*B.histo(dEm,range=(-0.002,0.002),bins=100)

sEm_histo = B.histo(sEm,range=(-0.05,0.05),bins=100,weights=WEIGHTS,
                     calc_w2=True)
sdEm_histo = B.histo(sdEm,range=(-0.002,0.002),bins=100,weights=WEIGHTS,
                     calc_w2=True)

h_to_plot = {'$E_m$':Em_histo}

sh_to_plot = {'$E_m$':sEm_histo}

#%% plot and fit

for plot in h_to_plot:

    B.pl.figure()
    s = RUN.setting
    h = h_to_plot[plot]
    hsim = sh_to_plot[plot]
    
    hsim.plot_exp(c='#ed6792',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='SIMC')
    hsim.fit(-0.015,0.03,plot_fit=False)
    hsim.fit(hsim.mean.value-hsim.sigma.value,hsim.mean.value+hsim.sigma.value)
    hsim.plot_fit(color='#de425b')
    
    mu_s = hsim.mean.value
    
    h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='Data') 
    h.fit(-0.035,0.01,plot_fit=False)
    h.fit(h.mean.value-h.sigma.value,h.mean.value+h.sigma.value)
    h.plot_fit(color='#1f77b4')
    
    mu = h.mean.value
    
    # if mu_s > mu:  
    #     diffs = mu_s - mu
    # else:
    #     diffs = mu - mu_s
    # diffs is SIMC-DATA
    diffs = mu_s - mu
    
    han, lab = fit_legend([h.fit_par,hsim.fit_par])
    
    #print(han,lab)
    
    diff_lab = f'SIMC - DATA = {diffs:.3e}'
    diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
    
    myhan = [han[1],han[2],han[0],han[3],diff_han]
    mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
    B.pl.legend(handles=myhan,labels=mylab,loc='best',fontsize=12)
    
    bot,top = plt.ylim()
    
    ### plot formatting
    # B.pl.title(f'$\Delta th_p$ for {s}')
    # B.pl.xlabel('')
    # B.pl.ylabel('')
    # B.pl.ylim(0,top)
    # #B.pl.legend(handles=han,labels=lab,loc='best')
    # B.pl.xlabel('[rad]')
    B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                colors=['#de425b','#1f77b4'],linestyles='--')
    B.pl.title(f'{plot} {RUN.many}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    ax = B.pl.gca()
    ax.set_autoscaley_on(True)    


