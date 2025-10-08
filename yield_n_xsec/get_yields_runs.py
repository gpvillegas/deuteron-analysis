#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 21:16:43 2025

@author: gvill

Extract Deuteron Cross Sections (PRELIMINARY)
"""

import data_init as D
import LT.box as B
import cut_handler as C

import numpy as np
from matplotlib import colormaps as cmp

#%% conversion factors
radtodeg = 180./np.pi

#%% select appropriate HCANA variables 
kin_var = ['H.kin.secondary.th_bq','H.kin.secondary.pmiss',
           'H.kin.secondary.emiss_nuc','P.kin.primary.Q2']

coinTime_var = ['CTime.epCoinTime_ROC2']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

br_sel = kin_var + coinTime_var + acc_var + calPID_var

#%% load runs with selected branches
run_list = [20871,20872]

deep_pm120 = D.DATA_INIT(data_type='deut23_data',run=run_list,
                         select_branches=br_sel)

#%% load SIMC files

deep_pm120_SIMC = D.DATA_INIT(data_type='SIMC',setting = 'pm_120')

#%% get normalization factors for data

NORM = {}
for run in run_list:
    NORM[run] = D.get_norm(run)

#%% get SIMC weights

WEIGHT = D.calc_weights(deep_pm120_SIMC.Branches)
WEIGHT_PS = D.calc_weights_PS(deep_pm120_SIMC.Branches)

#%% declare variables to make histos
Pm = {}
Th_rq = {}
Em = {}

for run in run_list:
    Pm[run] = deep_pm120.Branches[run]['H.kin.secondary.pmiss']
    Th_rq[run] = deep_pm120.Branches[run]['H.kin.secondary.th_bq']*radtodeg
    Em[run] = deep_pm120.Branches[run]['H.kin.secondary.emiss_nuc']

Pm_SIMC = deep_pm120_SIMC.Branches['Pm']
Th_rq_SIMC = deep_pm120_SIMC.Branches['theta_rq']*radtodeg
Em_SIMC = deep_pm120_SIMC.Branches['Em']

#%% define cuts for data

## add Coin Time corrected variable for cut

for run in run_list:
    coinTime = deep_pm120.Branches[run]['CTime.epCoinTime_ROC2']
    coinTime_corr = D.get_coinTimeOffset(coinTime)

    deep_pm120.Branches[run].update({'CTime.epCoinTime_ROC2_corr':coinTime_corr})

# =============================================================================
# cuts_to_apply = {'H.gtr.dp':(-10.,10.),'P.gtr.dp':(-10.,22.),
#                  'H.gtr.th':(-0.08,0.08),'H.gtr.ph':(-0.03,0.03),
#                  'P.gtr.th':(-0.03,0.03),'P.gtr.ph':(-0.03,0.03),
#                 'P.cal.etottracknorm':(0.7,1.3),
#                 'CTime.epCoinTime_ROC2_corr':(-2.,2.),
#                 'P.kin.primary.Q2':(4.,6.),
#                 'H.kin.secondary.emiss_nuc':(-0.04,0.04)}
# 
# cuts = deep_pm120.CUT(cut_dict=cuts_to_apply,
#                       where_branches=deep_pm120.Branches[20871])
# =============================================================================

ct_cut = C.WCUT(xmin=-2.,xmax=2.,name='epCoinTime')
cal_cut = C.WCUT(xmin=0.7,xmax=1.3,name='calPID')
Q2_cut = C.WCUT(xmin=4.,xmax=np.infty,name='Q2')
em_cut = C.WCUT(xmin=-np.infty,xmax=0.04,name='emiss')

C.HCANA_names.update({'epCoinTime':'CTime.epCoinTime_ROC2_corr',
                      'calPID':'P.cal.etottracknorm',
                      'Q2':'P.kin.primary.Q2',
                      'emiss':'H.kin.secondary.emiss_nuc'})

cuts_list = C.acceptance_cuts + [ct_cut,cal_cut,Q2_cut,em_cut]

cuts_to_apply = {}
for run in run_list:
    c_list = []
    for cut in cuts_list:
        br = deep_pm120.Branches[run][C.HCANA_names[cut.name]] 
        #cut.stats()
        c_list.append(cut(br))  
        
    cuts_to_apply[run] = c_list

all_cuts = {}
for run in run_list:
    c_all = cuts_to_apply[run][0]
    for arr in cuts_to_apply[run]:    
        c_all = c_all & arr      
    all_cuts[run] = c_all    

#%% define cuts for SIMC
# =============================================================================
# cuts_to_apply_SIMC = {'h_delta':(-10.,10.),'e_delta':(-10.,22.),
#             'h_xptar':(-0.08,0.08),'h_yptar':(-0.03,0.03),
#             'e_xptar':(-0.03,0.03),'e_yptar':(-0.03,0.03),
#             'Q2':(4.,6.),'Em':(-0.04,0.04)}
# 
# cuts_SIMC = deep_pm120_SIMC.CUT(cut_dict=cuts_to_apply_SIMC,
#                                 where_branches=deep_pm120_SIMC.Branches)
# =============================================================================
C.SIMC_names.update({'Q2':'Q2','emiss':'Em'})

cuts_list = C.acceptance_cuts + [Q2_cut,em_cut]

cuts_to_apply_sim = []
for cut in cuts_list:
    br = deep_pm120_SIMC.Branches[C.SIMC_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply_sim.append(cut_array)

all_cuts_sim = cuts_to_apply_sim[0]
for arr in cuts_to_apply_sim:    
    all_cuts_sim = all_cuts_sim & arr  
#%% apply cuts
Pm_cut = {}
Th_rq_cut = {}
Em_cut = {}

for run in run_list:  
    Pm_cut[run] = Pm[run][all_cuts[run]]   
    Th_rq_cut[run] = Th_rq[run][all_cuts[run]]
    Em_cut[run] = Em[run][all_cuts[run]]

Pm_SIMC_cut = Pm_SIMC[all_cuts_sim]
Th_rq_SIMC_cut = Th_rq_SIMC[all_cuts_sim]
Em_SIMC_cut = Em_SIMC[all_cuts_sim]
WEIGHT_cut = WEIGHT[all_cuts_sim]

#%% make 2d histos
Pm_vs_Thrq_histo ={}
Pm_vs_Em_histo={}

for run in run_list:
    Pm_vs_Thrq_histo[run] = B.histo2d(Th_rq_cut[run],Pm_cut[run],
                                      range=[(0.,180.),(-0.1,1.0)],
                                      bins=100)
    Pm_vs_Em_histo[run] = B.histo2d(Em_cut[run],Pm_cut[run],
                                    range=[(-0.02,0.05),(-0.1,1.0)],
                                    bins=100)
    
    
Pm_vs_Thrq_SIMC_histo = B.histo2d(Th_rq_SIMC_cut,Pm_SIMC_cut,
                                  weights=WEIGHT_cut, calc_w2=True,
                                  range=[(0,180),(-0.1,1.0)],
                                  bins=100)

Pm_vs_Em_SIMC_histo = B.histo2d(Em_SIMC_cut,Pm_SIMC_cut,
                                  weights=WEIGHT_cut, calc_w2=True,
                                  range=[(-0.02,0.05),(-0.1,1.0)],
                                  bins=100)
#%% plot histos

histos_to_plot = {'$P_{m}$ vs. $\\theta_{rq}$':Pm_vs_Thrq_histo,
                  '$P_{m}$ vs. $E_{m}$':Pm_vs_Em_histo}
sim_histos_to_plot = {'$P_{m}$ vs. $\\theta_{rq}$ SIMC':Pm_vs_Thrq_SIMC_histo,
                      '$P_{m}$ vs. $E_{m}$ SIMC':Pm_vs_Em_SIMC_histo}
for run in run_list:
    for plot in histos_to_plot:
        B.pl.figure()
        h = histos_to_plot[plot][run]
        h.plot(colormap=cmp['viridis'])
    
        B.pl.title(f'{plot} run {run}')
        B.pl.xlabel('')
        B.pl.ylabel('')
        
for plot in sim_histos_to_plot:
    B.pl.figure()
    h = sim_histos_to_plot[plot]
    h.plot(colormap=cmp['viridis'])

    B.pl.title(plot)
    B.pl.xlabel('')
    B.pl.ylabel('')        

#%% set up histogram projections

def get_yield(h2,norm,bin_range,x_min,x_max):
    counts = []
    for nx in bin_range[:]:
        h = norm*h2.project_y(bins = [nx])
        
        counts.append(h.sum(x_min,x_max))
    
    counts = np.array(counts)
       
    y_values = h2.y_bin_center[bin_range]
    yield_values = counts[:,0]
    yield_errors = counts[:,1]
    
    return (y_values,yield_values,yield_errors)

def get_yield_runs(rnum,h2,norm,bin_range,x_min,x_max):
    counts = []
    
    first=True
    
    for nx in bin_range[:]:
        hl = []
        for r in rnum:
            h = norm[r]*h2[r].project_y(bins = [nx])
            hl.append(h)
        if first:
            hall = h
            first = False
        else:    
            hall = hall + h
        hall = hall/len(rnum)
        counts.append(hall.sum(x_min,x_max))
    
    counts = np.array(counts)
       
    y_values = hall.bin_center[bin_range]
    yield_values = counts[:,0]
    yield_errors = counts[:,1]
    
    return (y_values,yield_values,yield_errors)

def get_maxbin_loc(histo):
    h = histo
    
    pos = np.where(h.bin_content==h.bin_content.max())
     
    bin_x = pos[0][0]
    bin_y= pos[1][0]
    
    return bin_x, bin_y

x,y = get_maxbin_loc(Pm_vs_Em_histo[20871])    
brange = np.arange(x-15,x+10)

pm_values, pmEm_yield, pmEm_error =\
    get_yield_runs(rnum=run_list,h2=Pm_vs_Em_histo, norm=NORM, 
              bin_range=brange, x_min=-0.05, x_max=0.05)

x,y = get_maxbin_loc(Pm_vs_Em_SIMC_histo)    
brange = np.arange(x-15,x+10)
    
pm_values_sim, pmEm_yield_sim, pmEm_error_sim =\
    get_yield(h2=Pm_vs_Em_SIMC_histo, norm=1, 
              bin_range=brange, x_min=-0.05, x_max=0.05)
    
# B.plot_exp(c='#7b8fd4',marker='+',markersize=8,
#                      capsize=0,mew=1,elinewidth=1,label='Data')



