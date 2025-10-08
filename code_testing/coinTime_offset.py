#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 09:59:34 2025

@author: gvill

Get coin time offset for deuteron runs and store them in db table
"""
import numpy as np
import matplotlib.pyplot as plt

import LT.box as B
import database_operations as db
import data_init as D
import cut_handler2 as cu2

#%%
deut_db = 'deuteron_db.db' 

gtr = D.get_list(db.retrieve(deut_db, 'Branches', 'TTree_Branches', 
                            where='Type = \'gtr\''))


br_sel = gtr + ['CTime.epCoinTime_ROC2','P.cal.etottracknorm']

#%% load root files
heep = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                    select_branches=br_sel)
# Run = D.DATA_INIT(data_type='deut23_data', run=20851, 
#                    select_branches=br_sel)
#%%
CTime = Run.Branches['CTime.epCoinTime_ROC2']

CTime_histo = B.histo(CTime,range=(64,68),bins=200,
                          title = f'Run {Run.many} no cuts')
CTime_histo.fit(plot_fit=False)

B.pl.figure()
CTime_histo.plot()
CTime_histo.plot_fit()
#%% plot ctime histos no cuts
CTime = {}
CTime_histos = {}

for r in heep.many:
    CTime[r] = heep.Branches[r]['CTime.epCoinTime_ROC2']
    
    CTime_histos[r] = B.histo(CTime[r],range=(64,68),bins=200,
                              title = f'Run {r} no cuts')
    CTime_histos[r].fit(plot_fit=False)
    
    # B.pl.figure()
    # CTime_histos[r].plot()
    # CTime_histos[r].plot_fit()

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('CTime.epCoinTime_ROC2 (no cuts)', fontsize=16)
rows = 2
cols = 4
i = 1
for r in heep.many:
    s = heep.setting[r]
    ax = plt.subplot(rows,cols,i)
    h = CTime_histos[r]
    h.plot()
    h.plot_fit()
    i+=1
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'Run {r} {s}')

#%% define cuts for data
std_cuts = {'H.gtr.dp':(-10.,10.),'P.gtr.dp':(-10.,20.),
            'H.gtr.th':(-0.08,0.08),'H.gtr.ph':(-0.03,0.03),
            'P.gtr.th':(-0.03,0.03),'P.gtr.ph':(-0.03,0.03),
            'P.cal.etottracknorm':(0.7,10)}
#%%
my_cuts = []

for k in std_cuts:
    vmin,vmax = std_cuts[k]
    my_cuts.append(cu2.wcut(vmin, vmax, name = k)) 
     
#%%
cuts = {}
for r in heep.many:
    flag = True
    for vc in std_cuts:
        cut_min = heep.Branches[r][vc] >= std_cuts[vc][0] 
        cut_max = heep.Branches[r][vc] <= std_cuts[vc][1]
        this_cut = cut_min & cut_max
        this_cut = np.array(this_cut)
        
        if flag:
            cuts[r] = this_cut
            flag = False
        else:
            cuts[r] = cuts[r] & this_cut  


# cuts = heep.CUT(iterator = None, cuts = my_cuts, where_branches = heep.Branches[20851])

#%% apply cuts

CTime_c = {}
CTime_histos_wcuts = {}
for r in heep.many:
    CTime_c[r] = CTime[r][cuts[r]]
    
    CTime_histos_wcuts[r] = B.histo(CTime_c[r],range=(64,68),bins=200,
                                    title=f'Run {r}') 
    
    
    CTime_histos_wcuts[r].fit(plot_fit=False)
    
    # B.pl.figure()
    # CTime_histos_wcuts[r].plot()
    # CTime_histos_wcuts[r].plot_fit()
    
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('CTime.epCoinTime_ROC2 (with cuts)', fontsize=16)
rows = 2
cols = 4
i = 1
for r in heep.many:
    s = heep.setting[r]
    ax = plt.subplot(rows,cols,i)
    h = CTime_histos_wcuts[r]
    h.plot()
    h.plot_fit()
    i+=1
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'Run {r} {s}')    
    
#%%

mu = {}

for r in heep.many:
    mu[r] = CTime_histos_wcuts[r].mean.value 
    

CTime_histos_wcuts_woffset = {}
for r in heep.many:
    CTime_corr = CTime_c[r] - mu[r]
    CTime_histos_wcuts_woffset[r] = B.histo(CTime_corr,range=(-5,5),bins=200,
                                    title=f'Run {r}') 
    
    
    CTime_histos_wcuts_woffset[r].fit(plot_fit=False)
    
    # B.pl.figure()
    # CTime_histos_wcuts_woffset[r].plot()
    #CTime_histos_wcuts_woffset[r].plot_fit()    
    
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('CTime.epCoinTime_ROC2 (w/cuts w/ offset)', fontsize=16)
rows = 2
cols = 4
i = 1
for r in heep.many:
    s = heep.setting[r]
    ax = plt.subplot(rows,cols,i)
    h = CTime_histos_wcuts_woffset[r]
    h.plot()
    h.plot_fit()
    i+=1
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'Run {r} {s}')  

#%%

def get_coinTimeOffset(coinTime):

    CTime = coinTime
    
    CTime_histo_0 = B.histo(CTime,range=(-100,100),bins=200,
                              title = 'CTime no cuts')
    CTime_histo_0.fit(plot_fit=False)
    
    CTime_histo_0.plot()
    CTime_histo_0.plot_fit()
    
    mu_0 = CTime_histo_0.mean.value
    
    sigma_0 = CTime_histo_0.sigma.value
    
    x_min = -mu_0 + 2.*sigma_0
    x_max = mu_0 + 2.*sigma_0
    
    CTime_histo = B.histo(CTime,range=(x_min,x_max),bins=200,
                              title = 'CTime no cuts')
    
    CTime_histo.fit(plot_fit=False)
    
    B.pl.figure()
    CTime_histo.plot()
    CTime_histo.plot_fit()
    
    mu = CTime_histo.mean.value
    
    CTime_corr = CTime - mu
    
    return CTime_corr
    