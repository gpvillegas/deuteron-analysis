#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 10:14:41 2025

@author: gvill

Create correlation plots for optics variables:
    xtar, ytar, xptar, yptar,
    xfp, yfp, xpfp, ypfp,
    W, Em, dp
"""

import data_init as D
import database_operations as db
import LT.box as B
import cut_handler as C

#import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cmp
from matplotlib.transforms import Bbox
import os

#%% helper functions
def full_extent(ax,padw=0.,padh=0.,tx=0.,ty=0.):
    """Get the full extent of an axes, including axes labels, tick labels, and
    titles."""
    # For text objects, we need to draw the figure first, otherwise the extents
    # are undefined.
    ax.figure.canvas.draw()
    items = ax.get_xticklabels() + ax.get_yticklabels()
    items += [ax, ax.title, ax.xaxis.label, ax.yaxis.label]
    # items += [ax, ax.title]
    bbox = Bbox.union([item.get_window_extent() for item in items])
    bbox1 = bbox.expanded(1.0 + padw, 1.0 + padh)
    return bbox1.translated(tx,ty) 

# %%

runs = [20840,20841,20846,20851,20858,20861,20868,20869] #heep_coin runs
# runs = [20840,20851,20868,20869]
# runs = [3283,3284,3285,3286,3287]
settings = ['delta_scan_0','delta_scan_+12']
deut_db_name = 'deuteron_db.db'

fp_branches = D.get_list(db.retrieve(deut_db_name, 'Branches', 'TTree_Branches',
                                     where='Branches like \'%.dc.%fp\''))

tar_branches = D.get_list(db.retrieve(deut_db_name, 'Branches', 'TTree_Branches',
                          where='Type = \'gtr\''))

kin_branches = ['H.kin.secondary.emiss', 'P.kin.primary.W']

sieve_branches = ['H.extcor.xsieve', 'H.extcor.ysieve',
                  'P.extcor.xsieve', 'P.extcor.ysieve']

selection = tar_branches + kin_branches + fp_branches + sieve_branches

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','W','Em','e_xptar',
               'e_yptar','h_xptar','h_yptar','h_xfp','h_xpfp','h_yfp','h_ypfp',
               'e_xfp','e_xpfp','e_yfp','e_ypfp']

# %%
# data = D.DATA_INIT(data_type='deut23_data', run=runs,
#                    select_branches=selection)

data = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin',
                   select_branches={'T':selection})


# %%
# simc = D.DATA_INIT(data_type='SIMC', setting=settings)

simc = D.DATA_INIT(data_type='SIMC', kin_study='heep_coin',
                   select_branches={'SNT':br_sel_SIMC},simc_type='-')

#%% define cuts for data
cuts_list = C.acceptance_cuts 

cuts_to_apply = {}
for r in data.many:
    c_list = []
    print(f'Cuts for Run {r}')
    for cut in cuts_list:
        cut.init()
        br = data.Branches[r][C.HCANA_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
        
    # add collimator cut
    hxc = data.Branches[r]['H.extcor.xsieve']
    hyc = data.Branches[r]['H.extcor.ysieve']
    
    hcoll_cut = C.coll_cut(hxc, hyc, spec='HMS')
    c_list.append(hcoll_cut)   
    
    cuts_to_apply[r] = c_list
    
all_cuts = {}
for r in data.many:
    all_cuts_arr = cuts_to_apply[r][0]
    for arr in cuts_to_apply[r]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts[r] = all_cuts_arr
#%% define cuts for simc
cuts_to_apply_sim = {}
for s in simc.many:
    c_list = []
    print(f'Cuts for SIMC setting {s}')
    for cut in cuts_list:
        cut.init()
        br = simc.Branches[s][C.SIMC_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
        
    cuts_to_apply_sim[s] = c_list

all_cuts_sim = {}
for s in simc.many:
    all_cuts_arr = cuts_to_apply_sim[s][0]
    for arr in cuts_to_apply_sim[s]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts_sim[s] = all_cuts_arr
    
# %% cuts single run
cuts_list = C.acceptance_cuts

cuts_to_apply = []
for cut in cuts_list:
    br = data.Branches[20851][C.HCANA_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply.append(cut_array)

all_cuts = cuts_to_apply[0]
for arr in cuts_to_apply:    
    all_cuts = all_cuts & arr  
    
cuts_to_apply_sim = []
for cut in cuts_list:
    cut.init()
    br = simc.Branches['delta_scan_0'][C.SIMC_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply_sim.append(cut_array)

all_cuts_sim = cuts_to_apply_sim[0]
for arr in cuts_to_apply_sim:    
    all_cuts_sim = all_cuts_sim & arr 


# %%
# HMS
hxptar = {}
hyptar = {}
hxfp = {}
hxpfp = {}
hyfp = {}
hypfp = {}
hdp = {}
Emiss = {}

#SHMS
pxtar = {}
pytar = {}
pxptar = {}
pyptar = {}
pxfp = {}
pxpfp = {}
pyfp = {}
pypfp = {}
pdp = {}
W = {}

for r in data.many:
    #HMS
    # target variables
    hxptar[r] = data.Branches[r]['H.gtr.th']
    hyptar[r] = data.Branches[r]['H.gtr.ph']
    # focal plane variables
    hxfp[r] = data.Branches[r]['H.dc.x_fp']
    hxpfp[r] = data.Branches[r]['H.dc.xp_fp']
    hyfp[r] = data.Branches[r]['H.dc.y_fp']
    hypfp[r]= data.Branches[r]['H.dc.yp_fp']
    # delta
    hdp[r] = data.Branches[r]['H.gtr.dp']
    # kinematic variables
    Emiss[r] = data.Branches[r]['H.kin.secondary.emiss']

    # SHMS
    # target variables
    pxtar[r] = data.Branches[r]['P.gtr.x']
    pytar[r] = data.Branches[r]['P.gtr.y']
    pxptar[r] = data.Branches[r]['P.gtr.th']
    pyptar[r] = data.Branches[r]['P.gtr.ph']
    # focal plane variables
    pxfp[r] = data.Branches[r]['P.dc.x_fp']
    pxpfp[r] = data.Branches[r]['P.dc.xp_fp']
    pyfp[r]= data.Branches[r]['P.dc.y_fp']
    pypfp[r] = data.Branches[r]['P.dc.yp_fp']
    # delta
    pdp[r] = data.Branches[r]['P.gtr.dp']
    # kinematic variables
    W[r]= data.Branches[r]['P.kin.primary.W']

data_var_hx = {'hxptar':[hxptar,(-0.08,0.08)],
               'hyptar':[hyptar,(-0.03,0.03)],
               'hxfp':[hxfp,(-30.,30.)],
               'hxpfp':[hxpfp,(-0.05,0.05)],
               'hyfp':[hyfp,(-10.,20.)],
               'hypfp':[hypfp,(-0.02,0.02)]}

data_var_px = {'pxptar':[pxptar,(-0.045,0.045)],
               'pyptar':[pyptar,(-0.04,0.04)],
               'pxfp':[pxfp,(-20.,22.)],
               'pxpfp':[pxpfp,(-0.05,0.05)],
               'pyfp':[pyfp,(-10.,10.)],
               'pypfp':[pypfp,(-0.02,0.02)]}

data_var_y = {'hdp':[hdp,(-15.,15.)],
                  'Emiss':[Emiss,(-0.1,0.1)],
                      'pdp':[pdp,(-22.,22.)],
                          'W':[W,(0.5,1.5)]}

# data_var = [hxptar,hyptar,hdp,Emiss,pxptar,pyptar,pdp,W]

#%% SIMC variables go here


# HMS SIMC
shxptar = {}
shyptar = {}
shxfp = {}
shxpfp = {}
shyfp = {}
shypfp = {}
shdp = {}
sEmiss = {}

#SHMS SIMC
spxptar = {}
spyptar = {}
spxfp = {}
spxpfp = {}
spyfp = {}
spypfp = {}
spdp = {}
sW = {}

for s in simc.many:
    # HMS SIMC
    # target variables
    shxptar[s] = simc.Branches[s]['h_xptar']
    shyptar[s] = simc.Branches[s]['h_yptar']
    # focal plane variables
    shxfp[s] = simc.Branches[s]['h_xfp']
    shxpfp[s] = simc.Branches[s]['h_xpfp']
    shyfp[s] = simc.Branches[s]['h_yfp']
    shypfp[s] = simc.Branches[s]['h_ypfp']
    # delta
    shdp[s] = simc.Branches[s]['h_delta']
    # kinematic variables
    sEmiss[s] = simc.Branches[s]['Em']
    
    # SHMS SIMC
    # target variables
    spxptar[s] = simc.Branches[s]['e_xptar']
    spyptar[s] = simc.Branches[s]['e_yptar']
    # focal plane variables
    spxfp[s] = simc.Branches[s]['e_xfp']
    spxpfp[s] = simc.Branches[s]['e_xpfp']
    spyfp[s] = simc.Branches[s]['e_yfp']
    spypfp[s] = simc.Branches[s]['e_ypfp']
    # delta
    spdp[s] = simc.Branches[s]['e_delta']
    # kinematic variables
    sW[s] = simc.Branches[s]['W']

simc_var_hx = {'hxptar':[shxptar,(-0.08,0.08)],
               'hyptar':[shyptar,(-0.03,0.03)],
               'hxfp':[shxfp,(-30.,30.)],
               'hxpfp':[shxpfp,(-0.05,0.05)],
               'hyfp':[shyfp,(-10.,20.)],
               'hypfp':[shypfp,(-0.02,0.02)]}

simc_var_px = {'pxptar':[spxptar,(-0.03,0.03)],
               'pyptar':[spyptar,(-0.015,0.015)],
               'pxfp':[spxfp,(-30.,30.)],
               'pxpfp':[spxpfp,(-0.05,0.05)],
               'pyfp':[spyfp,(-10.,20.)],
               'pypfp':[spypfp,(-0.02,0.02)]}

simc_var_y = {'hdp':[shdp,(-10.,10.)],
                  'Emiss':[sEmiss,(-0.1,0.1)],
                      'pdp':[spdp,(-22.,22.)],
                          'W':[sW,(0.5,1.5)]}


# %% assign variables single run
# HMS
# target variables
hxtar = data.Branches[20851]['H.gtr.x']
hxptar = data.Branches[20851]['H.gtr.th']
hytar = data.Branches[20851]['H.gtr.y']
hyptar = data.Branches[20851]['H.gtr.ph']
# focal plane variables
hxfp = data.Branches[20851]['H.dc.x_fp']
hxpfp = data.Branches[20851]['H.dc.xp_fp']
hyfp = data.Branches[20851]['H.dc.y_fp']
hypfp = data.Branches[20851]['H.dc.yp_fp']
# delta
hdp = data.Branches[20851]['H.gtr.dp']
# kinematic variables
Emiss = data.Branches[20851]['H.kin.secondary.emiss']

# HMS SIMC
# target variables
shxptar = simc.Branches['delta_scan_0']['h_xptar']
shytar = simc.Branches['delta_scan_0']['h_ytar']
shyptar = simc.Branches['delta_scan_0']['h_yptar']
# focal plane variables
shxfp = simc.Branches['delta_scan_0']['h_xfp']
shxpfp = simc.Branches['delta_scan_0']['h_xpfp']
shyfp = simc.Branches['delta_scan_0']['h_yfp']
shypfp = simc.Branches['delta_scan_0']['h_ypfp']
# delta
shdp = simc.Branches['delta_scan_0']['h_delta']
# kinematic variables
sEmiss = simc.Branches['delta_scan_0']['Em']

# SHMS
# target variables
pxtar = data.Branches[20851]['P.gtr.x']
pxptar = data.Branches[20851]['P.gtr.th']
pytar = data.Branches[20851]['P.gtr.y']
pyptar = data.Branches[20851]['P.gtr.ph']
# focal plane variables
pxfp = data.Branches[20851]['P.dc.x_fp']
pxpfp = data.Branches[20851]['P.dc.xp_fp']
pyfp = data.Branches[20851]['P.dc.y_fp']
pypfp = data.Branches[20851]['P.dc.yp_fp']
# delta
pdp = data.Branches[20851]['P.gtr.dp']
# kinematic variables
W = data.Branches[20851]['P.kin.primary.W']

# SHMS SIMC
# target variables
spxptar = simc.Branches['delta_scan_0']['e_xptar']
spytar = simc.Branches['delta_scan_0']['e_ytar']
spyptar = simc.Branches['delta_scan_0']['e_yptar']
# focal plane variables
spxfp = simc.Branches['delta_scan_0']['e_xfp']
spxpfp = simc.Branches['delta_scan_0']['e_xpfp']
spyfp = simc.Branches['delta_scan_0']['e_yfp']
spypfp = simc.Branches['delta_scan_0']['e_ypfp']
# delta
spdp = simc.Branches['delta_scan_0']['e_delta']
# kinematic variables
sW = simc.Branches['delta_scan_0']['W']
#%% apply cuts single run
# HMS
# target variables
hxptar = hxptar[all_cuts]
hytar = hytar[all_cuts]
hyptar = hyptar[all_cuts]
# focal plane variables
hxfp = hxfp[all_cuts]
hxpfp = hxpfp[all_cuts]
hyfp = hyfp[all_cuts]
hypfp = hypfp[all_cuts]
# delta
hdp = hdp[all_cuts]
# kinematic variables
Emiss = Emiss[all_cuts]

# HMS SIMC
# target variables
shxptar = shxptar[all_cuts_sim]
shytar = shytar[all_cuts_sim]
shyptar = shyptar[all_cuts_sim]
# focal plane variables
shxfp = shxfp[all_cuts_sim]
shxpfp = shxpfp[all_cuts_sim]
shyfp = shyfp[all_cuts_sim]
shypfp = shypfp[all_cuts_sim]
# delta
shdp = shdp[all_cuts_sim]
# kinematic variables
sEmiss = sEmiss[all_cuts_sim]

# SHMS
# target variables
pxptar = pxptar[all_cuts]
pytar = pytar[all_cuts]
pyptar = pyptar[all_cuts]
# focal plane variables
pxfp = pxfp[all_cuts]
pxpfp = pxpfp[all_cuts]
pyfp = pyfp[all_cuts]
pypfp = pypfp[all_cuts]
# delta
pdp = pdp[all_cuts]
# kinematic variables
W = W[all_cuts]

# SHMS SIMC
# target variables
spxptar = spxptar[all_cuts_sim]
spytar = spytar[all_cuts_sim]
spyptar = spyptar[all_cuts_sim]
# focal plane variables
spxfp = spxfp[all_cuts_sim]
spxpfp = spxpfp[all_cuts_sim]
spyfp = spyfp[all_cuts_sim]
spypfp = spypfp[all_cuts_sim]
# delta
spdp = spdp[all_cuts_sim]
# kinematic variables
sW = sW[all_cuts_sim]

#%% apply cuts
for y in data_var_y:
    for r in runs:
        var = data_var_y[y][0][r]
        var = var[all_cuts[r]]
        data_var_y[y][0][r] = var
        
for x in data_var_hx:
    for r in runs:
        var = data_var_hx[x][0][r]
        var = var[all_cuts[r]]
        data_var_hx[x][0][r] = var
        
for x in data_var_px:
    for r in runs:
        var = data_var_px[x][0][r]
        var = var[all_cuts[r]]       
        data_var_px[x][0][r] = var

for y in simc_var_y:
    for s in settings:
        var = simc_var_y[y][0][s]
        var = var[all_cuts_sim[s]]
        simc_var_y[y][0][s] = var
        
for x in simc_var_hx:
    for s in settings:
        var = simc_var_hx[x][0][s]
        var = var[all_cuts_sim[s]]
        simc_var_hx[x][0][s] = var
        
for x in simc_var_px:
    for s in settings:
        var = simc_var_px[x][0][s]
        var = var[all_cuts_sim[s]]
        simc_var_px[x][0][s] = var               

#%% make data histos with or without cuts
 
data_hist_to_plot = {}
for y in data_var_y:
    hy = data_var_y[y][0]
    ylim = data_var_y[y][1]
    for x in data_var_hx:
        hx = data_var_hx[x][0]
        xlim = data_var_hx[x][1]
        h = {}
        
        if x[0] == y[0] or y == 'Emiss' or y == 'W':
            for r in hx:
                    h[r] = B.histo2d(hx[r],hy[r],
                            range=[xlim,ylim],bins=100,logz=True,
                                xlabel=x,ylabel=y,title=f'{y}_vs_{x} Run {r}')           
            data_hist_to_plot[f'HMS {y}_vs_{x}'] = h           
        
    for x in data_var_px:
        hx = data_var_px[x][0]
        xlim = data_var_px[x][1]
        h = {}
        
        if x[0] == y[0] or y == 'Emiss' or y == 'W':
            for r in hx:
                    h[r] = B.histo2d(hx[r],hy[r],
                            range=[xlim,ylim],bins=100,logz=True,
                            xlabel=x,ylabel=y,title=f'{y}_vs_{x} Run {r}')
            
            data_hist_to_plot[f'SHMS {y}_vs_{x}'] = h

#%%
# create simc histos
simc_hist_to_plot = {}
for y in simc_var_y:
    hy = simc_var_y[y][0]
    ylim = simc_var_y[y][1]
    for x in simc_var_hx:
        hx = simc_var_hx[x][0]
        xlim = simc_var_hx[x][1]
        h = {}
        
        if x[0] == y[0] or y == 'Emiss' or y == 'W':
            for s in hx:
                    h[s] = B.histo2d(hx[s],
                                     hy[s],
                                     range=[xlim,ylim],bins=100,logz=True,
                                     xlabel=x,ylabel=y,
                                     title=f'{y}_vs_{x} setting {s}')
                    
            simc_hist_to_plot[f'HMS {y}_vs_{x}'] = h           
        
    for x in simc_var_px:
        hx = simc_var_px[x][0]
        xlim = simc_var_px[x][1]
        h = {}
        
        if x[0] == y[0] or y == 'Emiss' or y == 'W':
            for s in hx:
                    h[s] = B.histo2d(hx[s],
                                     hy[s],
                                     range=[xlim,ylim],bins=100,logz=True,
                                     xlabel=x,ylabel=y,
                                     title=f'{y}_vs_{x} setting {s}')
            
            simc_hist_to_plot[f'SHMS {y}_vs_{x}'] = h
#%% make histograms single run
#HMS histos
# W with target variables
h_xptar_W_histos =\
    B.histo2d(hxptar,W,range=[(-0.08,0.08),(0.5,1.5)],bins=100,logz=True,
              xlabel='xptar',ylabel='W')
h_yptar_W_histos =\
    B.histo2d(hyptar,W,range=[(-0.03,0.03),(0.5,1.5)],bins=100,logz=True,
              xlabel='yptar',ylabel='W')
# Emiss with target variables
h_xptar_Em_histos =\
    B.histo2d(hxptar,Emiss,range=[(-0.08,0.08),(-0.1,0.1)],bins=100,logz=True,
              xlabel='xptar',ylabel='Em')
h_yptar_Em_histos =\
    B.histo2d(hyptar,Emiss,range=[(-0.03,0.03),(-0.1,0.1)],bins=100,logz=True,
              xlabel='yptar',ylabel='Em')
# hdelta with target variables
h_xptar_hdp_histos =\
    B.histo2d(hxptar,hdp,range=[(-0.08,0.08),(-10.,10.)],bins=100,logz=True,
              xlabel='xptar',ylabel='hdp')
h_yptar_hdp_histos =\
    B.histo2d(hyptar,hdp,range=[(-0.03,0.03),(-10.,10.)],bins=100,logz=True,
              xlabel='yptar',ylabel='hdp')    
# W with focal plane variables
h_yfp_W_histos =\
    B.histo2d(hyfp,W,range=[(-10.,20.),(0.5,1.5)],bins=100,logz=True,
              xlabel='yfp',ylabel='W')
h_xfp_W_histos =\
    B.histo2d(hxfp,W,range=[(-30.,30.),(0.5,1.5)],bins=100,logz=True,
              xlabel='xfp',ylabel='W')
h_ypfp_W_histos =\
    B.histo2d(hypfp,W,range=[(-0.02,0.02),(0.5,1.5)],bins=100,logz=True,
              xlabel='ypfp',ylabel='W')
h_xpfp_W_histos =\
        B.histo2d(hxpfp,W,range=[(-0.04,0.04),(0.5,1.5)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='W')
# Em with focal plane variables
h_yfp_Em_histos =\
    B.histo2d(hyfp,Emiss,range=[(-10.,20.),(-0.1,1.)],bins=100,logz=True,
              xlabel='yfp',ylabel='Em')
h_xfp_Em_histos =\
    B.histo2d(hxfp,Emiss,range=[(-30.,30.),(-0.1,1.)],bins=100,logz=True,
              xlabel='xfp',ylabel='Em')
h_ypfp_Em_histos =\
    B.histo2d(hypfp,Emiss,range=[(-0.02,0.02),(-0.1,1.)],bins=100,logz=True,
              xlabel='ypfp',ylabel='Em')
h_xpfp_Em_histos =\
        B.histo2d(hxpfp,Emiss,range=[(-0.04,0.04),(-0.1,1.)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='Em')
# hdp with focal plane variables
h_yfp_hdp_histos =\
    B.histo2d(hyfp,hdp,range=[(-10.,20.),(-10.,10.)],bins=100,logz=True,
              xlabel='yfp',ylabel='hdp')
h_xfp_hdp_histos =\
    B.histo2d(hxfp,hdp,range=[(-30.,30.),(-10.,10.)],bins=100,logz=True,
              xlabel='xfp',ylabel='hdp')
h_ypfp_hdp_histos =\
    B.histo2d(hypfp,hdp,range=[(-0.02,0.02),(-10.,10.)],bins=100,logz=True,
              xlabel='ypfp',ylabel='hdp')
h_xpfp_hdp_histos =\
        B.histo2d(hxpfp,hdp,range=[(-0.04,0.04),(-10.,10.)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='hdp')
# hms xfp_yfp xpfp_ypfp
h_xfp_yfp_histos =\
    B.histo2d(hxfp,hyfp,range=[(-30.,30.),(-10.,20.)],bins=100,logz=True,
              xlabel='xfp',ylabel='yfp')
h_xpfp_ypfp_histos =\
        B.histo2d(hxpfp,hypfp,range=[(-0.04,0.04),(-0.02,0.02)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='ypfp')          
    
#HMS histos SIMC
# W with target variables
sh_xptar_W_histos =\
    B.histo2d(shxptar,sW,range=[(-0.08,0.08),(0.5,1.5)],bins=100,logz=True,
              xlabel='xptar',ylabel='W')
sh_yptar_W_histos =\
    B.histo2d(shyptar,sW,range=[(-0.03,0.03),(0.5,1.5)],bins=100,logz=True,
              xlabel='yptar',ylabel='W')
# Emiss with target variables
sh_xptar_Em_histos =\
    B.histo2d(shxptar,sEmiss,range=[(-0.08,0.08),(-0.1,0.1)],bins=100,logz=True,
              xlabel='xptar',ylabel='Em')
sh_yptar_Em_histos =\
    B.histo2d(shyptar,sEmiss,range=[(-0.03,0.03),(-0.1,0.1)],bins=100,logz=True,
              xlabel='yptar',ylabel='Em')
# hdelta with target variables
sh_xptar_hdp_histos =\
    B.histo2d(shxptar,shdp,range=[(-0.08,0.08),(-10.,10.)],bins=100,logz=True,
              xlabel='xptar',ylabel='hdp')
sh_yptar_hdp_histos =\
    B.histo2d(shyptar,shdp,range=[(-0.03,0.03),(-10.,10.)],bins=100,logz=True,
              xlabel='yptar',ylabel='hdp')   
# W with focal plane variables
sh_yfp_W_histos =\
    B.histo2d(shyfp,sW,range=[(-10.,20.),(0.5,1.5)],bins=100,logz=True,
              xlabel='yfp',ylabel='W')
sh_xfp_W_histos =\
    B.histo2d(shxfp,sW,range=[(-30.,30.),(0.5,1.5)],bins=100,logz=True,
              xlabel='xfp',ylabel='W')
sh_ypfp_W_histos =\
    B.histo2d(shypfp,sW,range=[(-0.02,0.02),(0.5,1.5)],bins=100,logz=True,
              xlabel='ypfp',ylabel='W')
sh_xpfp_W_histos =\
        B.histo2d(shxpfp,sW,range=[(-0.04,0.04),(0.5,1.5)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='W')
# Em with focal plane variables
sh_yfp_Em_histos =\
    B.histo2d(shyfp,sEmiss,range=[(-10.,20.),(-0.1,1.)],bins=100,logz=True,
              xlabel='yfp',ylabel='Em')
sh_xfp_Em_histos =\
    B.histo2d(shxfp,sEmiss,range=[(-30.,30.),(-0.1,1.)],bins=100,logz=True,
              xlabel='xfp',ylabel='Em')
sh_ypfp_Em_histos =\
    B.histo2d(shypfp,sEmiss,range=[(-0.02,0.02),(-0.1,1.)],bins=100,logz=True,
              xlabel='ypfp',ylabel='Em')
sh_xpfp_Em_histos =\
        B.histo2d(shxpfp,sEmiss,range=[(-0.04,0.04),(-0.1,1.)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='Em')
# hdp with focal plane variables
sh_yfp_hdp_histos =\
    B.histo2d(shyfp,shdp,range=[(-10.,20.),(-10.,10.)],bins=100,logz=True,
              xlabel='yfp',ylabel='hdp')
sh_xfp_hdp_histos =\
    B.histo2d(shxfp,shdp,range=[(-30.,30.),(-10.,10.)],bins=100,logz=True,
              xlabel='xfp',ylabel='hdp')
sh_ypfp_hdp_histos =\
    B.histo2d(shypfp,shdp,range=[(-0.02,0.02),(-10.,10.)],bins=100,logz=True,
              xlabel='ypfp',ylabel='hdp')
sh_xpfp_hdp_histos =\
            B.histo2d(shxpfp,shdp,range=[(-0.04,0.04),(-10.,10.)],bins=100,logz=True,
                      xlabel='xpfp',ylabel='hdp')
# hms xfp_yfp xpfp_ypfp
sh_xfp_yfp_histos =\
    B.histo2d(shxfp,shyfp,range=[(-30.,30.),(-10.,20.)],bins=100,logz=True,
              xlabel='xfp',ylabel='yfp')
sh_xpfp_ypfp_histos =\
        B.histo2d(shxpfp,shypfp,range=[(-0.04,0.04),(-0.02,0.02)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='ypfp') 
            
#SHMS histos
# W with target variables
p_xptar_W_histos =\
    B.histo2d(pxptar,W,range=[(-0.03,0.03),(0.5,1.5)],bins=100,logz=True,
              xlabel='xptar',ylabel='W')
p_yptar_W_histos =\
    B.histo2d(pyptar,W,range=[(-0.015,0.015),(0.5,1.5)],bins=100,logz=True,
              xlabel='yptar',ylabel='W')
# Emiss with target variables
p_xptar_Em_histos =\
    B.histo2d(pxptar,Emiss,range=[(-0.03,0.03),(-0.1,0.1)],bins=100,logz=True,
              xlabel='xptar',ylabel='Em')
p_yptar_Em_histos =\
    B.histo2d(pyptar,Emiss,range=[(-0.015,0.015),(-0.1,0.1)],bins=100,logz=True,
              xlabel='yptar',ylabel='Em')
# pdelta with target variables
p_xptar_pdp_histos =\
    B.histo2d(pxptar,pdp,range=[(-0.03,0.03),(-10.,5.)],bins=100,logz=True,
              xlabel='xptar',ylabel='pdp')
p_yptar_pdp_histos =\
    B.histo2d(pyptar,pdp,range=[(-0.015,0.015),(-10.,5.)],bins=100,logz=True,
              xlabel='yptar',ylabel='pdp')    
# W with focal plane variables
p_yfp_W_histos =\
    B.histo2d(pyfp,W,range=[(-10.,5.),(0.5,1.5)],bins=100,logz=True,
              xlabel='yfp',ylabel='W')
p_xfp_W_histos =\
    B.histo2d(pxfp,W,range=[(-10.,10.),(0.5,1.5)],bins=100,logz=True,
              xlabel='xfp',ylabel='W')
p_ypfp_W_histos =\
    B.histo2d(pypfp,W,range=[(-0.02,0.02),(0.5,1.5)],bins=100,logz=True,
              xlabel='ypfp',ylabel='W')
p_xpfp_W_histos =\
        B.histo2d(pxpfp,W,range=[(-0.04,0.04),(0.5,1.5)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='W')
# Em with focal plane variables
p_yfp_Em_histos =\
    B.histo2d(pyfp,Emiss,range=[(-10.,10.),(-0.1,0.1)],bins=100,logz=True,
              xlabel='yfp',ylabel='Em')
p_xfp_Em_histos =\
    B.histo2d(pxfp,Emiss,range=[(-10.,10.),(-0.1,0.1)],bins=100,logz=True,
              xlabel='xfp',ylabel='Em')
p_ypfp_Em_histos =\
    B.histo2d(pypfp,Emiss,range=[(-0.02,0.02),(-0.1,0.1)],bins=100,logz=True,
              xlabel='ypfp',ylabel='Em')
p_xpfp_Em_histos =\
        B.histo2d(pxpfp,Emiss,range=[(-0.04,0.04),(-0.1,0.1)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='Em')
# pdp with focal plane variables
p_yfp_pdp_histos =\
    B.histo2d(pyfp,pdp,range=[(-10.,20.),(-10.,5.)],bins=100,logz=True,
              xlabel='yfp',ylabel='pdp')
p_xfp_pdp_histos =\
    B.histo2d(pxfp,pdp,range=[(-20.,10.),(-10.,5.)],bins=100,logz=True,
              xlabel='xfp',ylabel='pdp')
p_ypfp_pdp_histos =\
    B.histo2d(pypfp,pdp,range=[(-0.02,0.02),(-10.,5.)],bins=100,logz=True,
              xlabel='ypfp',ylabel='pdp')
p_xpfp_pdp_histos =\
        B.histo2d(pxpfp,pdp,range=[(-0.04,0.04),(-10.,5.)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='pdp')        
# shms xfp_yfp xpfp_ypfp
p_xfp_yfp_histos =\
    B.histo2d(pxfp,pyfp,range=[(-30.,30.),(-10.,20.)],bins=100,logz=True,
              xlabel='xfp',ylabel='yfp')
p_xpfp_ypfp_histos =\
        B.histo2d(pxpfp,pypfp,range=[(-0.04,0.04),(-0.02,0.02)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='ypfp') 
    
#SHMS histos SIMC
# W with target variables
sp_xptar_W_histos =\
    B.histo2d(spxptar,sW,range=[(-0.03,0.03),(0.5,1.5)],bins=100,logz=True,
              xlabel='xptar',ylabel='W')
sp_yptar_W_histos =\
    B.histo2d(spyptar,sW,range=[(-0.015,0.015),(0.5,1.5)],bins=100,logz=True,
              xlabel='yptar',ylabel='W')
# Emiss with target variables
sp_ytar_Em_histos =\
    B.histo2d(spytar,sEmiss,range=[(-2.,2.),(-0.1,0.1)],bins=100,logz=True,
              xlabel='ytar',ylabel='Em')
sp_xptar_Em_histos =\
    B.histo2d(spxptar,sEmiss,range=[(-0.03,0.03),(-0.1,0.1)],bins=100,logz=True,
              xlabel='xptar',ylabel='Em')
sp_yptar_Em_histos =\
    B.histo2d(spyptar,sEmiss,range=[(-0.015,0.015),(-0.1,0.1)],bins=100,logz=True,
              xlabel='yptar',ylabel='Em')
# pdelta with target variables
sp_xptar_pdp_histos =\
    B.histo2d(spxptar,spdp,range=[(-0.03,0.03),(-10,5.)],bins=100,logz=True,
              xlabel='xptar',ylabel='pdp')
sp_yptar_pdp_histos =\
    B.histo2d(spyptar,spdp,range=[(-0.015,0.015),(-10,5.)],bins=100,logz=True,
              xlabel='yptar',ylabel='pdp') 
# W with focal plane variables
sp_yfp_W_histos =\
    B.histo2d(spyfp,sW,range=[(-10.,5.),(0.5,1.5)],bins=100,logz=True,
              xlabel='yfp',ylabel='W')
sp_xfp_W_histos =\
    B.histo2d(spxfp,sW,range=[(-10.,10.),(0.5,1.5)],bins=100,logz=True,
              xlabel='xfp',ylabel='W')
sp_ypfp_W_histos =\
    B.histo2d(spypfp,sW,range=[(-0.02,0.02),(0.5,1.5)],bins=100,logz=True,
              xlabel='ypfp',ylabel='W')
sp_xpfp_W_histos =\
        B.histo2d(spxpfp,sW,range=[(-0.04,0.04),(0.5,1.5)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='W')
# Em with focal plane variables
sp_yfp_Em_histos =\
    B.histo2d(spyfp,sEmiss,range=[(-10.,10.),(-0.1,0.1)],bins=100,logz=True,
              xlabel='yfp',ylabel='Em')
sp_xfp_Em_histos =\
    B.histo2d(spxfp,sEmiss,range=[(-10.,10.),(-0.1,0.1)],bins=100,logz=True,
              xlabel='xfp',ylabel='Em')
sp_ypfp_Em_histos =\
    B.histo2d(spypfp,sEmiss,range=[(-0.02,0.02),(-0.1,0.1)],bins=100,logz=True,
              xlabel='ypfp',ylabel='Em')
sp_xpfp_Em_histos =\
        B.histo2d(spxpfp,sEmiss,range=[(-0.04,0.04),(-0.1,0.1)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='Em')
# pdp with focal plane variables
sp_yfp_pdp_histos =\
    B.histo2d(spyfp,spdp,range=[(-10.,20.),(-10,5.)],bins=100,logz=True,
              xlabel='yfp',ylabel='pdp')
sp_xfp_pdp_histos =\
    B.histo2d(spxfp,spdp,range=[(-20.,10.),(-10,5.)],bins=100,logz=True,
              xlabel='xfp',ylabel='pdp')
sp_ypfp_pdp_histos =\
    B.histo2d(spypfp,spdp,range=[(-0.02,0.02),(-10,5.)],bins=100,logz=True,
              xlabel='ypfp',ylabel='pdp')
sp_xpfp_pdp_histos =\
            B.histo2d(spxpfp,spdp,range=[(-0.04,0.04),(-10.,5.)],bins=100,logz=True,
                      xlabel='xpfp',ylabel='pdp')             
# shms xfp_yfp xpfp_ypfp
sp_xfp_yfp_histos =\
    B.histo2d(spxfp,spyfp,range=[(-30.,30.),(-10.,20.)],bins=100,logz=True,
              xlabel='xfp',ylabel='yfp')
sp_xpfp_ypfp_histos =\
        B.histo2d(spxpfp,spypfp,range=[(-0.04,0.04),(-0.02,0.02)],bins=100,logz=True,
                  xlabel='xpfp',ylabel='ypfp') 

#HMS figure groups
h_kins_tar = [h_xptar_W_histos,h_yptar_W_histos,
                  h_xptar_Em_histos,h_yptar_Em_histos]

sh_kins_tar = [sh_xptar_W_histos,sh_yptar_W_histos,
                   sh_xptar_Em_histos,sh_yptar_Em_histos]

h_kins_fp = [h_yfp_W_histos,h_xfp_W_histos,h_ypfp_W_histos,h_xpfp_W_histos,
               h_yfp_Em_histos,h_xfp_Em_histos,h_ypfp_Em_histos,h_xpfp_Em_histos]

sh_kins_fp = [sh_yfp_W_histos,sh_xfp_W_histos,sh_ypfp_W_histos,sh_xpfp_W_histos,
               sh_yfp_Em_histos,sh_xfp_Em_histos,sh_ypfp_Em_histos,sh_xpfp_Em_histos]

#SHMS figure groups
p_kins_tar = [p_xptar_W_histos,p_yptar_W_histos,
                  p_xptar_Em_histos,p_yptar_Em_histos]

sp_kins_tar = [sp_xptar_W_histos,sp_yptar_W_histos,sp_xptar_Em_histos,sp_yptar_Em_histos]

p_kins_fp = [p_yfp_W_histos,p_xfp_W_histos,p_ypfp_W_histos,p_xpfp_W_histos,
               p_yfp_Em_histos,p_xfp_Em_histos,p_ypfp_Em_histos,p_xpfp_Em_histos,]

sp_kins_fp = [sp_yfp_W_histos,sp_xfp_W_histos,sp_ypfp_W_histos,sp_xpfp_W_histos,
               sp_yfp_Em_histos,sp_xfp_Em_histos,sp_ypfp_Em_histos,sp_xpfp_Em_histos]

#Delta plots
dp_tar = [h_xptar_hdp_histos,h_yptar_hdp_histos,
          p_xptar_pdp_histos,p_yptar_pdp_histos]

sdp_tar = [sh_xptar_hdp_histos,sh_yptar_hdp_histos,
          sp_xptar_pdp_histos,sp_yptar_pdp_histos]

dp_fp = [h_yfp_hdp_histos,h_xfp_hdp_histos,h_ypfp_hdp_histos,h_xpfp_hdp_histos,
         h_xfp_yfp_histos,h_xpfp_ypfp_histos,
         p_yfp_pdp_histos,p_xfp_pdp_histos,p_ypfp_pdp_histos,p_xpfp_pdp_histos,
         p_xfp_yfp_histos,p_xpfp_ypfp_histos]

sdp_fp = [sh_yfp_hdp_histos,sh_xfp_hdp_histos,sh_ypfp_hdp_histos,sh_xpfp_hdp_histos,
         sh_xfp_yfp_histos,sh_xpfp_ypfp_histos,
         sp_yfp_pdp_histos,sp_xfp_pdp_histos,sp_ypfp_pdp_histos,sp_xpfp_pdp_histos,
         sp_xfp_yfp_histos,sp_xpfp_ypfp_histos]

#%% make DATA-SIMC 2D histograms

diff_histos = {}
for h in data_hist_to_plot:
    dh = {}
    for r in data.many:       
        s = data.setting[r]
        histo_name = f'{h} $\Delta$ \nrun {r} setting {s}'
        histo = data_hist_to_plot[h][r]/simc_hist_to_plot[h][s]
        histo.title = histo_name
        
        dh[r] = histo     
    diff_histos[histo_name] = dh
        
#%% plot kins and delta vs xtar/ytar and fp var per run/setting and spectrometer
# first select a run and spectrometer arm ('H' -> HMS or 'S' -> SHMS)
run = 20840 # heep_coin runs [20840,20841,20846,20851,20858,20861,20868,20869]
setting = data.setting[run]
spec = 'SHMS'

h_list = [data_hist_to_plot[plot][run] for plot in data_hist_to_plot 
          if plot[0]==spec[0]]

rows = 2
cols = 4

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'Run {run}', fontsize=16)
i = 1
j = 1   
flag = False 
#create directory
path = f'./optics/plots/{run}-newfit_cafe4_deltaoptim2/{spec}/'
os.makedirs(path, exist_ok=True)
for histo in h_list:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(h.xlabel)
    ax.set_ylabel(h.ylabel)
    ax.set_title('')
    extent = full_extent(ax,0.20,0.,45.,0.).transformed(fig.dpi_scale_trans.inverted())
    # B.pl.savefig(path + f'{h.title}.png',bbox_inches=extent)
    if flag:
        break
    if i > 8:
        B.pl.savefig(path + f'{run} {spec} summary {j}.png')
        if j < 3:
            fig = B.pl.figure(layout='constrained',figsize=(20,10))
            fig.suptitle(f'Run {run}', fontsize=16)
            i = 1
            j+=1
   
sh_list = [simc_hist_to_plot[plot][setting] for plot in simc_hist_to_plot 
           if plot[0]==spec[0]]
       
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'SIMC setting {setting}', fontsize=16)
i = 1
j = 1
path = f'./optics/plots/{run}-newfit_cafe4_deltaoptim2/{spec}/SIMC/'
os.makedirs(path, exist_ok=True)
for histo in sh_list:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')
    if i > 6: 
        B.pl.savefig(path + f'{setting} {spec} summary {j}.png')
        if j < 3:
            fig = B.pl.figure(layout='constrained',figsize=(20,10))
            fig.suptitle(f'SIMC setting {setting}', fontsize=16)
            i = 1
            j+=1
        
#%% plot a desired kin (y) for a chosen spectrometer ('HMS'/'SHMS') 
# for all runs/settings 

#first pick a spectrometer and and kin you want to look at     
spec = 'SHMS'  # 5 for SHMS, 4 for HMS
y = 'W' # options are 'W', 'Emiss', and 'hdp'/'pdp'

h_list = [data_hist_to_plot[plot] for plot in data_hist_to_plot 
          if plot[5]==y]
sh_list = [simc_hist_to_plot[plot] for plot in simc_hist_to_plot 
           if plot[5]==y]
# diff_list = [diff_histos[plot] for plot in diff_histos 
#              if plot[5]==y]

#%% plot data 
rows = 2
cols = 4

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'{spec}', fontsize=16)
i = 1
j = 1

for histo in h_list:
    for r in data.many:
        ax = plt.subplot(rows,cols,i)
        h = histo[r]
        h.plot(colormap=cmp['viridis'],logz=True)
        i+=1
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(h.title)
        if i > (rows*cols) and j < 8:
            path = f'./optics/plots/{r}/{spec}/{y}/'
            #os.makedirs(path, exist_ok=True)
            # B.pl.savefig(path + f'{j}') 
            fig = B.pl.figure(layout='constrained',figsize=(20,10))
            fig.suptitle(f'{spec}', fontsize=16)
            i = 1
            j+=1
           
        
#%% plot simc   
rows = 2
cols = 3

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'SIMC {spec}', fontsize=16)
i = 1
j = 1
for histo in sh_list:
    for s in simc.many:
        ax = plt.subplot(rows,cols,i)
        h = histo[s]
        h.plot(colormap=cmp['viridis'])
        i+=1
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(h.title)
        if i > 6 and j < 6:
            fig = B.pl.figure(layout='constrained',figsize=(20,10))
            fig.suptitle(f'SIMC {spec}', fontsize=16)
            i = 1
            j+=1

#%% plot diffs
rows = 2
cols = 4

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'{spec}', fontsize=16)
i = 1
j = 1

for histo in diff_list:
    for r in data.many:
        ax = plt.subplot(rows,cols,i)
        h = histo[r]
        h.plot(colormap=cmp['viridis'])
        i+=1
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(h.title)
        if i > (rows*cols) and j < 6:
            break
            path = f'./optics/plots/{r}/{spec}/{y}/'
            #os.makedirs(path, exist_ok=True)
            # B.pl.savefig(path + f'{j}') 
            fig = B.pl.figure(layout='constrained',figsize=(20,10))
            fig.suptitle(f'{spec}', fontsize=16)
            i = 1
            j+=1   


#%% HMS kins vs tar 
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'HMS Run {runs}\nKinematic Variables vs. Target Variables', fontsize=16)
rows = 2
cols = 3
i = 1
for histo in h_kins_tar:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('SIMC HMS\nKinematic Variables vs. Target Variables', fontsize=16)
rows = 2
cols = 3
i = 1
for histo in sh_kins_tar:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

#%% HMS kins vs fp
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'HMS Run {runs}\nKinematic Variables vs. Focal Plane Variables', fontsize=16)
rows = 2
cols = 4
i = 1
for histo in h_kins_fp:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('HMS SIMC\nKinematic Variables vs. Focal Plane Variables', fontsize=16)
rows = 2
cols = 4
i = 1
for histo in sh_kins_fp:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

#%% SHMS kins vs tar 
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'SHMS Run {runs}\nKinematic Variables vs. Target Variables', fontsize=16)
rows = 2
cols = 3
i = 1
for histo in p_kins_tar:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('SIMC SHMS\nKinematic Variables vs. Target Variables', fontsize=16)
rows = 2
cols = 3
i = 1
for histo in sp_kins_tar:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

#%% SHMS kins vs fp
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'SHMS Run {runs}\nKinematic Variables vs. Focal Plane Variables', fontsize=16)
rows = 2
cols = 4
i = 1
for histo in p_kins_fp:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('SHMS SIMC\nKinematic Variables vs. Focal Plane Variables', fontsize=16)
rows = 2
cols = 4
i = 1
for histo in sp_kins_fp:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

#%% delta vs tar
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'Run {runs}\n$\delta$ vs. Target Variables', fontsize=16)
rows = 2
cols = 3
i = 1
for histo in dp_tar:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('SIMC\n$\delta$ vs. Target Variables', fontsize=16)
rows = 2
cols = 3
i = 1
for histo in sdp_tar:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

#%% delta vs fp
fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle(f'Run {runs}\n$\delta$ vs. Focal Plane Variables', fontsize=16)
rows = 3
cols = 4
i = 1
for histo in dp_fp:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('SIMC\n$\delta$ vs. Focal Plane Variables', fontsize=16)
rows = 3
cols = 4
i = 1
for histo in sdp_fp:
    ax = plt.subplot(rows,cols,i)
    h = histo
    h.plot(colormap=cmp['viridis'])
    i+=1
    ax.set_xlabel(histo.xlabel)
    ax.set_ylabel(histo.ylabel)
    ax.set_title('')              