#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:08:55 2025

@author: gvill

Make DC time vs wirenum plots for SHMS
"""
import data_init as D
import database_operations as db
import LT.box as B

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cmp

#%% helper functions live here

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
#%%
runs = [20840,20841,20846,20851,20858,20861,20868,20869] #heep_coin runs
deut_db_name = 'deuteron_db.db'

# declare dc plane names to loop over
dc1_plane_names = ['1u1','1u2','1x1','1x2','1v2','1v1']
dc2_plane_names = ['2v2','2v1','2x2','2x1','2u2','2u1']

dc_plane_names = dc1_plane_names + dc2_plane_names

# select branches names to load : .dc.

# dc_branches = D.get_list(db.retrieve(deut_db_name, 'Branches', 'TTree_Branches',
#                           where='Type = \'hdc\''))

dc_branches = D.get_list(db.retrieve(deut_db_name, 'Branches', 'TTree_Branches',
                          where='Branches like \'P.dc.%.time\'')) +\
                D.get_list(db.retrieve(deut_db_name, 'Branches', 'TTree_Branches',
                          where='Branches like \'P.dc.%.wirenum\''))
#%%
data = D.DATA_INIT(data_type='deut23_data',run=runs,
                   select_branches=dc_branches)

#%%
time = {}
wirenum = {}

for r in runs:
    thisrun_time = {}
    thisrun_wirenum = {}
    for plane in dc_plane_names:
        time_name = f'P.dc.{plane}.time'
        wirenum_name = f'P.dc.{plane}.wirenum'
        
        thisrun_time[plane] = get_array(data.Branches[r][time_name])    
        thisrun_wirenum[plane] = get_array(data.Branches[r][wirenum_name])
    time[r] = thisrun_time
    wirenum[r] = thisrun_wirenum
        
#%% remove zeroes (no nonsense)

for r in runs:
    t = time[r]
    w = wirenum[r]
    for plane in dc_plane_names:
        nozeros = t[plane] > 0.        
        
        t[plane] = t[plane][nozeros]
        w[plane] = w[plane][nozeros]
    time[r] = t 
    wirenum[r] = w
    
#%% make histograms

time_wirenum_histos = {}
for r in runs:
    thishistos = {}
    for plane in dc_plane_names:
        t = time[r][plane]
        w = wirenum[r][plane]
        
        thishistos[plane] = B.histo2d(w,t,range=[(0,100),(-50,350)],
                                               bins=100,logz=True)
    time_wirenum_histos[r] = thishistos
#%% plot histos in a 6x6 figure

for r in runs:    
    fig = B.pl.figure(layout='constrained',figsize=(20,10))
    fig.suptitle(f'SHMS Drift time vs. Wire number\n Run {r}', fontsize=16)
    rows = 2
    cols = 6
    i = 1
    for plane in dc_plane_names:
        ax = plt.subplot(rows,cols,i)
        h = time_wirenum_histos[r][plane]
        h.plot(colormap=cmp['viridis'])
        i+=1
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'{plane}')

    #B.pl.savefig(f'./shms_dc_calib/plots/htime_wirenum_{r}.png')
