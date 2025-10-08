#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 11:40:52 2025

@author: gvill

HMS/SHMS Collimator cut 

"""
import LT.box as B
import data_init as D
import cut_handler as C
import numpy as np
from matplotlib import colormaps as cmp

#%% functions live here

def is_inside(n, x, y, nlim, xlim, ylim):
    is_point_inside = []
    for i in range(n):
        ncross = 0
        for j in range(nlim-1):
            if (ylim[j] < y[i]) and (ylim[j+1] < y[i]): continue
            if (x[i] == xlim[j]): continue        
            t = x[i] - xlim[j]
            s = xlim[j+1] - x[i]
            if (t*s < 0.): continue       
            di = (ylim[j+1] - ylim[j])/(xlim[j+1] - xlim[j])
            f = ylim[j] + di*(x[i] - xlim[j])
            if (f < y[i]): continue       
            ncross += 1   
        is_point_inside.append((ncross % 2) == 1)    
    return np.array(is_point_inside)

#%% Load hcana variables

br_sel = ['H.extcor.xsieve','H.extcor.ysieve',
          'P.extcor.xsieve','P.extcor.ysieve']

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','theta_e','theta_p',
               'e_pf','h_pf','e_xptar','e_yptar','h_xptar','h_yptar']
#%% load data root files
# make sure files are where they should be!
# RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
#                   select_branches={'T':br_sel})
RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                  setting='delta_scan_0',select_branches={'T':br_sel})
#%%
# load simc root files
# make sure files are where they should be!
SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                   select_branches={'SNT':br_sel_SIMC},simc_type='-')

#%%

x_coll = RUN.Branches['H.extcor.xsieve']
y_coll = RUN.Branches['H.extcor.ysieve']

px_coll = RUN.Branches['P.extcor.xsieve']
py_coll = RUN.Branches['P.extcor.ysieve']

h_xcoll_ycoll = B.histo2d(y_coll,x_coll,range=[(-15,15),(-15,15)],bins=100)

ph_xcoll_ycoll = B.histo2d(py_coll,px_coll,range=[(-15,15),(-15,15)],bins=100)

coll_lim = np.array([(-2.2875,11.646),(2.2875,11.646),(4.575,5.823),
                        (4.575,-5.823),(2.2875,-11.646),(-2.2875,-11.646),
                        (-4.575,-5.823),(-4.575,5.823),(-2.2875,11.646)])

pcoll_lim = np.array([(-4.25,12.5),(4.25,12.5),(8.5,6.25),(8.5,-6.25),
                      (4.25,-12.5),(-4.25,-12.5),(-8.5,-6.25),(-8.5,6.25),
                      (-4.25,12.5)])

h_xcoll_ycoll.plot(colormap=cmp['viridis'])
B.pl.plot(coll_lim[:,0],coll_lim[:,1],color='r')
B.pl.figure()
ph_xcoll_ycoll.plot(colormap=cmp['viridis'])
B.pl.plot(pcoll_lim[:,0],pcoll_lim[:,1],color='r')

hex_cut = h_xcoll_ycoll.poly_cut([(-2.2875,11.646),(2.2875,11.646),(4.575,5.823),
                        (4.575,-5.823),(2.2875,-11.646),(-2.2875,-11.646),
                        (-4.575,-5.823),(-4.575,5.823)])
phex_cut = ph_xcoll_ycoll.poly_cut([(-4.25,12.5),(4.25,12.5),(8.5,6.25),(8.5,-6.25),
                      (4.25,-12.5),(-4.25,-12.5),(-8.5,-6.25),(-8.5,6.25),
                      (-4.25,12.5)])


h_cut = h_xcoll_ycoll
h_cut.bin_content = np.where(hex_cut,h_cut.bin_content,0)
h_cut.bin_error = np.where(hex_cut,h_cut.bin_error,0)

ph_cut = ph_xcoll_ycoll
bothcut = phex_cut & hex_cut
ph_cut.bin_content = np.where(bothcut,ph_cut.bin_content,0)
ph_cut.bin_error = np.where(bothcut,ph_cut.bin_error,0)

B.pl.figure()
h_cut.plot(colormap=cmp['viridis'])
B.pl.plot(coll_lim[:,0],coll_lim[:,1],color='r')
B.pl.figure()
ph_cut.plot(colormap=cmp['viridis'])
B.pl.plot(pcoll_lim[:,0],pcoll_lim[:,1],color='r')
#%%
no_nonsense = (x_coll<20.) & (x_coll>-20.) & (y_coll<20.) & (y_coll>-20)

xcoll = x_coll[no_nonsense]
ycoll = y_coll[no_nonsense]

B.pl.figure()
B.pl.plot(ycoll,xcoll,'.')
B.pl.plot(coll_lim[:,0],coll_lim[:,1])

hres = is_inside(x_coll.size,y_coll,x_coll,coll_lim[:,0].size,
                coll_lim[:,0],coll_lim[:,1])
# B.pl.plot(y_coll[res],x_coll[res],'.')
#%%
no_nonsense = (px_coll<20.) & (px_coll>-20.) & (py_coll<20.) & (py_coll>-20)

pxcoll = px_coll[no_nonsense]
pycoll = py_coll[no_nonsense]

B.pl.figure()
B.pl.plot(pycoll,pxcoll,'.')
B.pl.plot(pcoll_lim[:,0],pcoll_lim[:,1])

res = is_inside(px_coll.size,py_coll,px_coll,pcoll_lim[:,0].size,
                pcoll_lim[:,0],pcoll_lim[:,1])

r = hres & res
B.pl.plot(py_coll[res],px_coll[res],'.')