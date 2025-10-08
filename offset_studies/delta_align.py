#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 10:05:36 2025

@author: gvill

This code calculates (up to some order), the 0th matrix element
from the difference true_delta - delta 

"""

import numpy as np
import matplotlib.lines as mlines
from matplotlib import colormaps as cmp
import LT.box as B
import cut_handler as C
import data_init as D

#%% helper functions
def fit_legend(par):
    ax = B.pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    for i in range(len(par)):
        fit_info = (
               f"par {i} = {par[i].value:.5e} $\pm$ {par[i].err:.5e}"
           )
        # if i == len(par)-1:
        #     fit_info = (
        #            f"par {i} = {par[i].value:.5f} $\pm$ {par[i].err:.5f}"
        #        )
        # else:    
        #     fit_info = (
        #            f"par {i} = {par[i].value:.5f} $\pm$ {par[i].err:.5f}\n"
        #        )
        handles.append(mlines.Line2D([], [], color='None', label=fit_info))
        labels.append(fit_info)
    
    ax.legend(handles=handles, labels=labels, loc='lower right')
    return

def combine_histos(histos,many=[]):
    loop1 = True 
    for m in many:
        if loop1:
            hsum = histos[m]
            loop1 = False
        else:
            hsum += histos[m]
    
    return hsum

#%% Constants
# All masses/energy/momenta in GeV/c2

Mp_rad = 0.944  #mass of the proton in SIMC (using radcorr?)
Mp = 0.93827
Me = 0.000511
Eb = 10.542     #beam energy
Pc = 8.55       #SHMS central momentum setting (electron side)

#%% select hcana variables

kin_var = ['P.kin.primary.scat_ang_rad','P.kin.primary.x_bj','P.dc.x_fp']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

br_sel = kin_var + acc_var + calPID_var

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','theta_e','theta_p',
               'e_pf','h_pf','e_xptar','e_yptar','h_xptar','h_yptar']

#%% load data root files
# make sure files are where they should be!
RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                  select_branches={'T':br_sel})

# RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin',
#                   setting = 'delta_scan_0', select_branches={'T':br_sel})
#%% assign needed variables

# the = RUN.Branches['P.kin.primary.scat_ang_rad']
# xbj = RUN.Branches['P.kin.primary.x_bj']
# x_fp = RUN.Branches['P.dc.x_fp']
# delta = RUN.Branches['P.gtr.dp']

the = {}
xbj = {}
x_fp = {}
delta = {}
for m in RUN.many:
    the[m] = RUN.Branches[m]['P.kin.primary.scat_ang_rad']
    xbj[m] = RUN.Branches[m]['P.kin.primary.x_bj']
    x_fp[m] = RUN.Branches[m]['P.dc.x_fp']
    delta[m] = RUN.Branches[m]['P.gtr.dp']

#%% calculate true quantities

# nom = Mp_rad*Mp_rad - Mp*Mp - 2*Mp*Eb
# denom = -2*Mp - 4*Eb*np.sin(the*0.5)*np.sin(the*0.5)
# Eptrue = nom/denom

# Ptrue = np.sqrt(Eptrue*Eptrue - Me*Me)

# delta_true = 100*((Ptrue/Pc) - 1)   # [%]

# delta_diff = delta_true - delta

delta_diff = {}
for m in RUN.many:
    nom = Mp_rad*Mp_rad - Mp*Mp - 2*Mp*Eb
    denom = -2*Mp - 4*Eb*np.sin(the[m]*0.5)*np.sin(the[m]*0.5)
    Eptrue = nom/denom
    
    Ptrue = np.sqrt(Eptrue*Eptrue - Me*Me)
    
    delta_true = 100*((Ptrue/Pc) - 1)   # [%]
    
    delta_diff[m] = delta_true - delta[m]

#%% define cuts and apply them

# xbj_cut = C.WCUT(-1.025,1.025,name='xbj_cut')

# cuts_list = C.acceptance_cuts + [xbj_cut]

# cuts_to_apply = []
# for cut in cuts_list:
#     cvar = RUN.Branches[C.HCANA_names[cut.name]]
#     cuts_to_apply.append(cut(cvar))
        
# all_cuts = cuts_to_apply[0]
# for cut in cuts_to_apply:
#     all_cuts = all_cuts & cut

# delta_diff_cut = delta_diff[all_cuts]
# x_fp_cut = x_fp[all_cuts]

xbj_cut = C.WCUT(-1.025,1.025,name='xbj_cut')

cuts_list = C.acceptance_cuts + [xbj_cut]

cuts_to_apply = {}
for m in RUN.many:
    clist = []
    for cut in cuts_list:
        cvar = RUN.Branches[m][C.HCANA_names[cut.name]]
        clist.append(cut(cvar))
    cuts_to_apply[m] = clist    

all_cuts = {}
for m in RUN.many:        
    all_cuts[m] = cuts_to_apply[m][0]
    for cut in cuts_to_apply[m]:
        all_cuts[m] = all_cuts[m] & cut

delta_diff_cut = {}
x_fp_cut = {}
for m in RUN.many:
    delta_diff_cut[m] = delta_diff[m][all_cuts[m]]
    x_fp_cut[m] = x_fp[m][all_cuts[m]]

#%% plot delta_diff vs xfp

# h = B.histo2d(x_fp_cut, delta_diff_cut,range=[(-10,10),(-5,5)],bins=100)

# range_delta = (delta_diff_cut >= -0.8) & (delta_diff_cut <= 0.2)
# # range_delta = (delta_diff_cut >= -0.7) & (delta_diff_cut <= 0.1)
# range_x = (x_fp_cut >= -5) & (x_fp_cut <= 4.5)
# range_cut = range_delta & range_x

# fit = B.polyfit(x_fp_cut[range_cut], delta_diff_cut[range_cut],order=2)
# fit.plot(color='r',label='polyfit')

# h.plot(colormap = cmp['viridis'])

# B.pl.title('delta_diff = delta_true - delta vs. x_fp_cut')
# B.pl.xlabel('[cm]')
# B.pl.ylabel('[%]')
# fit_legend(fit.parameters)

delta_diff_xfp_h = {}
for m in RUN.many:

    h = B.histo2d(x_fp_cut[m], delta_diff_cut[m],
                     range=[(-20,22),(-5,5)],bins=100)
    # B.pl.figure()
    # h.plot()
    delta_diff_xfp_h[m] = h
    
    if m == 20840:
        xfp_fit = x_fp_cut[m]
        delta_fit = delta_diff_cut[m]
    else:
        xfp_fit = np.append(xfp_fit,x_fp_cut[m])
        delta_fit = np.append(delta_fit,delta_diff_cut[m])

xfp_fit = np.array(xfp_fit)
delta_fit = np.array(delta_fit)    

delta_diff_xfp_h_all = combine_histos(delta_diff_xfp_h,RUN.many)

range_delta = (delta_fit >= -0.8) & (delta_fit <= 0.3)
range_x = (xfp_fit >= -20) & (xfp_fit <= 22)
range_cut = range_delta & range_x

fit = B.polyfit(xfp_fit[range_cut], delta_fit[range_cut],order=2)
fit1 = B.polyfit(xfp_fit[range_cut], delta_fit[range_cut],order=1)
fit.plot(color='r',label='polyfit')
fit1.plot(color='b',label='linefit')

delta_diff_xfp_h_all.plot(colormap = cmp['viridis'],logz=True)

B.pl.title('delta_diff = delta_true - delta vs. x_fp_cut')
B.pl.xlabel('[cm]')
B.pl.ylabel('[%]')
fit_legend(fit.parameters)



