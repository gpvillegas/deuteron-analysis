#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 14:20:52 2025

@author: gvill
"""

import LT.box as B
import numpy as np
import matplotlib.lines as mlines

def fit_legend(fit):
    ax = B.pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    fit_info = (
           f"slope = {fit.slope:.3e} $\pm$ {fit.sigma_s:.3e}\n"
           f"offset = {fit.offset:.3e} $\pm$ {fit.sigma_o:.3e}"
       )
    handles.append(mlines.Line2D([], [], color='None', label=fit_info))
    labels.append(fit_info)
    
    #ax.legend(handles=handles, labels=labels)
    return (handles, labels)

def plot_n_fit(x,y,title=''):
    B.pl.figure()
    plot = B.plot_exp(x, y, marker='+',markersize=10,label=title)
    lfit = B.linefit(x[:4], y[:4])
    lfit.plot()
    B.pl.title(title)
    han, lab = fit_legend(lfit)
    B.pl.legend(handles=han,labels=lab,loc='best',fontsize='large')

    corr = [lfit(-8),lfit(-4),lfit(0),lfit(4),y[4],y[5]]
    return plot, corr

Eb = 10.542
kf_nom = 8.55 #[GeV]
pf_nom = np.array([3.499,3.145,2.783,2.417])
the_nom = np.array([0.2469,0.2258,0.2042,0.1820,0.1592,0.1344])
the_nom_deg = the_nom*180/np.pi

f = B.get_file('./offsets.txt')

dp = B.get_data(f, 'delta')
dthp_0 = B.get_data(f, 'step0_thp')
dkf_0 = B.get_data(f, 'step0_kf')
dPf_0 = B.get_data(f, 'step0_Pf')

dthp_3 = B.get_data(f, 'step3_thp')
dkf_3 = B.get_data(f, 'step3_kf')
dPf_3 = B.get_data(f, 'step3_Pf')

dthp_5 = B.get_data(f, 'step5_thp')

pmx6 = B.get_data(f, 'step6_pmx')
pmy6 = B.get_data(f, 'step6_pmy')

# dp_dkf = B.plot_exp(dp, dthp_3, marker='+',markersize=10)

# B.pl.figure()
# dp_dkf = B.plot_exp(dp, dkf_3, marker='+',markersize=10)
# lfit = B.linefit(dp[:4], dkf_3[:4])
# lfit.plot()

# B.pl.figure()
# dp_dkf = B.plot_exp(dp, dPf_3, marker='+',markersize=10)
# lfit = B.linefit(dp[:4], dPf_3[:4])
# lfit.plot()

# B.pl.figure()
# dp_dthp = B.plot_exp(dp, dthp_0, marker='+',markersize=10,
#                      label = 'dthp')
# lfit = B.linefit(dp[:4], dthp_0[:4])
# lfit.plot()
# B.pl.title('dthp')

# B.pl.figure()
# dp_dkf = B.plot_exp(dp, dkf_0, marker='+',markersize=10,
#                     label = 'dkf')
# lfit = B.linefit(dp[:4], dkf_0[:4])
# lfit.plot()
# B.pl.title('dkf')

# B.pl.figure()
# dp_dPf = B.plot_exp(dp, dPf_0, marker='+',markersize=10,
#                     label = 'dPf')
# lfit_Pf = B.linefit(dp[:4], dPf_0[:4])
# lfit_Pf.plot()
# B.pl.title('dPf')
# han, lab = fit_legend(lfit_Pf)
# B.pl.legend(handles=han,labels=lab,loc='best',fontsize='large')

# pf_nom = np.array([3.499,3.145,2.783,2.417])
# corr = [lfit_Pf(-8),lfit_Pf(-4),lfit_Pf(0),lfit_Pf(4)]

# corr_factor = 1-(corr/pf_nom)

# B.pl.figure()
# dp_dPf = B.plot_exp(dp, dthp_5, marker='+',markersize=10,
#                     label = 'dthp')
# lfit_thp = B.linefit(dp[:4], dthp_5[:4])
# lfit_thp.plot()
# B.pl.title('dthp')
# han, lab = fit_legend(lfit_thp)
# B.pl.legend(handles=han,labels=lab,loc='best',fontsize='large')

# corr = [lfit_thp(-8),lfit_thp(-4),lfit_thp(0),lfit_thp(4),dthp_5[4],dthp_5[5]]

# dthp_dthe_1 = kf_nom*(kf_nom - Eb*np.cos(the_nom))
# dthp_dthe_2 = Eb**2 - 2.*Eb*kf_nom*np.cos(the_nom) + kf_nom**2
# dthp_dthe = dthp_dthe_1/dthp_dthe_2

# corr_offset = corr/dthp_dthe

#### pmx, pmy

pmx_plot, pmx_corr = plot_n_fit(dp, pmx6, title='Pmx')
pmy_plot, pmy_corr = plot_n_fit(dp, pmy6, title='Pmy')

th_corr = pmx_corr[:4]/pf_nom

ph_corr = pmy_corr[:4]/pf_nom