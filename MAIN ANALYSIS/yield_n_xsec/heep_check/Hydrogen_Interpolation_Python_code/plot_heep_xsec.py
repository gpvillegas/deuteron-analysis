#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 15:02:44 2025

@author: gvill
"""

import LT.box as B

import pickle as P
from matplotlib import colormaps as cmp
# import matplotlib.lines as mlines
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline as cspl
import get_Xsec_err as M
import interpolate_ff as ff

#%% load pickle file with histograms
filename = 'yield_n_xsec/pickle/delta_scan_-8_histos_pass4.pcl'
with open(filename, 'rb') as f:
    hm8 = P.load(f)

filename = 'yield_n_xsec/pickle/delta_scan_-4_histos_pass4.pcl'
with open(filename, 'rb') as f:
    hm4 = P.load(f)

filename = 'yield_n_xsec/pickle/delta_scan_0_histos_pass4.pcl'
with open(filename, 'rb') as f:
    h0 = P.load(f)

filename = 'yield_n_xsec/pickle/delta_scan_+4_histos_pass4.pcl'
with open(filename, 'rb') as f:
    hp4 = P.load(f)    
    
filename = 'yield_n_xsec/pickle/delta_scan_+8_histos_pass4.pcl'
with open(filename, 'rb') as f:
    hp8 = P.load(f)  
    
# World ep cross sections from J.Arrington's paper
filename = 'yield_n_xsec/supplemental_materials/data/proton/World_ep_CrossSections.dat'    
w = B.get_file(filename)

#%%
# plot all delta cross sections in a continuos plot

def set_nans(array,set_nan_to=0.):
    new_array = np.zeros(len(array))
    for i in range(len(array)):
        if array[i] == array[i]:
            new_array[i] = array[i]
        else:
            new_array[i] = set_nan_to
            array[i] = set_nan_to
            continue
    return new_array

B.pl.figure()

a0 = B.histo(bin_center = h0['heep_xsec_W'].bin_center,
             bin_content = set_nans(h0['heep_xsec_W'].bin_content),
             bin_error = set_nans(h0['heep_xsec_W'].bin_error))
am4 = B.histo(bin_center = hm4['heep_xsec_W'].bin_center,
             bin_content = set_nans(hm4['heep_xsec_W'].bin_content),
             bin_error = set_nans(hm4['heep_xsec_W'].bin_error))
am8 = B.histo(bin_center = hm8['heep_xsec_W'].bin_center,
             bin_content = set_nans(hm8['heep_xsec_W'].bin_content),
             bin_error = set_nans(hm8['heep_xsec_W'].bin_error))
ap4 = B.histo(bin_center = hp4['heep_xsec_W'].bin_center,
             bin_content = set_nans(hp4['heep_xsec_W'].bin_content),
             bin_error = set_nans(hp4['heep_xsec_W'].bin_error))
ap8 = B.histo(bin_center = hp8['heep_xsec_W'].bin_center,
             bin_content = set_nans(hp8['heep_xsec_W'].bin_content),
             bin_error = set_nans(hp8['heep_xsec_W'].bin_error))

# aworld = B.plot_exp(,world_sig,world_dsig)
am8.plot_exp(ignore_zeros=True,label='$\delta=-8$')
am4.plot_exp(ignore_zeros=True,label='$\delta=-4$')
a0.plot_exp(ignore_zeros=True,label='$\delta=0$')
ap4.plot_exp(ignore_zeros=True,label='$\delta=+4$')
ap8.plot_exp(ignore_zeros=True,label='$\delta=+8$')
B.pl.legend()
B.pl.yscale('log')
# B.pl.xlim(3.3,5.5)
B.pl.title('H(e,e\'p) Cross Section')
B.pl.xlabel('W [GeV]')
B.pl.ylabel('')

#%% xsec q2
B.pl.figure()

a0 = B.histo(bin_center = h0['heep_xsec_Q2'].bin_center,
             bin_content = set_nans(h0['heep_xsec_Q2'].bin_content),
             bin_error = set_nans(h0['heep_xsec_Q2'].bin_error))
am4 = B.histo(bin_center = hm4['heep_xsec_Q2'].bin_center,
             bin_content = set_nans(hm4['heep_xsec_Q2'].bin_content),
             bin_error = set_nans(hm4['heep_xsec_Q2'].bin_error))
am8 = B.histo(bin_center = hm8['heep_xsec_Q2'].bin_center,
             bin_content = set_nans(hm8['heep_xsec_Q2'].bin_content),
             bin_error = set_nans(hm8['heep_xsec_Q2'].bin_error))
ap4 = B.histo(bin_center = hp4['heep_xsec_Q2'].bin_center,
             bin_content = set_nans(hp4['heep_xsec_Q2'].bin_content),
             bin_error = set_nans(hp4['heep_xsec_Q2'].bin_error))
ap8 = B.histo(bin_center = hp8['heep_xsec_Q2'].bin_center,
             bin_content = set_nans(hp8['heep_xsec_Q2'].bin_content),
             bin_error = set_nans(hp8['heep_xsec_Q2'].bin_error))

#get world data
# get right energy range for deuteron
world_e0 = w['E0']
e0_range = (world_e0 > 1000) & (world_e0 < 1100)

# q2 is x, sort in increasing q2
world_q2 = pd.Series(w['Q^2'])[e0_range].sort_values()
sort_by = world_q2.index

world_sig = pd.Series(w['Sig'])[sort_by]
world_dsig = pd.Series(w['dSig'])[sort_by]

Eb = 10.542
Q2 = bin_center = h0['heep_xsec_Q2'].bin_center
kf = np.linspace(7,9,Q2.size)
the = np.linspace(10,14,Q2.size)

w_xsec, dw_xsec = M.sig(Eb,kf,the,Q2)

# f_world_sig = cspl(world_q2,world_sig)
B.plot_exp(Q2,w_xsec,dw_xsec)
aworld = B.plot_exp(world_q2,world_sig,world_dsig)
am8.plot_exp(ignore_zeros=True,label='$\delta=-8$')
am4.plot_exp(ignore_zeros=True,label='$\delta=-4$')
a0.plot_exp(ignore_zeros=True,label='$\delta=0$')
ap4.plot_exp(ignore_zeros=True,label='$\delta=+4$')
ap8.plot_exp(ignore_zeros=True,label='$\delta=+8$')
B.pl.legend()
B.pl.yscale('log')
# B.pl.xlim(3.3,5.5)
B.pl.title('H(e,e\'p) Cross Section')
B.pl.xlabel('$Q^2$ [$GeV^2$]')
B.pl.ylabel('')

#%% yield ratio

B.pl.figure()

a0 = B.histo(bin_center = h0['heep_ratio_Q2'].bin_center,
             bin_content = set_nans(h0['heep_ratio_Q2'].bin_content),
             bin_error = set_nans(h0['heep_ratio'].bin_error))
am4 = B.histo(bin_center = hm4['heep_ratio'].bin_center,
             bin_content = set_nans(hm4['heep_ratio'].bin_content),
             bin_error = set_nans(hm4['heep_ratio'].bin_error),
             label='$\delta$ -4')
am8 = B.histo(bin_center = hm8['heep_ratio'].bin_center,
             bin_content = set_nans(hm8['heep_ratio'].bin_content),
             bin_error = set_nans(hm8['heep_ratio'].bin_error),
             label='$\delta$ -8')

a0.plot_exp(ignore_zeros=True,label='$\delta$ 0')
am4.plot_exp(ignore_zeros=True)
am8.plot_exp(ignore_zeros=True)

B.pl.ylim(0.5,1.5)
B.pl.title('H(e,e\'p) Yield Ratios')
B.pl.xlabel('$Q^2$ $[GeV^2]$')
B.pl.ylabel('')
B.pl.legend()