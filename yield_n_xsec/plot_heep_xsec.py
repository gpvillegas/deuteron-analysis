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

#%% load pickle file with histograms
filename = 'yield_n_xsec/pickle/delta_scan_-8_histos.pcl'
with open(filename, 'rb') as f:
    hm8 = P.load(f)

filename = 'yield_n_xsec/pickle/delta_scan_-4_histos.pcl'
with open(filename, 'rb') as f:
    hm4 = P.load(f)

filename = 'yield_n_xsec/pickle/delta_scan_0_histos.pcl'
with open(filename, 'rb') as f:
    h0 = P.load(f)

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

a0 = B.histo(bin_center = h0['heep_xsec'].bin_center,
             bin_content = set_nans(h0['heep_xsec'].bin_content),
             bin_error = set_nans(h0['heep_xsec'].bin_error))
am4 = B.histo(bin_center = hm4['heep_xsec'].bin_center,
             bin_content = set_nans(hm4['heep_xsec'].bin_content),
             bin_error = set_nans(hm4['heep_xsec'].bin_error))
am8 = B.histo(bin_center = hm8['heep_xsec'].bin_center,
             bin_content = set_nans(hm8['heep_xsec'].bin_content),
             bin_error = set_nans(hm8['heep_xsec'].bin_error))

#get world data
world_q2 = w['Q^2']
world_th = w['Theta']
world_sig = w['Sig']
world_dsig = w['dSig']

#get only necessary points
# cut = world_q2 > 3.0

# aworld = B.plot_exp(world_th,world_sig,world_dsig)

a0.plot_exp(ignore_zeros=True)
am4.plot_exp(ignore_zeros=True)
am8.plot_exp(ignore_zeros=True)
B.pl.yscale('log')
# B.pl.xlim(3.3,5.5)
B.pl.title('H(e,e\'p) Cross Section')
B.pl.xlabel('$Q^2$ $[GeV^2]$')
B.pl.ylabel('')

#%% yield ratio

B.pl.figure()

a0 = B.histo(bin_center = h0['heep_ratio'].bin_center,
             bin_content = set_nans(h0['heep_ratio'].bin_content),
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