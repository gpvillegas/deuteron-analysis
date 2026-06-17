#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:05:19 2026

@author: gvill
"""

import LT.box as B
import pickle as P
import bin_info2 as BI
import numpy as np

plot_setting = input('Select setting to plot:\npm_120, pm_580, pm_800, pm_900\n')

filename = f'yield_n_xsec/{plot_setting}_histos_comm_bin.pcl'
with open(filename, 'rb') as f:
    histos = P.load(f)

theory_filename = f'yield_n_xsec/jml_theory_xsec/{plot_setting}_jml_theory.txt'
ft = B.get_file(theory_filename)

meas_filename = f'yield_n_xsec/average_xsec/{plot_setting}_jml_average_xsec_set_1.txt'
fm = B.get_file(meas_filename)

red_filename = f'yield_n_xsec/bin_centering/{plot_setting}_jml_bc_corr_set_1.txt'
fr = B.get_file(red_filename)

xb = ft['xb'] # th_nq
yb = ft['yb'] # pm    

sel = (xb==45.) & (yb==0.58)
model_Xsec_atAvgKin = ft['fsiXsec'][sel]

y_data_ecorr = histos[f'{plot_setting}_yield_norm']['set_1']
y_model_rad = histos[f'{plot_setting}_JML_Yield_FSI_rad']

yd_ecorr_bi = BI.get_histo_data_arrays(y_data_ecorr, LT_histo=True)
ym_rad_bi = BI.get_histo_data_arrays(y_model_rad, LT_histo=True)

yield_data_ecorr = yd_ecorr_bi.cont[sel]
yield_model_rad = ym_rad_bi.cont[sel]

Xsec_vertexCorr_JRmethod = model_Xsec_atAvgKin*(yield_data_ecorr/yield_model_rad)

#%%
pm_sett = 'pm_580'

avgKin_filename = f'yield_n_xsec/average_kin/{pm_sett}_jmlfsi_norad_avgkin.txt'
fak = B.get_file(avgKin_filename)

yero_avgKin_filename = 'yield_n_xsec/pm580_fsi_norad_avgkin_set1.txt'
fak_yero = B.get_file(yero_avgKin_filename)

sel = (fak['xb']==45.) & (fak['yb']>0.58)

print('Avg Kin Comparison for thnq = 45. pm = 0.58')
print('pm bin = ', fak['yb'][sel])
print('Ei (mine) = ', fak['Ei'][sel])
print('Ei (yero) = ', fak_yero['Ei'][sel])
print('Q2 (mine) = ', fak['Q2'][sel])
print('Q2 (yero) = ', fak_yero['Q2'][sel])
print('omega (mine) = ', fak['omega'][sel])
print('omega (yero) = ', fak_yero['omega'][sel])
print('theta_pq (mine) = ', fak['theta_pq'][sel])
print('theta_pq (yero) = ', fak_yero['theta_pq'][sel])
print('phi_pq (mine) = ', np.acos(fak['cos_phi'][sel])*180/np.pi)
print('phi_pq (yero) = ', np.acos(fak_yero['cos_phi'][sel])*180/np.pi)
print('phi_pq (histo) = ', fak['phi_pq_histo'][sel])

# how to get the average of angles:
# https://scientificgems.wordpress.com/2014/08/18/mathematics-in-action-averaging-angles/

phi_pq = np.atan2(fak['sin_phi'][sel],fak['cos_phi'][sel]) 
print('phi_pq (new avg) = ', phi_pq*180/np.pi)  

