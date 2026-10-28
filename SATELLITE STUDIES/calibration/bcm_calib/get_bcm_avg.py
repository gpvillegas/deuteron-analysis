#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 11:34:49 2026

@author: gvill
"""

import numpy as np
import LT.box as B
import matplotlib.pyplot as plt
from matplotlib import colormaps as cmp
from matplotlib import rc

fname = 'calibration/bcm_calib/BCM_calib_constants.txt'

f = B.get_file(fname)

BCM1_Gain = np.average(f['BCM1_Gain'])
BCM1_Gain_err = np.sqrt(np.sum(f['BCM1_Gain_err']**2)/2)
BCM1_Offset = np.average(f['BCM1_Offset'])
BCM1_Offset_err = np.sqrt(np.sum(f['BCM1_Offset_err']**2)/2)

BCM2_Gain = np.average(f['BCM2_Gain'])
BCM2_Gain_err = np.sqrt(np.sum(f['BCM2_Gain_err']**2)/2)
BCM2_Offset = np.average(f['BCM2_Offset'])
BCM2_Offset_err = np.sqrt(np.sum(f['BCM2_Offset_err']**2)/2)

BCM4A_Gain = np.average(f['BCM4A_Gain'])
BCM4A_Gain_err = np.sqrt(np.sum(f['BCM4A_Gain_err']**2)/2)
BCM4A_Offset = np.average(f['BCM4A_Offset'])
BCM4A_Offset_err = np.sqrt(np.sum(f['BCM4A_Offset_err']**2)/2)

BCM4C_Gain = np.average(f['BCM4C_Gain'])
BCM4C_Gain_err = np.sqrt(np.sum(f['BCM4C_Gain_err']**2)/2)
BCM4C_Offset = np.average(f['BCM4C_Offset'])
BCM4C_Offset_err = np.sqrt(np.sum(f['BCM4C_Offset_err']**2)/2)

print('BCM Calib Range1 and Range2 AVERAGE')
print('BCM1 Gain = ',BCM1_Gain,'+/-',BCM1_Gain_err)
print('BCM1 Offset = ',BCM1_Offset,'+/-',BCM1_Offset_err)
print('BCM2 Gain = ',BCM2_Gain,'+/-',BCM2_Gain_err)
print('BCM2 Offset = ',BCM2_Offset,'+/-',BCM2_Offset_err)
print('BCM4A Gain = ',BCM4A_Gain,'+/-',BCM4A_Gain_err)
print('BCM4A Offset = ',BCM4A_Offset,'+/-',BCM4A_Offset_err)
print('BCM4C Gain = ',BCM4C_Gain,'+/-',BCM4C_Gain_err)
print('BCM4C Offset = ',BCM4C_Offset,'+/-',BCM4C_Offset_err)

#%% plot BCM1 error model

#Use latex commands (e.g. \textit ot \textbf)
# rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 18
}
plt.rc('font', **font)

fname = 'calibration/bcm_calib/BCM1_error_model.txt'
f = B.get_file(fname)

I = f['I']
# convert relative uncertainty to percent
dI_gain = f['dI_gain']*100
dI_offset = f['dI_offset']*100
dI_wire = f['dI_wire']*100
dI_nonlin = f['dI_nonlin']*100
dI_tot = f['dI_tot']*100

B.pl.figure(figsize=(10,7))
B.plot_line(I,dI_tot,label='Total',color=cmp['Paired'](5))
B.plot_line(I,dI_gain,label='Gain',linestyle='--',color=cmp['Paired'](7))
B.plot_line(I,dI_offset,label='Offset',linestyle='--',color=cmp['Paired'](3))
B.plot_line(I,dI_wire,label='Wire',linestyle='-.',color=cmp['Paired'](9))
B.plot_line(I,dI_nonlin,label='Non-Linearity',linestyle='-.',color=cmp['Paired'](1))

B.pl.ylabel('Relative Uncertainty [%]')
B.pl.xlabel('Current [$\mu$A]')
B.pl.legend()

B.pl.savefig('BCM1_error_model.pdf')