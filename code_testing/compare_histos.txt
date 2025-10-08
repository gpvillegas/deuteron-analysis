#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 16:02:48 2024

analysis his 2d histos

@author: boeglinw
"""
import numpy as np
import LT.box as B


#%%
hist_dir = './histos/'

run2 = 20871
run1 = 20872
exp_file1 = f'EmPm_run{run1}.data'
exp_file2 = f'EmPm_run{run2}.data'


#%% load histos

# exp1 data
h2_exp1 = B.histo2d(file = hist_dir + exp_file1)
h2_exp1.xlabel='Pm'
h2_exp1.ylabel='Em'
h2_exp1.title = 'Exp. Data'

# exp2 data
h2_exp2 = B.histo2d(file = hist_dir + exp_file2)
h2_exp2.xlabel='Pm'
h2_exp2.ylabel='Em'
h2_exp2.title = 'Exp. Data'

#%% setup projections and plot

# xbin number

bin_range = np.arange(10,28)


Em_min = -0.05
Em_max_a = np.array([0.05, 0.07, 0.09, 0.11])
Em_max =  0.09

B.pl.close('all')
plot_all = False
plot_separate = False

#%% perfomr projections and integrate
def Em_integration(bin_range, Emin = Em_min, Emax = Em_max, plot_all = False):
    counts_exp1 = []
    counts_exp2 = []
    for nx in bin_range[:]:
        hexp1 = h2_exp1.project_y(bins = [nx])
        hexp2 = h2_exp2.project_y(bins = [nx])
        
        counts_exp1.append(hexp1.sum(Emin, Emax))
        counts_exp2.append(hexp2.sum(Emin, Emax))
        
        if plot_all:
            B.pl.figure()
            hexp1.plot(filled = False)
            hexp2.plot(filled = False)
        
    counts_exp1 = np.array(counts_exp1)
    counts_exp2 = np.array(counts_exp2)
    
    c1 = counts_exp1[:,0]
    dc1 = counts_exp1[:,1]
    
    c2 = counts_exp2[:,0]
    dc2 = counts_exp2[:,1]
    
    # calculate ratios
    
    pm_values = h2_exp1.x_bin_center[bin_range]
    
    R = c1/c2
    
    # assuming the errors in the exp. histos are correct
    dR = np.sqrt((dc1/c2)**2 + (c1*dc2/c2**2)**2)
    
    return pm_values, R, dR
    
#%% 

for Em_max in Em_max_a:
    pm_values, R, dR = Em_integration(bin_range, Em_min, Em_max, plot_all = False)
    B.plot_exp(pm_values, R, dR, label = f'Em_max = {Em_max:.2e} (GeV)')
    
B.pl.legend()
B.pl.xlabel('Pm (GeV/c)')
B.pl.ylabel('Y_exp/Y_simc')
B.pl.title(f'Yield Ratios for {run1}/{run2}')


B.pl.savefig(f'Yield_ratios_{run1}_{run2}.pdf')
