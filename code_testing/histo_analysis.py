#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 16:02:48 2024

analysis his 2d histos

@author: boeglinw
"""
import numpy as np
import LT.box as B

def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    B.pl.legend(fontsize = fsize)

#%%
hist_dir = './'

run = 20871
setting = 1
exp_file = f'EmPm_run{run}.data'
simc_file = f'EmPm_SIMC_{setting}.data'


#%% load histos

# exp data
h2_exp = B.histo2d(file = hist_dir + exp_file)
h2_exp.xlabel='Pm'
h2_exp.ylabel='Em'
h2_exp.title = 'Exp. Data'

# SIMC data
h2_simc = B.histo2d(file = hist_dir + simc_file)
h2_simc.xlabel='Pm'
h2_simc.ylabel='Em'
h2_simc.title = 'SIMC Data'

#%% setup projections and plot

# xbin number

bin_range = np.arange(10,28)


Em_min = -0.05
Em_max_a120 = np.array([0.11,0.09,0.07,0.05])
Em_max =  0.09

Em_min_580 = 0.3
Em_max_a580 = np.array([0.75, 0.85, 1.0, 1.5, 2.5])
Em_max_580 = 2.5

Em_min_800 = 0.3
Em_max_a800 = np.array([0.75, 0.85, 1.0, 1.5, 2.5])
Em_max_800 = 2.5

Em_min_900 = 0.3
Em_max_a900 = np.array([0.75, 0.85, 1.0, 1.5, 2.5])
Em_max_900 = 2.5


B.pl.close('all')
plot_all = True
plot_separate = False

#%% perfomr projections and integrate
def Em_integration(bin_range, Emin = Em_min, Emax = Em_max, 
                   plot_all = False):
    counts_exp = []
    counts_simc = []
    for nx in bin_range[:]:
        hexp = h2_exp.project_y(bins = [nx])
        hsimc = h2_simc.project_y(bins = [nx])
        
        counts_exp.append(hexp.sum(Emin, Emax))
        counts_simc.append(hsimc.sum(Emin, Emax))
        
        if plot_all:
            B.pl.figure()
            hexp.plot(filled = False)
            hsimc.plot(filled = False)
        
    counts_exp = np.array(counts_exp)
    counts_simc = np.array(counts_simc)
    
    # calculate ratios
    
    pm_values = h2_exp.x_bin_center[bin_range]
    
    R = counts_exp[:,0]/counts_simc[:,0]
    
    # assuming the errors in the exp. histo are correct
    dR = counts_exp[:,1]/counts_simc[:,0]
    
    return pm_values, R, dR
    
#%% 
colors = ['#ffa600','#ff5f6b','#d953b8','#3d6bd4']
i=0
for Em_max in Em_max_a120:
    pm_values, R, dR = Em_integration(bin_range, Em_min, Em_max, plot_all = False)
    B.plot_exp(pm_values, R, dR, label = f' $E_m$ max = {Em_max:.2f} (GeV)',
               c=colors[i],marker='+',markersize=16,
               capsize=0,mew=2,elinewidth=3)
    i+=1
    
# B.pl.legend()
# B.pl.xlabel('Pm (GeV/c)')
# B.pl.ylabel('Y_exp/Y_simc')
# B.pl.title(f'Yield Ratios for {run}')
pptxify('Yield Ratios for $p_m = 120$ MeV setting',
        '$P_m$ (GeV/c)','$Y_{exp}/Y_{SIMC}$')

#B.pl.savefig(f'Yield_ratios_{run}.pdf')
