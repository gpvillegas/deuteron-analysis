#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 10:29:05 2026

@author: gvill

This code gets selected kinematics, hitograms them and fits them to get the 
means and std deviations for use in offset optimization code

"""
import LT.box as B
import data_init as D
import numpy as np

# heep rootfiles locations
DATA_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/pass_3/"
SIMC_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/heep_alloffsets/"

# data and simc branch names
DATA_branches = ['P.kin.primary.W','H.kin.secondary.emiss',
                 'H.kin.secondary.pmiss_x','H.kin.secondary.pmiss_z']
SIMC_branches = ['W','Em','Pmx','Pmz']

# load data and simc rootfiles
data = D.DATA_INIT(data_type='deut23_data',
                         kin_study='heep_coin',
                         select_branches={'T':DATA_branches},
                         ROOTfiles_path=DATA_baseDIR)
simc = D.DATA_INIT(data_type='SIMC',
                         kin_study='heep_coin',
                         select_branches={'SNT':SIMC_branches},
                         SIMC_ROOTfiles_path=SIMC_baseDIR,
                         simc_type='-')

#%%
# write results to a file for DATA
header = """# Means and Std. Deviations for Kinematic variables:
# W, Em, Pmx, Pmz
# From DATA
#! Run[i,0]/ W_mean[f,1]/ W_mean_err[f,2]/ W_sigma[f,3]/ W_sigma_err[f,4]/ Em_mean[f,5]/ Em_mean_err[f,6]/ Em_sigma[f,7]/ Em_sigma_err[f,8]/ Pmx_mean[f,9]/ Pmx_mean_err[f,10]/ Pmx_sigma[f,11]/ Pmx_sigma_err[f,12]/ Pmz_mean[f,13]/ Pmz_mean_err[f,14]/ Pmz_sigma[f,15]/ Pmz_sigma_err[f,16]/
"""

fname = 'SYSTEMATIC STUDIES/kin_systematics/kin_DATA_info3.txt'
f = open(fname,'w')
f.write(header)

# get kinematic variables, histogram and fit, get means and sigma
for m in [20841, 20846, 20851, 20858, 20861]:
    
    # create DATA kinematic histograms
    W_h = B.histo(data.Branches[m]['P.kin.primary.W'],
                  range=(0.8,1.2),bins=100,
                  title='Invariant Mass',xlabel='W [GeV]',ylabel='')
    Em_h = B.histo(data.Branches[m]['H.kin.secondary.emiss'],
                   range=(-0.1,0.1),bins=100,
                   title='Missing Energy',xlabel='Em [GeV]',ylabel='')
    Pmx_h = B.histo(data.Branches[m]['H.kin.secondary.pmiss_x'],
                    range=(-0.1,0.1),bins=100,
                    title='Mising Momentum X Component',xlabel='Pmx [GeV]',ylabel='')
    Pmz_h = B.histo(data.Branches[m]['H.kin.secondary.pmiss_z'],
                    range=(-0.1,0.1),bins=100,
                    title='Missing Momentum Z Component',xlabel='Pmz [GeV]',ylabel='')
    
    # fit the histograms with a gaussian
    W_h.fit(plot_fit=False)
    mu = W_h.mean.value
    s = W_h.sigma.value
    W_h.fit(mu-0.6*s,mu+0.5*s,plot_fit=False)
    
    Em_h.fit(plot_fit=False)
    mu = Em_h.mean.value
    s = Em_h.sigma.value
    if m == 20841:
        Em_h.fit(mu-0.8*s,mu+0.5*s,plot_fit=False)
    else:
        Em_h.fit(mu-s,mu+0.5*s,plot_fit=False)
    
    Pmx_h.fit(plot_fit=False)
    mu = Pmx_h.mean.value
    s = Pmx_h.sigma.value
    Pmx_h.fit(mu-0.5*s,mu+0.8*s,plot_fit=False)
    
    Pmz_h.fit(plot_fit=False)
    mu = Pmz_h.mean.value
    s = Pmz_h.sigma.value
    Pmz_h.fit(mu-0.5*s,mu+0.8*s,plot_fit=False)
    
    # show plots
    h_to_plot = [W_h,Em_h,Pmx_h,Pmz_h]
    B.pl.figure(figsize=(10,10))
    B.pl.suptitle(f'Run {m}\n'+\
                  f'Chi2/DoF W, Em, Pmx, Pmz = {W_h.chi2_red:.2f} {Em_h.chi2_red:.2f} {Pmx_h.chi2_red:.2f} {Pmz_h.chi2_red:.2f}')
    i=1
    for h in h_to_plot:
        B.pl.subplot(2,2,i)
        h.plot()
        h.plot_fit()
        i+=1
        
    # write means and sigma to file
    f.write(f'{m} {W_h.mean.value:0.4e} {W_h.mean.err:0.4e} {W_h.sigma.value:0.4e} {W_h.sigma.err:0.4e}'+\
            f' {Em_h.mean.value:0.4e} {Em_h.mean.err:0.4e} {Em_h.sigma.value:0.4e} {Em_h.sigma.err:0.4e}'+\
            f' {Pmx_h.mean.value:0.4e} {Pmx_h.mean.err:0.4e} {Pmx_h.sigma.value:0.4e} {Pmx_h.sigma.err:0.4e}'+\
            f' {Pmz_h.mean.value:0.4e} {Pmz_h.mean.err:0.4e} {Pmz_h.sigma.value:0.4e} {Pmz_h.sigma.err:0.4e}\n')
    
    print(f'Chi2/DoF W, Em, Pmx, Pmz = {W_h.chi2_red:.2f} {Em_h.chi2_red:.2f} {Pmx_h.chi2_red:.2f} {Pmz_h.chi2_red:.2f}')
f.close()

#%%
header = """
# Means and Std. Deviations for Kinematic variables:
# W, Em, Pmx, Pmz
# From SIMC
#! Setting[s,0]/ W_sim_mean[f,1]/ W_sim_mean_err[f,2]/ W_sim_sigma[f,3]/ W_sim_sigma_err[f,4]/ Em_sim_mean[f,5]/ Em_sim_mean_err[f,6]/ Em_sim_sigma[f,7]/ Em_sim_sigma_err[f,8]/ Pmx_sim_mean[f,9]/ Pmx_sim_mean_err[f,10]/ Pmx_sim_sigma[f,11]/ Pmx_sim_sigma_err[f,12]/ Pmz_sim_mean[f,13]/ Pmz_sim_mean_err[f,14]/ Pmz_sim_sigma[f,15]/ Pmz_sim_sigma_err[f,16]/
"""

fname = 'SYSTEMATIC STUDIES/kin_systematics/kin_SIMC_info3.txt'
f = open(fname,'w')
f.write(header)

for s in ['delta_scan_-8','delta_scan_-4','delta_scan_0','delta_scan_+4','delta_scan_+8']:
    # create SIMC kinematic histograms
    W_h = B.histo(simc.Branches[s]['W'],
                  range=(0.8,1.2),bins=100,
                  title='Invariant Mass',xlabel='W [GeV]',ylabel='')
    Em_h = B.histo(simc.Branches[s]['Em'],
                   range=(-0.1,0.1),bins=100,
                   title='Missing Energy',xlabel='Em [GeV]',ylabel='')
    Pmx_h = B.histo(simc.Branches[s]['Pmx'],
                    range=(-0.1,0.1),bins=100,
                    title='Mising Momentum X Component',xlabel='Pmx [GeV]',ylabel='')
    Pmz_h = B.histo(simc.Branches[s]['Pmz'],
                    range=(-0.1,0.1),bins=100,
                    title='Missing Momentum Z Component',xlabel='Pmz [GeV]',ylabel='')
    
    # fit the histograms with a gaussian
    W_h.fit(plot_fit=False)
    mu = W_h.mean.value
    si = W_h.sigma.value
    W_h.fit(mu-0.8*si,mu+0.6*si,plot_fit=False)
    
    Em_h.fit(plot_fit=False)
    mu = Em_h.mean.value
    si = Em_h.sigma.value
    if s == 'delta_scan_-8':
        Em_h.fit(mu-0.9*si,mu+0.1*si,plot_fit=False)
    elif s == 'delta_scan_-4':
        Em_h.fit(mu-si,mu,plot_fit=False)
    elif s == 'delta_scan_0':
        Em_h.fit(mu-1.1*si,mu+0.1*si,plot_fit=False)
    elif s == 'delta_scan_+4':
        Em_h.fit(mu-0.7*si,mu+0.7*si,plot_fit=False)
    elif s == 'delta_scan_+8':
        Em_h.fit(mu-0.7*si,mu+0.7*si,plot_fit=False)
    
    Pmx_h.fit(plot_fit=False)
    mu = Pmx_h.mean.value
    si = Pmx_h.sigma.value
    Pmx_h.fit(mu-0.7*si,mu+si,plot_fit=False)
    
    Pmz_h.fit(plot_fit=False)
    mu = Pmz_h.mean.value
    si = Pmz_h.sigma.value
    Pmz_h.fit(mu-0.8*si,mu+0.3*si,plot_fit=False)
    
    # show plots
    h_to_plot = [W_h,Em_h,Pmx_h,Pmz_h]
    B.pl.figure(figsize=(10,10))
    B.pl.suptitle(f'Setting {s}\n'+\
                  f'Chi2/DoF W, Em, Pmx, Pmz = {W_h.chi2_red:.2f} {Em_h.chi2_red:.2f} {Pmx_h.chi2_red:.2f} {Pmz_h.chi2_red:.2f}')

    i=1
    for h in h_to_plot:
        B.pl.subplot(2,2,i)
        h.plot()
        h.plot_fit()
        i+=1
        
    f.write(f'{s} {W_h.mean.value:0.4e} {W_h.mean.err:0.4e} {W_h.sigma.value:0.4e} {W_h.sigma.err:0.4e}'+\
            f' {Em_h.mean.value:0.4e} {Em_h.mean.err:0.4e} {Em_h.sigma.value:0.4e} {Em_h.sigma.err:0.4e}'+\
            f' {Pmx_h.mean.value:0.4e} {Pmx_h.mean.err:0.4e} {Pmx_h.sigma.value:0.4e} {Pmx_h.sigma.err:0.4e}'+\
            f' {Pmz_h.mean.value:0.4e} {Pmz_h.mean.err:0.4e} {Pmz_h.sigma.value:0.4e} {Pmz_h.sigma.err:0.4e}\n')

f.close()

#%%
fname = 'SYSTEMATIC STUDIES/kin_systematics/kin_DATA_info3.txt'
f = B.get_file(fname)

runs = f['Run']

W_mean = f['W_mean']
# W_sigma = f['W_sigma']
W_mean_err = f['W_mean_err']

Em_mean = f['Em_mean']
# Em_sigma = f['Em_sigma']
Em_mean_err = f['Em_mean_err']

Pmx_mean = f['Pmx_mean']
# Pmx_sigma = f['Pmx_sigma']
Pmx_mean_err = f['Pmx_mean_err']

Pmz_mean = f['Pmz_mean']
# Pmz_sigma = f['Pmz_sigma']
Pmz_mean_err = f['Pmz_mean_err']

fname = 'SYSTEMATIC STUDIES/kin_systematics/kin_SIMC_info2.txt'
f = B.get_file(fname)

W_sim_mean = f['W_sim_mean']
# W_sigma = f['W_sim_sigma']
W_sim_mean_err = f['W_sim_mean_err']

Em_sim_mean = f['Em_sim_mean']
# Em_sigma = f['Em_sim_sigma']
Em_sim_mean_err = f['Em_sim_mean_err']

Pmx_sim_mean = f['Pmx_sim_mean']
# Pmx_sigma = f['Pmx_sim_sigma']
Pmx_sim_mean_err = f['Pmx_sim_mean_err']

Pmz_sim_mean = f['Pmz_sim_mean']
# Pmz_sigma = f['Pmz_sim_sigma']
Pmz_sim_mean_err = f['Pmz_sim_mean_err']

print('W_mean =', np.array(W_mean).tolist())
# print('W_sigma =', np.array(W_sigma).tolist())
print('W_mean_err =', np.array(W_mean_err).tolist())
print('W_sim_mean =', np.array(W_sim_mean).tolist())
# print('W_sim_sigma =', np.array(W_sigma).tolist())
print('W_sim_mean_err =', np.array(W_sim_mean_err).tolist())

print('Em_mean =', np.array(Em_mean).tolist())
# print('Em_sigma =', np.array(Em_sigma).tolist())
print('Em_mean_err =', np.array(Em_mean_err).tolist())
print('Em_sim_mean =', np.array(Em_sim_mean).tolist())
# print('Em_sim_sigma =', np.array(Em_sigma).tolist())
print('Em_sim_mean_err =', np.array(Em_sim_mean_err).tolist())

print('Pmx_mean =', np.array(Pmx_mean).tolist())
# print('Pmx_sigma =', np.array(Pmx_sigma).tolist())
print('Pmx_mean_err =', np.array(Pmx_mean_err).tolist())
print('Pmx_sim_mean =', np.array(Pmx_sim_mean).tolist())
# print('Pmx_sim_sigma =', np.array(Pmx_sigma).tolist())
print('Pmx_sim_mean_err =', np.array(Pmx_sim_mean_err).tolist())

print('Pmz_mean =', np.array(Pmz_mean).tolist())
# print('Pmz_sigma =', np.array(Pmz_sigma).tolist())
print('Pmz_mean_err =', np.array(Pmz_mean_err).tolist())
print('Pmz_sim_mean =', np.array(Pmz_sim_mean).tolist())
# print('Pmz_sim_sigma =', np.array(Pmz_sigma).tolist())
print('Pmz_sim_mean_err =', np.array(Pmz_sim_mean_err).tolist())