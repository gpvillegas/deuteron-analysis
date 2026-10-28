#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 14:11:03 2025

@author: gvill
"""
import LT.box as B
import data_init as D
import numpy as np
from scipy.optimize import minimize

#%% HELPER FUNCTIONS
def rm_nan(arr):
    has_nan = np.isnan(arr)
    clean_arr = arr[~has_nan]
    num_nan = arr.size - clean_arr.size
    
    print(f'{num_nan} nan removed')
    return clean_arr

#%% Declare constants
Eb = 10.542         # GeV
Mp = 0.938272       # GeV
rtd = 180/np.pi     # rad to degree

#%% Load hcana variables
my_branches = ['P.kin.primary.scat_ang_rad','H.kin.secondary.xangle','P.gtr.p',
               'H.gtr.p','H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th',
               'P.gtr.ph','P.kin.primary.W','H.kin.secondary.emiss',
               'H.kin.secondary.pmiss', 'H.kin.secondary.pmiss_x',
               'H.kin.secondary.pmiss_y','H.kin.secondary.pmiss_z',
               'P.kin.primary.omega','P.kin.primary.q3m','P.kin.primary.q_x',
               'P.kin.primary.q_y','P.kin.primary.q_z','P.kin.primary.Q2',
               'H.kin.secondary.ph_xq','P.kin.primary.ph_q']

my_branches_sim = ['theta_e', 'theta_p', 'e_pf', 'h_pf', 'e_delta', 'h_delta',
                   'h_xptar','h_yptar','e_xptar','e_yptar','Weight','Normfac',
                   'W','Em','Pm','Pmx','Pmy','Pmz','Q2','nu',]

#%% load data root files
# make sure files are where they should be!
RUN = D.DATA_INIT(data_type='deut23_data', run = 20851, 
                  select_branches=my_branches)

# RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
#                   select_branches=my_branches)

# load simc root files
SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                   setting='delta_scan_0',
                   select_branches=my_branches_sim)

#%% assign wanted variables 

# kinematic observables
W = RUN.Branches['P.kin.primary.W']

Em = RUN.Branches['H.kin.secondary.emiss']

Pm = RUN.Branches['H.kin.secondary.pmiss']

Pmx = RUN.Branches['H.kin.secondary.pmiss_x']

Pmy = RUN.Branches['H.kin.secondary.pmiss_y']

Pmz = RUN.Branches['H.kin.secondary.pmiss_z']

Q2 = RUN.Branches['P.kin.primary.Q2']

# kinematic observables SIMC
W_sim = SIMC.Branches['W']

Em_sim = SIMC.Branches['Em']

Pm_sim = SIMC.Branches['Pm']

Pmx_sim = SIMC.Branches['Pmx']

Pmy_sim = SIMC.Branches['Pmy']

Pmz_sim = SIMC.Branches['Pmz']

Q2_sim = SIMC.Branches['Q2']

# measured quantities
kf = RUN.Branches['P.gtr.p'] #electron momentum (GeV)

the = RUN.Branches['P.kin.primary.scat_ang_rad'] #e scattering angle (rad)

pf = RUN.Branches['H.gtr.p'] #proton momentum (GeV)

thp = RUN.Branches['H.kin.secondary.xangle'] - the #p scattering angle (rad)

qx = RUN.Branches['P.kin.primary.q_x']

qy = RUN.Branches['P.kin.primary.q_y']

qz = RUN.Branches['P.kin.primary.q_z']

p_oop = RUN.Branches['P.gtr.th']

h_oop = RUN.Branches['H.gtr.th']

# measured quantities SIMC
kf_sim = SIMC.Branches['e_pf']/1e3 # GeV

the_sim = SIMC.Branches['theta_e'] # rad

pf_sim = SIMC.Branches['h_pf']/1e3 # GeV

thp_sim = SIMC.Branches['theta_p'] # rad

e_oop_sim = SIMC.Branches['e_xptar']

h_oop_sim = SIMC.Branches['h_xptar']

#%% plot measured quatities, fit and extract mean for central correction

hkf = B.histo(kf, range=(8.4,8.75), bins=100, title='kf')
hkf.fit(plot_fit=False)

hthe = B.histo(the, range=(0.195,0.215), bins=100, title='the')
hthe.fit(plot_fit=False)

hpf = B.histo(pf, range=(2.6,2.9), bins=100, title='pf')
hpf.fit(plot_fit=False)

hthp = B.histo(thp, range=(0.64,0.70), bins=100, title='thp')
hthp.fit(plot_fit=False)

hh_oop = B.histo(h_oop, range=(-0.06,0.06), bins=100, title='h_oop')
hh_oop.fit(plot_fit=False)

hp_oop = B.histo(p_oop, range=(-0.02,0.02), bins=100, title='p_oop')
hp_oop.fit(plot_fit=False)

h_to_plot = [hkf, hthe, hpf, hthp, hh_oop, hp_oop]

for h in h_to_plot:
    B.pl.figure()
    h.plot()
    h.plot_fit()
    
#%% plot measured quatities, fit and extract mean for central correction
# SIMC quantities

hkf_sim = B.histo(kf_sim, range=(8.4,8.75), bins=100, title='kf')
hkf_sim.fit(plot_fit=False)

hthe_sim = B.histo(the_sim, range=(0.195,0.215), bins=100, title='the')
hthe_sim.fit(plot_fit=False)

hpf_sim = B.histo(pf_sim, range=(2.6,2.9), bins=100, title='pf')
hpf_sim.fit(plot_fit=False)

hthp_sim = B.histo(thp_sim, range=(0.64,0.70), bins=100, title='thp')
hthp_sim.fit(plot_fit=False)

hh_oop_sim = B.histo(h_oop_sim, range=(-0.06,0.06), bins=100, title='h_oop')
hh_oop_sim.fit(plot_fit=False)

he_oop_sim = B.histo(e_oop_sim, range=(-0.02,0.02), bins=100, title='p_oop')
he_oop_sim.fit(plot_fit=False)

h_to_plot = [hkf_sim, hthe_sim, hpf_sim, hthp_sim, hh_oop_sim, he_oop_sim]

for h in h_to_plot:
    B.pl.figure()
    h.plot()
    h.plot_fit()    
#%% get the means of the distributions, these are the values I'm gonna use to get
# central correction factors

# DATA
x1 = hkf.mean.value
x2 = hthe.mean.value
x3 = hpf.mean.value
x4 = hthp.mean.value
x5 = hp_oop.mean.value
x6 = hh_oop.mean.value


#%%
my_vars = [W,Em,Pm,Pmx,Pmy,Pmz,Q2,kf,the,pf,thp,qx,qy,qz]

#remove nans
for v in my_vars:
    v = rm_nan(v)

#%% normalization

NORM = D.get_norm(20851)

WEIGHT = D.calc_weights(SIMC.Branches)
 
#%% calculated first derivatives

kfz = x1*np.cos(x2)

# W first derivatives 
num1 = -Mp + Eb*np.cos(x2) - Eb
denum1 = 2.*np.sqrt(Mp*(Mp + Eb- x1) + Eb*(kfz - x1))

dW1 = num1/denum1

num2 = -Eb*x1*np.sin(x2)

dW2 = num2/denum1

dW3 = 0

dW4 = 0

dW5 = 0

dW6 = 0

# Em first derivatives
dEm1 = -1

dEm2 = 0

dEm3 = x3/np.sqrt(Mp**2 + x3**2)

dEm4 = 0

dEm5 = 0

dEm6 = 0

# Pmx first derivatives
dPmx1 = np.sin(x5)*np.cos(x2)

dPmx2 = -x1*np.sin(x5)*np.sin(x2)

dPmx3 = np.sin(x6)*np.cos(x4)

dPmx4 = -x3*np.sin(x6)*np.sin(x4)

dPmx5 = x1*np.cos(x5)*np.cos(x2)

dPmx6 = x3*np.cos(x6)*np.cos(x4)

# Pmy first derivatives
dPmy1 = np.sin(x5)*np.sin(x2)

dPmy2 = x1*np.sin(x5)*np.cos(x2)

dPmy3 = np.sin(x6)*np.sin(x4)

dPmy4 = x3*np.sin(x6)*np.cos(x4)

dPmy5 = x1*np.cos(x5)*np.sin(x2)

dPmy6 = x3*np.cos(x6)*np.sin(x4)

#%% get difference in the observables from SIMC

hW = NORM*B.histo(W, range=(0.85,1.), bins=100)
hW_sim = B.histo(W_sim, range=(0.85,1.), bins=100, weights = WEIGHT, calc_w2=True)

B.pl.figure()
hW.plot(filled=False)
hW.fit(0.88,0.94)

hW_sim.plot(filled=False)
hW_sim.fit(0.92,0.98)

ax = B.pl.gca()
ax.set_autoscaley_on(True)
B.pl.title('W')

hEm = NORM*B.histo(Em, range=(-0.075,0.05), bins=100)
hEm_sim = B.histo(Em_sim, range=(-0.075,0.05), bins=100, weights = WEIGHT, calc_w2=True)

B.pl.figure()
hEm.plot(filled=False)
hEm.fit(-0.04,-0.02)

hEm_sim.plot(filled=False)
hEm_sim.fit(-0.005,0.007)

ax = B.pl.gca()
ax.set_autoscaley_on(True)
B.pl.title('Em')

hPmx = NORM*B.histo(Pmx, range=(-0.04,0.04), bins=100)
hPmx_sim = B.histo(Pmx_sim, range=(-0.04,0.04), bins=100, weights = WEIGHT, calc_w2=True)

B.pl.figure()
hPmx.plot(filled=False)
hPmx.fit(0.01,0.04)

hPmx_sim.plot(filled=False)
hPmx_sim.fit(-0.02,0.02)

ax = B.pl.gca()
ax.set_autoscaley_on(True)
B.pl.title('Pmx')

hPmy = NORM*B.histo(Pmy, range=(-0.04,0.04), bins=100)
hPmy_sim = B.histo(Pmy_sim, range=(-0.04,0.04), bins=100, weights = WEIGHT, calc_w2=True)

B.pl.figure()
hPmy.plot(filled=False)
hPmy.fit(-0.01,0.02)

hPmy_sim.plot(filled=False)
hPmy_sim.fit(-0.02,0.01)

ax = B.pl.gca()
ax.set_autoscaley_on(True)
B.pl.title('Pmy')
#%% difference to minimize
dW = hW.mean.value - hW_sim.mean.value
dEm = hEm.mean.value - hEm_sim.mean.value
dPmx = hPmx.mean.value - hPmx_sim.mean.value
dPmy = hPmy.mean.value - hPmy_sim.mean.value

#%% constructing the vectors and matrices
dR = np.array([dW, dEm, dPmx, dPmy])

j1 = np.array([dW1, dW2, dW3, dW4, dW5, dW6])
j2 = np.array([dEm1, dEm2, dEm3, dEm4, dEm5, dEm6])
j3 = np.array([dPmx1, dPmx2, dPmx3, dPmx4, dPmx5, dPmx6])
j4 = np.array([dPmy1, dPmy2, dPmy3, dPmy4, dPmy5, dPmy6])

J = np.array([j1, j2, j3, j4])

#%% matrix algebra

J_trans = J.transpose()

J_inv = np.linalg.inv(J_trans @ J)

#%% offset vector

dx = J_inv @ J_trans @ dR

#%% minimization setup
def chi_square_fun(delta):
    res = dR - J @ delta
    res_t = res.transpose()
    chi2 = res_t @ res
    
    return chi2

result = minimize(chi_square_fun,dx)

dx1, dx2, dx3, dx4, dx5, dx6 = result.x