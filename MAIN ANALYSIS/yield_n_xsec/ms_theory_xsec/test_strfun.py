#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 16:10:37 2025

test interpolation using respy:
    
    2 identical module with different names have been created:
        - respy_pwia  : used to load PWIA grid
        - respy_fsi : used to load fsi grid

@author: boeglinw
"""

import numpy as np
import LT.box as B

from respy_pwia import response_functions_mod as RSM_pwia
from respy_fsi import response_functions_mod as RSM_fsi

#%% helper functions
def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    B.pl.legend(fontsize = fsize)

#%%
dtr = np.pi/180


pm  = 0.938272
dm  = 1.875612
pmn = 0.939565
eb  = 0.00221

dtr = np.pi/180.

#%% solve quadratic eq.

def quadsol(a,b,c):
    # find locations where a is 0
    a_0 = a == 0.
    a_ok = ~a_0
    b_0 = b ==0
    b_ok = ~b_0
    discr = b**2 - 4*a*c
    x1 = np.zeros_like(a)
    x2 = np.zeros_like(a)
    # handle case where a is 0
    x1[a_0] = -(b/c)[a_0]
    x2[a_0] = x1[a_0]
    # handle case where b is 0
    x1[b_0] = -np.sqrt(-c[b_0]/a[b_0]) 
    x2[b_0] = np.sqrt(-c[b_0]/a[b_0]) 
    # handle the rest
    q = -0.5 * (b + np.sign(b)*np.sqrt(discr))
    q_0 = q==0
    q_ok = ~q_0
    all_ok = q_ok & a_ok & b_ok
    x1[all_ok] = q[all_ok]/a[all_ok]
    x2[all_ok] =c[all_ok]/q[all_ok]    
    return x1,x2

def quad_eval(a,b,c,x):
    return a*x**2 + b*x + c

#%%

def calc_kin (pr, thr):
    er = np.sqrt(pmn**2 + pr**2)
    er_tilde = np.sqrt(pm**2 + pr**2)

    delta2 = (dm - er)**2 - er_tilde**2

    alpha = 2.*pr*np.cos(thr)

    a = alpha**2 - 4*delta2 - 4*er_tilde**2
    b = 2*alpha*(delta2 - Q2)

    c = np.ones_like(a)*((delta2 + Q2)**2 + 4*er_tilde**2*Q2)

    qv1,qv2 = quadsol(a, b, c)

    # pick the solution that gives a positive q ma

    is_pos1 = qv1> 0.
    is_pos2 = qv2> 0.

    qvf = np.zeros_like(qv1)
    qvf.fill(np.nan)

    qvf[is_pos1] = qv1[is_pos1]   
    qvf[is_pos2] = qv2[is_pos2]

    nu = np.sqrt(qvf**2 - Q2)
    return nu

#%%
# PWIA
#resp_file_pwia = 'strfun_grid_CD-Bonn_PWIA.data'
resp_file_pwia = 'yield_n_xsec/ms_theory_xsec/grid/PAR/strfun_grid_1_1_60_50_72.data'


# FSI
resp_file_fsi = 'strfun_grid_CD-Bonn_FSI_old.data'

#%% initialize grid

iret = RSM_pwia.resp_initialize(resp_file_pwia)
iret1 = RSM_fsi.resp_initialize(resp_file_fsi)


#%% check interpolation

pr = np.linspace(0.01, 1., 1000)

q2_loc = 4.4 
thr_loc = 75.

#nu_loc = calc_kin(pr, thr_loc)

#%%
resp_pwia = np.array([RSM_pwia.resp_interp(prv, q2_loc, thr_loc) for prv in pr])
resp_fsi = np.array([RSM_fsi.resp_interp(prv, q2_loc, thr_loc) for prv in pr])

#%%
rl_pwia = resp_pwia[:,0]; rt_pwia = resp_pwia[:,1]; rlt_pwia = resp_pwia[:,2]; rtt_pwia = resp_pwia[:,3]
rl_fsi = resp_fsi[:,0]; rt_fsi = resp_fsi[:,1]; rlt_fsi = resp_fsi[:,2]; rtt_fsi = resp_fsi[:,3]

#%%
# remove zeroes (except first element)
def no_zero(arr):
    nozero_arr = arr[1:]
    nozero_arr = list(nozero_arr[nozero_arr > 0.])
    nozero_arr.insert(0,arr[0])
    
    return np.array(nozero_arr)

q2i = np.where(RSM_pwia.x1==q2_loc)[0][0]
thi = np.where(RSM_pwia.x2==thr_loc)[0][0]
prmax = RSM_pwia.max_x3

true_pr_pwia = no_zero(np.array([RSM_pwia.x3[i][q2i][thi] for i in range(prmax)]))

true_rl_pwia = no_zero(np.array([RSM_pwia.rl[i][q2i][thi] for i in range(prmax)]))
true_rt_pwia = no_zero(np.array([RSM_pwia.rt[i][q2i][thi] for i in range(prmax)]))
true_rlt_pwia = no_zero(np.array([RSM_pwia.rlt[i][q2i][thi] for i in range(prmax)]))
true_rtt_pwia = no_zero(np.array([RSM_pwia.rtt[i][q2i][thi] for i in range(prmax)]))

q2i = np.where(RSM_fsi.x1==q2_loc)[0][0]
thi = np.where(RSM_fsi.x2==thr_loc)[0][0]
prmax = RSM_fsi.max_x3

true_pr_fsi = no_zero(np.array([RSM_fsi.x3[i][q2i][thi] for i in range(prmax)]))
pr_rtt = np.array([RSM_fsi.x3[i][q2i][thi] for i in range(prmax)])

true_rl_fsi = no_zero(np.array([RSM_fsi.rl[i][q2i][thi] for i in range(prmax)]))
true_rt_fsi = no_zero(np.array([RSM_fsi.rt[i][q2i][thi] for i in range(prmax)]))
true_rlt_fsi = no_zero(np.array([RSM_fsi.rlt[i][q2i][thi] for i in range(prmax)]))
true_rtt_fsi = np.array([RSM_fsi.rtt[i][q2i][thi] for i in range(prmax)])


# compare momentum distributions
#%%
B.pl.figure()
B.pl.scatter(true_pr_pwia, true_rl_pwia,
             c='#7b8fd4',s=84,marker='+',label='PWIA input')
B.pl.scatter(true_pr_fsi, true_rl_fsi,
             c='#ed6792',s=84,marker='+',label='PWIA+FSI input')
B.pl.plot(pr, rl_pwia,c='#7b8fd4',label='PWIA interp')
B.pl.plot(pr, rl_fsi,c='#ed6792', label='PWIA+FSI interp')
B.pl.yscale('log')
#B.pl.title(r'$R_L\ \theta = $' +f'{thr_loc}')
pptxify(x='$p_r$',y='$R_L$',t=f'$R_L$\t$\\theta = {thr_loc}$\t$Q^2 = {q2_loc}$')

B.pl.figure()
B.pl.scatter(true_pr_pwia, true_rt_pwia,
             c='#7b8fd4',s=84,marker='+',label='PWIA input')
B.pl.scatter(true_pr_fsi, true_rt_fsi,
             c='#ed6792',s=84,marker='+',label='PWIA+FSI input')
B.pl.plot(pr, rt_pwia,c='#7b8fd4',label='PWIA interp')
B.pl.plot(pr, rt_fsi,c='#ed6792', label='PWIA+FSI interp')
B.pl.yscale('log')
# B.pl.title(r'$R_T\ \theta = $' +f'{thr_loc}')
pptxify(x='$p_r$',y='$R_T$',t=f'$R_T$\t$\\theta = {thr_loc}$\t$Q^2 = {q2_loc}$')

B.pl.figure()
B.pl.scatter(true_pr_pwia, true_rlt_pwia,
             c='#7b8fd4',s=84,marker='+',label='PWIA input')
B.pl.scatter(true_pr_fsi, true_rlt_fsi,
             c='#ed6792',s=84,marker='+',label='PWIA+FSI input')
B.pl.plot(pr, rlt_pwia,c='#7b8fd4',label='PWIA interp')
B.pl.plot(pr, rlt_fsi,c='#ed6792', label='PWIA+FSI interp')
B.pl.yscale('log')
# B.pl.title(r'$R_{LT}\ \theta = $' +f'{thr_loc}')
pptxify(x='$p_r$',y='$R_{LT}$',t='$R_{LT}$'+f'\t$\\theta = {thr_loc}$\t$Q^2 = {q2_loc}$')

B.pl.figure()
B.pl.scatter(true_pr_pwia, true_rtt_pwia,
             c='#7b8fd4',s=84,marker='+',label='PWIA input')
B.pl.scatter(pr_rtt, true_rtt_fsi,
             c='#ed6792',s=84,marker='+',label='PWIA+FSI input')
B.pl.plot(pr, rtt_pwia,c='#7b8fd4',label='PWIA interp')
B.pl.plot(pr, rtt_fsi,c='#ed6792', label='PWIA+FSI interp')
B.pl.yscale('log')
# B.pl.title(r'$R_{TT}\ \theta = $' +f'{thr_loc}')
pptxify(x='$p_r$',y='$R_{TT}$',t='$R_{TT}$'+f'\t$\\theta = {thr_loc}$\t$Q^2 = {q2_loc}$')
#%% angular distributions

pm_set = 0.8

theta_a = np.linspace(1,179, 178)


nu_loc = calc_kin(pm_set, theta_a*dtr)

is_ok = ~np.isnan(nu_loc)

theta = theta_a[is_ok]

pt = pm_set*np.sin(theta*dtr)

#%%
resp_t_pwia = np.array([RSM_pwia.resp_interp(pm_set, q2_loc, thv) for thv in theta])
resp_t_fsi = np.array([RSM_fsi.resp_interp(pm_set, q2_loc, thv) for thv in theta])

rl_th_pwia = resp_t_pwia[:,0]; rt_th_pwia = resp_t_pwia[:,1]; rlt_th_pwia = resp_t_pwia[:,2]; rtt_th_pwia = resp_t_pwia[:,3]
rl_th_fsi = resp_t_fsi[:,0]; rt_th_fsi = resp_t_fsi[:,1]; rlt_th_fsi = resp_t_fsi[:,2]; rtt_th_fsi = resp_t_fsi[:,3]


sel = rl_th_pwia > 0.
#%%
B.pl.figure()
B.pl.plot(theta[sel], rl_th_pwia[sel])
B.pl.plot(theta[sel], rl_th_fsi[sel])
B.pl.title(r'$R_L\ p_m = $' +f'{pm_set}')

B.pl.figure()
B.pl.plot(theta[sel], rt_th_pwia[sel])
B.pl.plot(theta[sel], rt_th_fsi[sel])
B.pl.title(r'$R_T\ p_m = $' +f'{pm_set}')

B.pl.figure()
B.pl.plot(theta[sel], rlt_th_pwia[sel])
B.pl.plot(theta[sel], rlt_th_fsi[sel])
B.pl.title(r'$R_{LT}\ p_m = $' +f'{pm_set}')

B.pl.figure()
B.pl.plot(theta[sel], rtt_th_pwia[sel])
B.pl.plot(theta[sel], rtt_th_fsi[sel])
B.pl.title(r'$R_{TT}\ p_m = $' +f'{pm_set}')

#%% calc cross sections
# calculate energy transfer


theta_r = theta*dtr
Q2 = q2_loc


er = np.sqrt(pmn**2 + pm_set**2)
er_tilde = np.sqrt(pm**2 + pm_set**2)

delta = dm - er
delta2 = delta**2 - er_tilde**2

alpha = 2.*pm_set*np.cos(theta_r)

an = 4*delta**2 - alpha**2
bn = 4*delta*(delta2 - Q2)*np.ones_like(an)
cn = (delta2 - Q2)**2 - alpha**2*Q2

nv1,nv2 = quadsol(an, bn, cn)

test1 = (2*delta*nv1 + (delta2 -Q2))/(-alpha*np.sqrt(nv1**2 + Q2))
test2 = (2*delta*nv2 + (delta2 -Q2))/(-alpha*np.sqrt(nv2**2 + Q2))

sel1 = (test1 > 0) & (nv1 > 0.)
sel2 = (test2 > 0) & (nv2 > 0.)

nuv = np.zeros_like(an)
nuv.fill(np.nan)

nuv[sel1] = nv1[sel1]
nuv[sel2] = nv2[sel2]



#%% calculate electron kinematics given q2, pr, thr



er = np.sqrt(pmn**2 + pm_set**2)
er_tilde = np.sqrt(pm**2 + pm_set**2)

delta2 = (dm - er)**2 - er_tilde**2

alpha = 2.*pm_set*np.cos(theta_r)

a = alpha**2 - 4*delta2 - 4*er_tilde**2
b = 2*alpha*(delta2 - Q2)

c = np.ones_like(a)*((delta2 + Q2)**2 + 4*er_tilde**2*Q2)

qv1,qv2 = quadsol(a, b, c)




# pick the solution that gives a positive q ma

is_pos1 = qv1> 0.
is_pos2 = qv2> 0.

qvf = np.zeros_like(qv1)
qvf.fill(np.nan)

qvf[is_pos1] = qv1[is_pos1]   
qvf[is_pos2] = qv2[is_pos2]

nu = np.sqrt(qvf**2 - Q2)

#%%
E_i = 10.2
E_f = E_i - nu

# pick only goo values: no nan's and no negative E_f's
sel = (E_f>0.) & (~np.isnan(E_f))


# calculate electron scattering angle
the = 2*np.arcsin(np.sqrt(q2_loc/(4*E_i*E_f[sel])))


phi_r = 0.
  

Ef = E_f[sel]
theta_nq = theta[sel] 


#%% calculate cross sections

sig_res_pwia = np.array([RSM_pwia.resp_calc_sigma(E_i, ee_f, the[i], pm_set, phi_r) for i,ee_f in enumerate(Ef)]) 
sig_res_fsi = np.array([RSM_fsi.resp_calc_sigma(E_i, ee_f, the[i], pm_set, phi_r) for i,ee_f in enumerate(Ef)])


sigma_pwia = sig_res_pwia[:,0]
sigma_fsi = sig_res_fsi[:,0]

#%%

B.pl.figure()
B.pl.plot(theta_nq, sigma_pwia)
B.pl.plot(theta_nq, sigma_fsi)
B.pl.title(r'$\sigma = $' +f'{pm_set}')


#%% ratio

R = sigma_fsi/sigma_pwia

B.pl.figure()
B.pl.plot(theta_nq, R)
B.pl.title(r'$\sigma_{fsi}/\sigma_{pwia}\ p_m = $' +f'{pm_set}')
B.pl.hlines(1., theta.min(), theta.max(), color = 'r', ls = '--')
#B.pl.ylim((0.5, 3.5))
