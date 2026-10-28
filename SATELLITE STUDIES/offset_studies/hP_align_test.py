#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:54:02 2025

@author: gvill

HMS check momentum alignment: 
        This code checks for central momentum alignment
        and plots HMS delta vs focal plane variables to 
        check for correlations

        The momentum aligment is done using HMS variables only,
        starting from the HMS angle assumption,
        - Calculate HMS momentum using HMS angle, Beam energy: Pf(Eb,thp)
        - Calculate SHMS momentum components using HMS momentum
        - Calculate SHMS momentum magnitude
        - Calculate SHMS angle using SHMS momentum components
        
        Then compare calculated quantities with HCANA quantities (measured)
        by computing the difference
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import LT.box as B

#import deut_analysis_tools as dat
import database_operations as db
import data_init as D
import cut_handler as C
#%% file locations
deut_db = '/home/gvill/deuteron/deuteron_db/deuteron_db.db' 


## FUNCTION DEFINITIONS
# To get rid of unhelpful groupings (lists of lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    return [k[index] for k in db_res]

def get_norm(run):
    (q, h_teff, p_teff, lt) =\
        db.retrieve('deuteron_db.db', 
                    'BCM4A_charge, HMS_TrkEff, SHMS_TrkEff, T6_tLT', 
                    'RUN_LIST', where = f"run=\'{run}\'")[0]
    norm = 1/(q*h_teff*p_teff*lt)
    return norm

def calc_weights(simc_branches):
    w = simc_branches['Weight']
    nf = simc_branches['Normfac']
    nevt = w.size
    
    return w*nf/nevt

def fit_legend(fit_params):
    ax = B.pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    for par in fit_params:
        fit_info = (
               f"A = {par['A'].value:.2f} $\pm$ {par['A'].err:.2f}\n"
               f"$\mu$ = {par['mean'].value:.4e} $\pm$ {par['mean'].err:.4e}\n"
               f"$\sigma$ = {par['sigma'].value:.4e} $\pm$ {par['sigma'].err:.4e}\n"
           )
        handles.append(mlines.Line2D([], [], color='None', label=fit_info))
        labels.append(fit_info)
    
    #ax.legend(handles=handles, labels=labels)
    return (handles, labels)
    #ax.legend(handles=handles, labels=labels, loc='best')
 
def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    #B.pl.legend(fontsize = fsize)
        
    
#%% Declare constants
Eb = 10.542         # GeV
Mp = 0.938272       # GeV

#%% Load hcana variables

# kin_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'kin\''))
# gtr_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'gtr\''))

my_branches = ['P.kin.primary.scat_ang_rad','H.kin.secondary.xangle','P.gtr.p',
               'H.gtr.p','H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th',
               'P.gtr.ph','H.kin.secondary.Prec_x','H.kin.secondary.Prec_y',
               'H.kin.secondary.Prec_z']

my_branches_sim = ['theta_e', 'theta_p', 'e_pf', 'h_pf', 'e_delta', 'h_delta',
                   'h_xptar','h_yptar','e_xptar','e_yptar','Weight','Normfac',
                   'Pmx','Pmy','Pmz']

branches_sel = my_branches
#%% load data root files
# make sure files are where they should be!
R = 20851
sett = 'delta_scan_0'
typ = '-'

RUN = D.DATA_INIT(data_type='deut23_data', run = R, 
                  select_branches={'T':branches_sel})
# load simc root files
# make sure files are where they should be!
SIMC = D.DATA_INIT(data_type='SIMC',
                   setting=sett,select_branches={'SNT':my_branches_sim},
                   simc_type=typ)     

#%% single run - assign measured quatities
# electron scattering angle
th_e = RUN.Branches['P.kin.primary.scat_ang_rad']
th_e_sim = SIMC.Branches['theta_e']

# angle of detected particle (proton) with scattered electron
xangle = RUN.Branches['H.kin.secondary.xangle']

# detected proton scattering angle
th_p_meas = xangle - th_e
th_p_meas_sim = SIMC.Branches['theta_p']

# measured electron final momentum
kf_meas = RUN.Branches['P.gtr.p']
kf_meas_sim = SIMC.Branches['e_pf']*(1/1e3)

# measured proton final momentum
Pf_meas = RUN.Branches['H.gtr.p']
Pf_meas_sim = SIMC.Branches['h_pf']*(1/1e3)

#%% single run - Declare calculated quatities

# Proton final momentum in GeV
Pf_calc = (2.*Mp*Eb*(Mp+Eb)*np.cos(th_p_meas))/\
            (Mp*Mp + 2.*Mp*Eb + Eb*Eb*np.sin(th_p_meas)*np.sin(th_p_meas))

Pf_calc_sim = (2.*Mp*Eb*(Mp+Eb)*np.cos(th_p_meas_sim))/\
    (Mp*Mp + 2.*Mp*Eb + Eb*Eb*np.sin(th_p_meas_sim)*np.sin(th_p_meas_sim))

# electron final momentum components in GeV
kf_x = Pf_calc*np.sin(th_p_meas)
kf_y = 0.
kf_z = Eb - Pf_calc*np.cos(th_p_meas)

kf_x_sim = Pf_calc_sim*np.sin(th_p_meas_sim)
kf_y_sim = 0.
kf_z_sim = Eb - Pf_calc_sim*np.cos(th_p_meas_sim)

# electron final momentum in GeV
kf_calc = np.sqrt(kf_x*kf_x + kf_y*kf_y + kf_z*kf_z)

kf_calc_sim = np.sqrt(kf_x_sim*kf_x_sim + kf_y_sim*kf_y_sim + kf_z_sim*kf_z_sim)

# electron recoil angle in rad
np.seterr(invalid='ignore') # ignore RunTimeWarning due to 0/0 division
        
th_e_calc = np.arctan(kf_x/kf_z)

th_e_calc_sim = np.arctan(kf_x_sim/kf_z_sim)

#%% single run - Calculated - Measured quatities   
dkf = kf_calc - kf_meas
dPf = Pf_calc - Pf_meas
dth_e = th_e_calc - th_e

dkf_sim = kf_calc_sim - kf_meas_sim
dPf_sim = Pf_calc_sim - Pf_meas_sim
dth_e_sim = th_e_calc_sim - th_e_sim

#%% single run - calculate normalization factors and weights
NORM = get_norm(R)
WEIGHTS = calc_weights(SIMC.Branches)

#%% single run - make histos of differences (no cuts) then plot 

dkf_hist = NORM*B.histo(dkf,range=(-0.05,0.05),bins=100)
dPf_hist = NORM*B.histo(dPf,range=(-0.04,0.04),bins=100)
dth_e_hist = NORM*B.histo(dth_e,range=(-0.005,0.005),bins=100)

# SIMC
dkf_hist_sim = B.histo(dkf_sim,range=(-0.05,0.05),bins=100,
                       weights=WEIGHTS,calc_w2=True)
dPf_hist_sim = B.histo(dPf_sim,range=(-0.04,0.04),bins=100,
                       weights=WEIGHTS,calc_w2=True)
dth_e_hist_sim = B.histo(dth_e_sim,range=(-0.005,0.005),bins=100,
                         weights=WEIGHTS,calc_w2=True)

hist_to_plot = {'$\Delta k_f$':dkf_hist,'$\Delta P_f$':dPf_hist,
                '$\Delta \\theta_e$':dth_e_hist}
sim_hist_to_plot = {'$\Delta k_f$':dkf_hist_sim,'$\Delta P_f$':dPf_hist_sim,
                    '$\Delta \\theta_e$':dth_e_hist_sim}

for plot in hist_to_plot:
    h = hist_to_plot[plot]
    hS = sim_hist_to_plot[plot]
    
    B.pl.figure()
    h.plot()
    hS.plot(filled=False,color='#ed6792')
    B.pl.title(f'{plot} {RUN.many}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    ax = B.pl.gca()
    ax.set_autoscaley_on(True)

#%% single run - Define cuts for data: MISSING coin time cut
# =============================================================================
# =============================================================================

# Note for the future, calling the cut class will change NaNs to 1e38
#   in the cut array directly, is there pass by reference/value in python?

cuts_list = C.acceptance_cuts

cuts_to_apply = []
for cut in cuts_list:
    br = RUN.Branches[C.HCANA_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply.append(cut_array)

all_cuts = cuts_to_apply[0]
for arr in cuts_to_apply:    
    all_cuts = all_cuts & arr    

#%% old way of doing cuts  
std_cuts = {'H.gtr.dp':(-10.,10.),'P.gtr.dp':(-10.,22.),
            'H.gtr.th':(-0.08,0.08),'H.gtr.ph':(-0.03,0.03),
            'P.gtr.th':(-0.03,0.03),'P.gtr.ph':(-0.03,0.03)}
     
cuts = []
flag = True
for vc in std_cuts:
    cut_min = RUN.Branches[vc] > std_cuts[vc][0] 
    cut_max = RUN.Branches[vc] < std_cuts[vc][1]
    this_cut = cut_min & cut_max
    this_cut = np.array(this_cut)
    
    if flag:
        cuts = this_cut
        flag = False
    else:
        cuts = cuts & this_cut   

#%% single run - Define cuts for simc         
# =============================================================================
# =============================================================================
cuts_to_apply_sim = []
for cut in cuts_list:
    cut.init()
    br = SIMC.Branches[C.SIMC_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply_sim.append(cut_array)

all_cuts_sim = cuts_to_apply_sim[0]
for arr in cuts_to_apply_sim:    
    all_cuts_sim = all_cuts_sim & arr 

#%% old way of doing cuts 
std_cuts_sim = {'h_delta':(-10.,10.),'e_delta':(-10.,22.),
            'h_xptar':(-0.08,0.08),'h_yptar':(-0.03,0.03),
            'e_xptar':(-0.03,0.03),'e_yptar':(-0.03,0.03)}
cuts_sim = []
flag = True
for vc in std_cuts_sim:
    cut_min = SIMC.Branches[vc] > std_cuts_sim[vc][0] 
    cut_max = SIMC.Branches[vc] < std_cuts_sim[vc][1]
    this_cut = cut_min & cut_max
    this_cut = np.array(this_cut)
    
    if flag:
        cuts_sim = this_cut
        flag = False
    else:
        cuts_sim = cuts_sim & this_cut
        
#%% single run - apply cuts
dkf_wcuts = dkf[all_cuts]
dPf_wcuts = dPf[all_cuts]
dth_e_wcuts = dth_e[all_cuts]

dkf_sim_wcuts = dkf_sim[all_cuts_sim]
dPf_sim_wcuts = dPf_sim[all_cuts_sim]
dth_e_sim_wcuts = dth_e_sim[all_cuts_sim]
#%%
W = WEIGHTS[all_cuts_sim]

#%% single run - delta theta_e; making the histos
dkf_hist_wcuts = NORM*B.histo(dkf_wcuts,range=(-0.05,0.05),bins=100)

dPf_hist_wcuts = NORM*B.histo(dPf_wcuts,range=(-0.04,0.04),bins=100)

dth_e_hist_wcuts = NORM*B.histo(dth_e_wcuts,range=(-0.005,0.005),bins=100)

# SIMC
dkf_hist_sim_wcuts = B.histo(dkf_sim_wcuts,range=(-0.05,0.05),bins=100,weights=W,calc_w2=True)

dPf_hist_sim_wcuts = B.histo(dPf_sim_wcuts,range=(-0.04,0.04),bins=100,weights=W,calc_w2=True)

dth_e_hist_sim_wcuts = B.histo(dth_e_sim_wcuts,range=(-0.005,0.005),bins=100,weights=W,calc_w2=True)

hist_to_plot = {'$\Delta k_f$':dkf_hist_wcuts,
                '$\Delta P_f$':dPf_hist_wcuts,
                '$\Delta \\theta_e$':dth_e_hist_wcuts}
sim_hist_to_plot = {'$\Delta k_f$':dkf_hist_sim_wcuts,
                    '$\Delta P_f$':dPf_hist_sim_wcuts,
                    '$\Delta \\theta_e$':dth_e_hist_sim_wcuts}

#%% single run - plotting and fitting the histos

for plot in hist_to_plot:

    B.pl.figure()
    s = RUN.setting
    h = hist_to_plot[plot]
    hsim = sim_hist_to_plot[plot]
    
    hsim.plot_exp(c='#ed6792',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='SIMC')
    hsim.fit(plot_fit=False)
    hsim.fit(hsim.mean.value-hsim.sigma.value,hsim.mean.value+hsim.sigma.value)
    hsim.plot_fit(color='#de425b')
    
    mu_s = hsim.mean.value
    
    h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='Data') 
    h.fit(plot_fit=False)
    h.fit(h.mean.value-h.sigma.value,h.mean.value+h.sigma.value)
    h.plot_fit(color='#1f77b4')
    
    mu = h.mean.value
    
    # if mu_s > mu:  
    #     diffs = mu_s - mu
    # else:
    #     diffs = mu - mu_s
    # diffs is SIMC-DATA
    diffs = mu_s - mu
    
    han, lab = fit_legend([h.fit_par,hsim.fit_par])
    
    #print(han,lab)
    
    diff_lab = f'SIMC - DATA = {diffs:.3e}'
    diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
    
    myhan = [han[1],han[2],han[0],han[3],diff_han]
    mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
    B.pl.legend(handles=myhan,labels=mylab,loc='best',fontsize=12)
    
    bot,top = plt.ylim()
    
    ### plot formatting
    # B.pl.title(f'$\Delta th_p$ for {s}')
    # B.pl.xlabel('')
    # B.pl.ylabel('')
    # B.pl.ylim(0,top)
    # #B.pl.legend(handles=han,labels=lab,loc='best')
    # B.pl.xlabel('[rad]')
    B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                colors=['#de425b','#1f77b4'],linestyles='--')
    B.pl.title(f'{plot} {RUN.many}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    ax = B.pl.gca()
    ax.set_autoscaley_on(True)    

#%% momentum component alignment

## in SIMC what is written is the missing momentum which is 
# NEGATIVE the recoil momentum: Pf + Pr = q -> Pr = q - Pf = -Pm
# Pm is also in the SIMC Lab coordinates (x down, y left), 
# in HCANA Pm is in the reference frame of q (q is z, x is qxe)
# so we need to compare the data with PmPer, PmPar, PmOop
pmx = -RUN.Branches['H.kin.secondary.Prec_x'][all_cuts]
pmy = -RUN.Branches['H.kin.secondary.Prec_y'][all_cuts]
pmz = RUN.Branches['H.kin.secondary.Prec_z'][all_cuts]

pmx_h = NORM*B.histo(pmx, range=(-0.1,0.1), bins=50)
pmy_h = NORM*B.histo(pmy, range=(-0.1,0.1), bins=50)
pmz_h = NORM*B.histo(pmz, range=(-0.025,0.2), bins=50)

spmx = SIMC.Branches['Pmx'][all_cuts_sim]
spmy = SIMC.Branches['Pmy'][all_cuts_sim]
spmz = SIMC.Branches['Pmz'][all_cuts_sim]

spmx_h = B.histo(spmx, range=(-0.5,0.5), bins=50,
                              weights=W, calc_w2=True)
spmy_h = B.histo(spmy, range=(-0.1,0.1), bins=50, 
                              weights=W, calc_w2=True)
spmz_h = B.histo(spmz, range=(-0.1,0.1), bins=50, 
                              weights=W, calc_w2=True)
 
histos_to_plot = {'Pmx':pmx_h, 'Pmy':pmy_h,
                  'Pmz':pmz_h}

simhistos_to_plot = {'Pmx':spmx_h, 'Pmy':spmy_h, 
                     'Pmz':spmz_h}

#%% plotting and fitting the histos

fit_lims = {'Pmx':(-0.5,0.5,-0.5,0.5),'Pmy':(-0.001,0.02,-0.001,0.02),
            'Pmz':(-0.05,0.05,-0.05,0.05)}
for plot in histos_to_plot:
    # fig = B.pl.figure(layout='constrained',figsize=(20,10))
    # fig.suptitle(f'{plot}', fontsize=16)
    B.pl.figure()
    s = RUN.setting
    h = histos_to_plot[plot]
    hsim = simhistos_to_plot[plot]
    
    h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='Data') 
    # h.fit(fit_lims[plot][2],fit_lims[plot][3],plot_fit=False)
    # h.fit(h.mean.value-h.sigma.value*0.5,h.mean.value+h.sigma.value*0.5 )
    # h.plot_fit(color='#1f77b4')
    
    # mu = h.mean.value

    hsim.plot_exp(c='#ed6792',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='SIMC')
    # hsim.fit(fit_lims[plot][0],fit_lims[plot][1],plot_fit=False)
    # hsim.fit(hsim.mean.value-hsim.sigma.value*0.5,hsim.mean.value+hsim.sigma.value*0.5 )
    # hsim.plot_fit(color='#de425b')
    
    # mu_s = hsim.mean.value

    # if mu_s > mu:  
    #     diffs = mu_s - mu
    # else:
    #     diffs = mu - mu_s
    
    #DATA-SIMC
    diffs = mu - mu_s
    
    #han, lab = fit_legend([h.fit_par,hsim.fit_par])
    
    if plot == 'Pmx':
        diffs_Pmx = diffs
    elif plot == 'Pmy':
        diffs_Pmy = diffs
    elif plot == 'Pmz':
        diffs_Pmz = diffs   
    
    # #print(han,lab)
    
    # diff_lab = f'DATA - SIMC = {diffs:.3e}'
    # diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
    
    # myhan = [han[1],han[2],han[0],han[3],diff_han]
    # mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
    # B.pl.legend(handles=myhan,labels=mylab,loc='best',fontsize='xx-small')
    
    bot,top = plt.ylim()

    # # plot formatting
    # B.pl.title(f'$\Delta th_p$ for {s}')
    # B.pl.xlabel('')
    # B.pl.ylabel('')
    # B.pl.ylim(0,top)
    # #B.pl.legend(handles=han,labels=lab,loc='best')
    # B.pl.xlabel('[rad]')
    # B.pl.vlines([hsim.mean.value,h.mean.value], 0,[hsim.A.value,h.A.value],
    #             colors=['#de425b','#1f77b4'],linestyles='--')

    ax = B.pl.gca()    
    ax.set_autoscaley_on(True)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')

    