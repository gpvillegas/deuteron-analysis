#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 10:19:47 2024

@author: gvill

Optics Optimization Procedure 

Based on 'centroid-alignment' procedure

Centroid alignment: check for angle/momentum offset in electron arm angle, 
by calculating proton arm angle/momentum from e arm quantities. Starting from
electron scattering angle.
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import LT.box as B
import cut_handler as C
import database_operations as db
import data_init as D
from matplotlib.transforms import Bbox

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
    
def full_extent(ax,padw=0.,padh=0.,tx=0.,ty=0.):
    """Get the full extent of an axes, including axes labels, tick labels, and
    titles."""
    # For text objects, we need to draw the figure first, otherwise the extents
    # are undefined.
    ax.figure.canvas.draw()
    items = ax.get_xticklabels() + ax.get_yticklabels()
    items += [ax, ax.title, ax.xaxis.label, ax.yaxis.label]
    # items += [ax, ax.title]
    bbox = Bbox.union([item.get_window_extent() for item in items])
    bbox1 = bbox.expanded(1.0 + padw, 1.0 + padh)
    return bbox1.translated(tx,ty) 

        
    
#%% Declare constants
Eb = 10.542         # GeV
Mp = 0.938272       # GeV

#%% Load hcana variables

# kin_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'kin\''))
# gtr_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'gtr\''))
# branches_sel = kin_branches + gtr_branches

kin_var = ['P.kin.primary.scat_ang_rad','H.kin.secondary.xangle','P.gtr.p',
           'H.gtr.p','P.kin.primary.x_bj','P.kin.primary.W',
           'H.kin.secondary.emiss','H.kin.secondary.pmiss_x',
           'H.kin.secondary.pmiss_y','H.kin.secondary.pmiss_z']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

br_sel = kin_var + acc_var + calPID_var

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','theta_e','theta_p',
               'e_pf','h_pf','e_xptar','e_yptar','h_xptar','h_yptar','W','Em',
               'Pmx','Pmy','Pmz']
#%% load data root files
# make sure files are where they should be!
run_num = 20851
sett = 'delta_scan_0'

RUN = D.DATA_INIT(data_type='deut23_data', run = run_num, 
                  select_branches={'T':br_sel})

SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                    setting=sett,select_branches={'SNT':br_sel_SIMC},simc_type='-')

#%% Assign measured quatities

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

#%% Declare calculated quatities

# electron final momentum in GeV
kf_calc = (Mp*Eb)/(Mp + 2.*Eb*np.sin(0.5*th_e)*np.sin(0.5*th_e))

kf_calc_sim = (Mp*Eb)/(Mp + 2.*Eb*np.sin(0.5*th_e_sim)*np.sin(0.5*th_e_sim))

# proton final momentum components in GeV
Pf_x = kf_calc*np.sin(th_e)
Pf_y = 0.
Pf_z = Eb - kf_calc*np.cos(th_e)

Pf_x_sim = kf_calc_sim*np.sin(th_e_sim)
Pf_y_sim = 0.
Pf_z_sim = Eb - kf_calc_sim*np.cos(th_e_sim)

# proton final momentum in GeV
Pf_calc = np.sqrt(Pf_x*Pf_x + Pf_y*Pf_y + Pf_z*Pf_z)

Pf_calc_sim = np.sqrt(Pf_x_sim*Pf_x_sim + Pf_y_sim*Pf_y_sim + Pf_z_sim*Pf_z_sim)

# proton recoil angle
np.seterr(invalid='ignore') # ignore RunTimeWarning due to 0/0 division
        
th_p_calc = np.arctan(Pf_x/Pf_z)

th_p_calc_sim = np.arctan(Pf_x_sim/Pf_z_sim)

#%% Calculated - Measured quatities

dkf = kf_calc - kf_meas
dPf = Pf_calc - Pf_meas
dth_p = th_p_calc - th_p_meas

dkf_sim = kf_calc_sim - kf_meas_sim
dPf_sim = Pf_calc_sim - Pf_meas_sim
dth_p_sim = th_p_calc_sim - th_p_meas_sim

#%%
# Plot differences (no cuts)
dkf_hist = B.histo(dkf,range=(-0.1,0.1),bins=100)
dkf_hist.plot()

B.pl.figure()
dPf_hist = B.histo(dPf,range=(-1.,1.),bins=200)
dPf_hist.plot(filled=False)

B.pl.figure()
dth_p_hist = B.histo(dth_p,range=(-1.,1.),bins=200)
dth_p_hist.plot(filled=False)

# SIMC
WEIGHTS = calc_weights(SIMC.Branches)
dkf_hist_sim = B.histo(dkf_sim,range=(-0.1,0.1),bins=100,weights=WEIGHTS)
dkf_hist_sim.plot()

B.pl.figure()
dPf_hist = B.histo(dPf,range=(-1.,1.),bins=200)
dPf_hist.plot(filled=False)

B.pl.figure()
dth_p_hist_sim = B.histo(dth_p_sim,range=(-1.,1.),bins=200,weights=WEIGHTS)
dth_p_hist_sim.plot(filled=False)


#%% define cuts
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

#cuts for SIMC     
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
#%% get normalization factor and weights
NORM = get_norm(run_num)
WEIGHTS = calc_weights(SIMC.Branches)
WEIGHTS = WEIGHTS[all_cuts_sim]

#%%
#plot kins

W = RUN.Branches['P.kin.primary.W'][all_cuts]
Em = RUN.Branches['H.kin.secondary.emiss'][all_cuts]

W_sim = SIMC.Branches['W'][all_cuts_sim]
Em_sim = SIMC.Branches['Em'][all_cuts_sim]

W_h = NORM*B.histo(W,range=(0.85,1.05),bins=100)
Em_h = NORM*B.histo(Em, range=(-0.1,0.1), bins=100)

W_sim_h = B.histo(W_sim,range=(0.85,1.05),bins=100, weights=WEIGHTS, calc_w2=True)
Em_sim_h = B.histo(Em_sim, range=(-0.1,0.1), bins=100, weights=WEIGHTS, calc_w2=True)

B.pl.figure()
W_h.plot(color='#7b8fd4')
W_h.fit(0.925,0.975)
W_h.plot_fit(color='#1f77b4')

W_sim_h.plot(filled=False, color='#ed6792')
W_sim_h.fit(0.91,0.975)
W_sim_h.plot_fit(color='#de425b')

dms_W = W_h.mean.value - W_sim_h.mean.value
ax = B.pl.gca()
ax.set_autoscaley_on(True)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'W Run {run_num} {sett}')

handles, labels = ax.get_legend_handles_labels()

diff_lab = f'DATA - SIMC = {dms_W:.3e}'
diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
handles.append(diff_han)
labels.append(diff_lab)

B.pl.legend(handles=handles,labels=labels,loc='best',fontsize='medium')

B.pl.figure()
Em_h.plot(color='#7b8fd4')
Em_h.fit(0,0.01)
Em_h.plot_fit(color='#1f77b4')

Em_sim_h.plot(filled=False, color='#ed6792')
Em_sim_h.fit(0,0.01)
Em_sim_h.plot_fit(color='#de425b')

dms_Em = Em_h.mean.value - Em_sim_h.mean.value

ax = B.pl.gca()
ax.set_autoscaley_on(True)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title(f'Em Run {run_num} {sett}')

handles, labels = ax.get_legend_handles_labels()

diff_lab = f'DATA - SIMC = {dms_Em:.3e}'
diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
handles.append(diff_han)
labels.append(diff_lab)

B.pl.legend(handles=handles,labels=labels,loc='best',fontsize='medium')


#%% Plot differences with cuts; apply cuts

dkf = dkf[all_cuts]
dPf = dPf[all_cuts]
dth_p = dth_p[all_cuts]

dkf_sim = dkf_sim[all_cuts_sim]
dPf_sim = dPf_sim[all_cuts_sim]
dth_p_sim = dth_p_sim[all_cuts_sim]

#%% delta theta_p; making the histos

dthp_wcuts_histos = NORM*B.histo(dth_p,(-0.015,0.015),100)
dkf_wcuts_histos = NORM*B.histo(dkf,(-0.1,0.1),100)
dPf_wcuts_histos = NORM*B.histo(dPf,(-0.1,0.1),100)

dthp_wcuts_simhistos = B.histo(dth_p_sim,(-0.015,0.015),100,weights=WEIGHTS,calc_w2=True)
dkf_wcuts_simhistos = B.histo(dkf_sim,(-0.1,0.1),100,weights=WEIGHTS,calc_w2=True)
dPf_wcuts_simhistos = B.histo(dPf_sim,(-0.1,0.1),100,weights=WEIGHTS,calc_w2=True)
#%% plotting and fitting the histos
hist_to_plot = {'$\Delta \\theta_p$':dthp_wcuts_histos,
                '$\Delta k_f$':dkf_wcuts_histos,
                '$\Delta P_f$':dPf_wcuts_histos}
sim_hist_to_plot = {'$\Delta \\theta_p$':dthp_wcuts_simhistos,
                '$\Delta k_f$':dkf_wcuts_simhistos,
                '$\Delta P_f$':dPf_wcuts_simhistos}

for plot in hist_to_plot:
    fig = B.pl.figure(layout='constrained')
    fig.suptitle(f'{plot}', fontsize=16)

    s = RUN.setting
    h = hist_to_plot[plot]
    hsim = sim_hist_to_plot[plot]

    hsim.plot_exp(c='#ed6792',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='SIMC')
    hsim.fit(plot_fit=False)
    hsim.fit(hsim.mean.value-hsim.sigma.value*0.5,hsim.mean.value+hsim.sigma.value*0.5 )
    hsim.plot_fit(color='#de425b')
    
    mu_s = hsim.mean.value
    
    h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                         capsize=0,mew=1,elinewidth=1,label='Data') 
    h.fit(plot_fit=False)
    h.fit(h.mean.value-h.sigma.value*0.5,h.mean.value+h.sigma.value*0.5 )
    h.plot_fit(color='#1f77b4')
    
    mu = h.mean.value

    if mu_s > mu:  
        diffs = mu_s - mu
    else:
        diffs = mu - mu_s
    
    han, lab = fit_legend([h.fit_par,hsim.fit_par])
    
    if plot == '$\Delta \\theta_p$':
        diffs_thp = diffs
    elif plot == '$\Delta k_f$':
        diffs_kf = diffs
    elif plot == '$\Delta P_f$':
        diffs_Pf = diffs   

    # #print(han,lab)
    
    diff_lab = f'DATA - SIMC = {diffs:.3e}'
    diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
    
    myhan = [han[1],han[2],han[0],han[3],diff_han]
    mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
    B.pl.legend(handles=myhan,labels=mylab,loc='best',fontsize='medium')
    
    bot,top = plt.ylim()

    # # plot formatting
    # B.pl.title(f'$\Delta th_p$ for {s}')
    # B.pl.xlabel('')
    # B.pl.ylabel('')
    # B.pl.ylim(0,top)
    # #B.pl.legend(handles=han,labels=lab,loc='best')
    # B.pl.xlabel('[rad]')
    B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                colors=['#de425b','#1f77b4'],linestyles='--')
    
    ax = B.pl.gca()
    ax.set_autoscaley_on(True)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'Run {run_num} {sett}')

    # pptxify(t = f'Run: {r} Setting: {s}',
    #         x = '$\Delta \\theta_p = \\theta_{p,calc} - \\theta_{p,meas}$ [rad]',
    #         y = '')
    
    # save_where = './offset_studies/step1/'
    # fname = save_where + f'{plot}_{r}_{s}_1.png'
    # extent = full_extent(ax).transformed(fig.dpi_scale_trans.inverted())
    # B.pl.savefig(fname,bbox_inches=extent)

    
#%% DATA - SIMC differences

delta = np.array([-8,-4,0,4,8,12])
diffs_sorted = np.array([diffs_kf[20846],diffs_kf[20851],diffs_kf[20858],
                diffs_kf[20861],diffs_kf[20840],diffs_kf[20869]])


plt.scatter(delta, diffs_sorted)
plt.title('DATA-SIMC for each $\delta$ setting (STEP 1)')
B.pl.savefig('STEP1_DATA-SIMC_vs_delta.png') 

#%% momentum component alignment

Pmx = RUN.Branches['H.kin.secondary.pmiss_x']
Pmy = RUN.Branches['H.kin.secondary.pmiss_y']
Pmz = RUN.Branches['H.kin.secondary.pmiss_z']

Pmx_sim = SIMC.Branches['Pmx']
Pmy_sim = SIMC.Branches['Pmy']
Pmz_sim = SIMC.Branches['Pmz']

# apply cuts
Pmx = Pmx[all_cuts]
Pmy = Pmy[all_cuts]
Pmz = Pmz[all_cuts]

Pmx_sim = Pmx_sim[all_cuts_sim]
Pmy_sim = Pmy_sim[all_cuts_sim]
Pmz_sim = Pmz_sim[all_cuts_sim]

# make histograms

Pmx_wcuts_histo = NORM*B.histo(Pmx, range=(-0.025,0.050), bins=100)
Pmy_wcuts_histo = NORM*B.histo(Pmy, range=(-0.025,0.050), bins=100)
Pmz_wcuts_histo = NORM*B.histo(Pmz, range=(-0.025,0.050), bins=100)

Pmx_wcuts_histo_sim = B.histo(Pmx_sim, range=(-0.025,0.050), bins=100,
                              weights=WEIGHTS, calc_w2=True)
Pmy_wcuts_histo_sim = B.histo(Pmy_sim, range=(-0.025,0.050), bins=100, 
                              weights=WEIGHTS, calc_w2=True)
Pmz_wcuts_histo_sim = B.histo(Pmz_sim, range=(-0.025,0.050), bins=100, 
                              weights=WEIGHTS, calc_w2=True)

histos_to_plot = {'Pmx':Pmx_wcuts_histo, 'Pmy':Pmy_wcuts_histo, 'Pmz':Pmz_wcuts_histo}

simhistos_to_plot = {'Pmx':Pmx_wcuts_histo_sim, 'Pmy':Pmy_wcuts_histo_sim, 'Pmz':Pmz_wcuts_histo_sim}

# plot histos and fit

for plot in histos_to_plot:

    B.pl.figure()
    s = RUN.setting
    h = histos_to_plot[plot]
    hsim = simhistos_to_plot[plot]
    
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
    B.pl.title(f'{plot}')
    B.pl.xlabel('')
    B.pl.ylabel('')
    ax = B.pl.gca()
    ax.set_autoscaley_on(True)


