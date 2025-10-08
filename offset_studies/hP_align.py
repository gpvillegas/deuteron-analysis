#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:54:02 2025

@author: gvill

HMS check momentum alignment: 
        This code checks for central momentum alignment
        and plots HMS delta vs focal plane variables to 
        check for correlations

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import colormaps as cmp
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

def fit_legend_line(fit):
    ax = B.pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    fit_info = (
           f"slope = {fit.slope:.3e} $\pm$ {fit.sigma_s:.3e}\n"
           f"offset = {fit.offset:.3e} $\pm$ {fit.sigma_o:.3e}"
       )
    handles.append(mlines.Line2D([], [], color='None', label=fit_info))
    labels.append(fit_info)
    
    ax.legend(handles=handles, labels=labels)
    return (handles, labels)
 
def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    #B.pl.legend(fontsize = fsize)
    
def combine_histos(histos,many=[]):
    loop1 = True 
    for m in many:
        if loop1:
            hsum = histos[m]
            loop1 = False
        else:
            hsum += histos[m]
    
    return hsum    

def set_nan(inpt,set_nan_to = 0):
    has_nan = np.isnan(inpt)
    inpt[has_nan] = set_nan_to
    print(f'Set {np.where(has_nan==True)[0].size} nans to {set_nan_to}\n')
    
    return

def set_zero(inpt,set_zero_to=1):
    is_zero = np.where(inpt == 0.)
    for i in is_zero[0]:
        inpt[i] = set_zero_to
    print(f'Set {is_zero[0].size} zeroes to {set_zero_to}\n')
    return        
    
#%% Declare constants
Eb = 10.542         # GeV
Mp = 0.938272       # GeV

dtr = np.pi/180.

#%% Load hcana variables

# kin_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'kin\''))
# gtr_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'gtr\''))
sieve_branches = ['H.extcor.xsieve', 'H.extcor.ysieve',
                  'P.extcor.xsieve', 'P.extcor.ysieve']

my_branches = ['P.kin.primary.scat_ang_rad','H.kin.secondary.xangle','P.gtr.p',
               'H.gtr.p','H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th',
               'P.gtr.ph','P.kin.primary.x_bj','H.kin.secondary.pmiss_x',
               'H.kin.secondary.pmiss_y','H.kin.secondary.pmiss_z',
               'H.kin.secondary.pmiss','H.kin.secondary.Prec_x',
               'H.kin.secondary.Prec_y','H.kin.secondary.Prec_z',
               'P.kin.primary.W']

my_branches_sim = ['theta_e', 'theta_p', 'e_pf', 'h_pf', 'e_delta', 'h_delta',
                   'h_xptar','h_yptar','e_xptar','e_yptar','Weight','Normfac',
                   'Pmx','Pmy','Pmz','PmPer','PmPar','PmOop','W']

branches_sel = my_branches + sieve_branches
#%% load data root files
# make sure files are where they should be!
# RUN = D.DATA_INIT(data_type='deut23_data', run = 20851, 
#                   select_branches=branches_sel)
# R = [20871,20873,20888,20958]

# RUN = D.DATA_INIT(data_type='deut23_data', run=R, 
#                     select_branches={'T':branches_sel})

RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                  select_branches={'T':branches_sel})

# RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_singles', 
#                   select_branches=branches_sel)
#%% load simc root files
# make sure files are where they should be!
# SIMC = D.DATA_INIT(data_type='SIMC', kin_study='deep',
#                    simc_type = 'jmlfsi_rad',select_branches={'SNT':my_branches_sim})


SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                   select_branches={'SNT':my_branches_sim},simc_type='-')

# SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
#                    setting='delta_scan_0')     
    
# SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_singles') 
#%% Assign measured quatities
#testing offsets
Pf_offset = 1.0

th_e = {}
xangle = {}
th_p_meas = {}
kf_meas = {}
Pf_meas = {}
for r in RUN.many:
    # electron scattering angle in rad
    th_e[r] = RUN.Branches[r]['P.kin.primary.scat_ang_rad']

    # angle of detected particle (proton) with scattered electron in rad
    xangle[r] = RUN.Branches[r]['H.kin.secondary.xangle']

    # detected proton scattering angle in rad
    th_p_meas[r] = xangle[r] - th_e[r]

    # measured electron final momentum in GeV
    kf_meas[r] = RUN.Branches[r]['P.gtr.p']

    # measured proton final momentum in GeV
    Pf_meas[r] = RUN.Branches[r]['H.gtr.p']*Pf_offset

th_e_sim = {}
th_p_meas_sim = {}
kf_meas_sim = {}
Pf_meas_sim = {}
for s in SIMC.many:
    th_e_sim[s] = SIMC.Branches[s]['theta_e']

    th_p_meas_sim[s] = SIMC.Branches[s]['theta_p']

    kf_meas_sim[s] = SIMC.Branches[s]['e_pf']*(1/1e3)

    Pf_meas_sim[s] = SIMC.Branches[s]['h_pf']*(1/1e3)

#%% Declare calculated quatities

Pf_calc = {}
kf_calc = {}
th_e_calc = {}
dPdthp = {}

for r in RUN.many:
    # proton angle in rad
    thp = th_p_meas[r]
    
    # proton final momentum in GeV
    pf = (2.*Mp*Eb*(Mp+Eb)*np.cos(thp))/\
                (Mp*Mp + 2.*Mp*Eb + Eb*Eb*np.sin(thp)*np.sin(thp))

    # electron final momentum components in GeV
    kfx = pf*np.sin(thp)
    kfy = 0.
    kfz = Eb - pf*np.cos(thp)

    # electron final momentum in GeV
    kf = np.sqrt(kfx*kfx + kfy*kfy + kfz*kfz)

    # electron recoil angle in rad
    np.seterr(invalid='ignore') # ignore RunTimeWarning due to 0/0 division
            
    the = np.arctan(kfx/kfz)
    
    Pf_calc[r] = pf
    kf_calc[r] = kf
    th_e_calc[r] = the

Pf_calc_sim = {}
kf_calc_sim = {}
th_e_calc_sim = {}

for s in SIMC.many:
    # proton angle in rad
    thp = th_p_meas_sim[s]
    
    # proton final momentum in GeV
    pf = (2.*Mp*Eb*(Mp+Eb)*np.cos(thp))/\
                (Mp*Mp + 2.*Mp*Eb + Eb*Eb*np.sin(thp)*np.sin(thp))

    # electron final momentum components in GeV
    kfx = pf*np.sin(thp)
    kfy = 0.
    kfz = Eb - pf*np.cos(thp)

    # electron final momentum in GeV
    kf = np.sqrt(kfx*kfx + kfy*kfy + kfz*kfz)

    # electron recoil angle in rad
    np.seterr(invalid='ignore') # ignore RunTimeWarning due to 0/0 division
            
    the = np.arctan(kfx/kfz)
    
    Pf_calc_sim[s] = pf
    kf_calc_sim[s] = kf
    th_e_calc_sim[s] = the

#%% Calculated - Measured quatities

dkf = {}
dPf = {}
dth_e = {}
for r in RUN.many:
    dkf[r] = kf_calc[r] - kf_meas[r]
    dPf[r] = Pf_calc[r] - Pf_meas[r]
    dth_e[r] = th_e_calc[r] - th_e[r]

dkf_sim = {}
dPf_sim = {}
dth_e_sim = {}
for s in SIMC.many:
    dkf_sim[s] = kf_calc_sim[s] - kf_meas_sim[s]
    dPf_sim[s] = Pf_calc_sim[s] - Pf_meas_sim[s]
    dth_e_sim[s] = th_e_calc_sim[s] - th_e_sim[s]
    
#%% Define cuts for data
cuts_list = C.acceptance_cuts

cuts_to_apply = {}
for r in RUN.many:
    c_list = []
    print(f'Cuts for Run {r}')
    for cut in cuts_list:
        cut.init()
        br = RUN.Branches[r][C.HCANA_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
    
    # add collimator cut
    hxc = RUN.Branches[r]['H.extcor.xsieve']
    hyc = RUN.Branches[r]['H.extcor.ysieve']
    
    hcoll_cut = C.coll_cut(hxc, hyc, spec='HMS')
    c_list.append(hcoll_cut) 
    
    cuts_to_apply[r] = c_list
    
all_cuts = {}
for r in RUN.many:
    all_cuts_arr = cuts_to_apply[r][0]
    for arr in cuts_to_apply[r]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts[r] = all_cuts_arr

#%% Define cuts for SIMC

cuts_to_apply_sim = {}
for s in SIMC.many:
    c_list = []
    print(f'Cuts for SIMC setting {s}')
    for cut in cuts_list:
        cut.init()
        br = SIMC.Branches[s][C.SIMC_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
        
    cuts_to_apply_sim[s] = c_list

all_cuts_sim = {}
for s in SIMC.many:
    all_cuts_arr = cuts_to_apply_sim[s][0]
    for arr in cuts_to_apply_sim[s]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts_sim[s] = all_cuts_arr
    
#%% apply cuts
for r in RUN.many:
    dkf[r] = dkf[r][all_cuts[r]]
    dPf[r] = dPf[r][all_cuts[r]]
    dth_e[r] = dth_e[r][all_cuts[r]]
   
for s in SIMC.many:
    dkf_sim[s] = dkf_sim[s][all_cuts_sim[s]]
    dPf_sim[s] = dPf_sim[s][all_cuts_sim[s]]
    dth_e_sim[s] = dth_e_sim[s][all_cuts_sim[s]]
    
#%% calculate normalization factor for data and locate weigths for SIMC
NORM = {}
for r in RUN.many:
    NORM[r] = get_norm(r)

WEIGHTS = {}    
for s in SIMC.many:    
    W = calc_weights(SIMC.Branches[s])
    WEIGHTS[s] = W[all_cuts_sim[s]]
#%% make histos
dkf_hist_wcuts = {}
dPf_hist_wcuts = {}
dth_e_hist_wcuts = {}

for r in RUN.many:    
    dkf_hist_wcuts[r] = NORM[r]*B.histo(dkf[r],range=(-0.05,0.05),
                                        bins=100)

    dPf_hist_wcuts[r] = NORM[r]*B.histo(dPf[r],range=(-0.04,0.04),
                                        bins=100)

    dth_e_hist_wcuts[r] = NORM[r]*B.histo(dth_e[r],range=(-0.005,0.005),
                                          bins=100)
     
dkf_hist_sim_wcuts = {}
dPf_hist_sim_wcuts = {}
dth_e_hist_sim_wcuts = {}

for s in SIMC.many:
    dkf_hist_sim_wcuts[s] = B.histo(dkf_sim[s],range=(-0.05,0.05),
                                    bins=100,weights=WEIGHTS[s],calc_w2=True)
    
    dPf_hist_sim_wcuts[s] = B.histo(dPf_sim[s],range=(-0.04,0.04),
                                    bins=100,weights=WEIGHTS[s],calc_w2=True)
    
    dth_e_hist_sim_wcuts[s] = B.histo(dth_e_sim[s],range=(-0.005,0.005),
                                      bins=100,weights=WEIGHTS[s],calc_w2=True)

hist_to_plot = {'$\Delta k_f$':dkf_hist_wcuts,
                '$\Delta P_f$':dPf_hist_wcuts,
                '$\Delta \\theta_e$':dth_e_hist_wcuts}
sim_hist_to_plot = {'$\Delta k_f$':dkf_hist_sim_wcuts,
                    '$\Delta P_f$':dPf_hist_sim_wcuts,
                    '$\Delta \\theta_e$':dth_e_hist_sim_wcuts}    
 
#%% plotting and fitting the histos
offsets = {}
diffs_the = {}
diffs_kf = {}
diffs_Pf = {}
for plot in hist_to_plot:
    fig = B.pl.figure(layout='constrained',figsize=(20,10))
    fig.suptitle(f'{plot}', fontsize=16)
    rows = 2
    cols = 4
    i = 1
    off = []
    for r in RUN.many:       
        s = RUN.setting[r]
        ax = plt.subplot(rows,cols,i)
        h = hist_to_plot[plot][r]
        hsim = sim_hist_to_plot[plot][s]
        
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
        
        diffs = mu_s - mu
        
        if plot == '$\Delta \\theta_e$':
            diffs_the[r] = diffs
        elif plot == '$\Delta k_f$':
            diffs_kf[r] = diffs
        elif plot == '$\Delta P_f$':
            diffs_Pf[r] = diffs 
        
        han, lab = fit_legend([h.fit_par,hsim.fit_par])
        
        #print(han,lab)
        
        diff_lab = f'SIMC - DATA = {diffs:.3e}'
        diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
        
        myhan = [han[1],han[2],han[0],han[3],diff_han]
        mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
        B.pl.legend(handles=myhan,labels=mylab,loc='best')
        
        bot,top = plt.ylim()
        
        B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                    colors=['#de425b','#1f77b4'],linestyles='--')
        
        i+=1
        ax.set_autoscaley_on(True)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'Run {r} {s}')  
        
        # save_where = '/home/gvill/Desktop/plots for pptx/'
        # fname = save_where + f'dthp_{r}_{s}_oldOffsets.png'
        # B.pl.savefig(fname) 
        
        off.append(diffs)
    offsets[plot] = off   
#%% get mean of diffs

pf_offset = np.array(offsets['$\\Delta P_f$'][1:5])
pf_offset = pf_offset.mean()        

print('pf_offset = ', pf_offset)

#%% calculate dpf offsets to thp
delta = np.array([-8,-4,0,4,8,12])
pf_cent = np.array([3.499,3.145,2.783,2.417,2.048,1.664])
deep_pf_cent = np.array([3.0523,2.2622,2.1210,2.0474])
thp_cent = np.array([33.34,35.75,38.54,41.795,45.645,50.5])

# derivative respect to thp
# t1 = 2*Mp*Eb*(Eb+Mp)*np.cos(thp_cent*dtr)
# dt1 = -2*Mp*Eb*(Eb+Mp)*np.sin(thp_cent*dtr)
# t2 = Mp**2 + 2*Mp*Eb + (Eb**2)*(np.sin(thp_cent*dtr)**2)
# dt2 = 2*(Eb**2)*np.sin(thp_cent*dtr)*np.cos(thp_cent*dtr)   
# dPdthp = (dt1/t2) - (t1*dt2)/t2**2

fac1 = 2*Mp*Eb*(Eb+Mp)*np.sin(thp_cent*dtr)
fac2 = Eb*Eb*np.sin(thp_cent*dtr)*np.sin(thp_cent*dtr) + Mp*Mp + 2*Mp*Eb
fac = -fac1/fac2
term = (2*Eb*Eb*np.cos(thp_cent*dtr)*np.cos(thp_cent*dtr))
dPdthp = fac*(1+term/fac2)

dpf = np.array([diffs_Pf[20841],diffs_Pf[20846],diffs_Pf[20851],
                diffs_Pf[20858],diffs_Pf[20861],diffs_Pf[20869]])

thp_offsets = dpf/dPdthp
print(thp_offsets/dtr)
B.pl.figure()
B.pl.plot(pf_cent,thp_offsets/dtr,'.',label='heep_coin')
line = B.linefit(pf_cent,thp_offsets/dtr,plot_fit=False)
line.plot(color='r')
B.pl.plot(deep_pf_cent,line(deep_pf_cent),'.',color='black', label ='deep')
print(line(deep_pf_cent))
B.pl.text(1.75,-0.060,f'slope = {line.slope:.4e}$\pm${line.sigma_s:.4e}')
B.pl.text(1.75,-0.063,f'offset = {line.offset:.4e}$\pm${line.sigma_o:.4e}')
B.pl.legend()
#%% plot pf_calc vs delta/ xbj vs delta

Pf_calc_cut = {}
Pf_meas_cut = {}
pf_nom = {20840:1.664,20841:3.499,20846:3.145,20851:2.783,20858:2.417,
          20861:2.048,20868:1.664,20869:1.664}
delta = {}
xbj = {}

for m in RUN.many:
    
    p = Pf_calc[m][all_cuts[m]]/pf_nom[m]
    Pf_calc_cut[m] = p
    d = RUN.Branches[m]['H.gtr.dp']
    delta[m] = d[all_cuts[m]]
    
    xbj[m] = RUN.Branches[m]['P.kin.primary.x_bj'][all_cuts[m]]  
                 
# make and store 2d histograms
Pf_calc_vs_delta_h = {}
xbj_vs_delta_h = {}
xbj_h = {}
xbj_peaks = {}
for m in RUN.many:
    pf = Pf_calc_cut[m]
    x = xbj[m]
    d = delta[m]

    h = B.histo2d(d,100*(pf-1)-d,range=[(-10.5,10.5),(-5,5)],bins=100)
    h2 = B.histo2d(d,x,range=[(-10.5,10.5),(0.75,1.25)],bins=100)
    h3 = B.histo(x,range=(0.9,1.1),bins=100)
    h3.fit(plot_fit=False)
    h3.fit(h3.mean.value-h3.sigma.value,h3.mean.value+h3.sigma.value,
           plot_fit=False)
    
    
    Pf_calc_vs_delta_h[m] = h
    xbj_vs_delta_h[m] = h2
    xbj_h[m] = h3
    xbj_peaks[m] = h3.mean.value
    
    if m == 20840:
        y_fit = 100*(pf-1)-d
        x_fit = d
    else:
        y_fit = np.append(y_fit,100*(pf-1)-d)
        x_fit = np.append(x_fit,d)
    
    # B.pl.figure()
    # h3.plot()
    # h3.plot_fit()
    
# combine all histos

Pf_calc_vs_delta_h_all = combine_histos(Pf_calc_vs_delta_h,RUN.many)
xbj_vs_delta_h_all = combine_histos(xbj_vs_delta_h,RUN.many) 

range_y = (y_fit >= -0.5) & (y_fit <= 0.4)
range_x = (x_fit >= -8) & (x_fit <= 8)
range_fit = range_x & range_y

fit = B.polyfit(x_fit[range_fit], y_fit[range_fit],order=1)

B.pl.figure()
Pf_calc_vs_delta_h_all.plot(colormap=cmp['viridis'])
fit.plot(color='r')
B.pl.title('((Pf_calc/P_cent)-1)-$\delta$ vs $\delta$')
B.pl.xlabel('[%]')
B.pl.ylabel('[%]')

#B.pl.figure()
#xbj_vs_delta_h_all.plot(colormap=cmp['viridis'])
# B.pl.title('x_bj vs $\delta$')
# B.pl.xlabel('[%]')
# B.pl.ylabel('[]')

                 
#%% DATA - SIMC differences

delta = np.array([-8,-4,0,4,8,12])
# diffs_sorted = np.array([diffs_kf[20841],diffs_kf[20846],diffs_kf[20851],
#                 diffs_kf[20858],diffs_kf[20861],diffs_kf[20869]])
# diffs_sorted = np.array([diffs_the[20841],diffs_the[20846],diffs_the[20851],
#                 diffs_the[20858],diffs_the[20861],diffs_the[20869]])
diffs_sorted = np.array([diffs_Pf[20841],diffs_Pf[20846],diffs_Pf[20851],
                diffs_Pf[20858],diffs_Pf[20861],diffs_Pf[20869]])
# diffs_sorted = np.array([xbj_peaks[20841],xbj_peaks[20846],xbj_peaks[20851],
#                 xbj_peaks[20858],xbj_peaks[20861],xbj_peaks[20869]])

B.pl.figure()
B.pl.scatter(delta, diffs_sorted)
line = B.linefit(delta, diffs_sorted)
    
plt.title('DATA-SIMC per delta setting')
han,lab = fit_legend_line(line)
B.pl.legend(handles=han,labels=lab)
# B.pl.savefig('STEP1_DATA-SIMC_vs_delta.png')

#%% momentum component alignment
Pmx = {}
Pmy = {}
Pmz = {}
Pmx_wcuts_histo = {}
Pmy_wcuts_histo = {}
Pmz_wcuts_histo = {}

for r in RUN.many:
    ## in SIMC what is written is the missing momentum which is 
    # NEGATIVE the recoil momentum: Pf + Pr = q -> Pr = q - Pf = -Pm
    # Pm is also in the SIMC Lab coordinates (x down, y left), 
    # in HCANA Pm is in the reference frame of q (q is z, x is qxe)
    # so we need to compare the data with PmPer, PmPar, PmOop
    pmx = -RUN.Branches[r]['H.kin.secondary.Prec_x'][all_cuts[r]]
    pmy = -RUN.Branches[r]['H.kin.secondary.Prec_y'][all_cuts[r]]
    pmz = RUN.Branches[r]['H.kin.secondary.Prec_z'][all_cuts[r]]
    
    pm = RUN.Branches[r]['H.kin.secondary.pmiss'][all_cuts[r]]
    
    PmPer = np.sqrt(pm**2 - pmz**2 - pmy**2)
    
    pmx_h = NORM[r]*B.histo(pmx, range=(-0.1,0.1), bins=50)
    pmy_h = NORM[r]*B.histo(pmy, range=(-0.1,0.1), bins=50)
    pmz_h = NORM[r]*B.histo(pmz, range=(-0.1,0.1), bins=50)
    
    Pmx[r] = pmx
    Pmy[r] = pmy
    Pmz[r] = pmz
    Pmx_wcuts_histo[r] = pmx_h
    Pmy_wcuts_histo[r] = pmy_h
    Pmz_wcuts_histo[r] = pmz_h
    
Pmx_sim = {}
Pmy_sim = {}
Pmz_sim = {}
Pmx_wcuts_histo_sim = {}
Pmy_wcuts_histo_sim = {}
Pmz_wcuts_histo_sim = {}

for s in SIMC.many:
    pmx = SIMC.Branches[s]['Pmx'][all_cuts_sim[s]]
    pmy = SIMC.Branches[s]['Pmy'][all_cuts_sim[s]]
    pmz = SIMC.Branches[s]['Pmz'][all_cuts_sim[s]]


    pmx_h = B.histo(pmx, range=(-0.1,0.1), bins=50,
                                  weights=WEIGHTS[s], calc_w2=True)
    pmy_h = B.histo(pmy, range=(-0.1,0.1), bins=50, 
                                  weights=WEIGHTS[s], calc_w2=True)
    pmz_h = B.histo(pmz, range=(-0.1,0.1), bins=50, 
                                  weights=WEIGHTS[s], calc_w2=True)
    Pmx_sim[s] = pmx
    Pmy_sim[s] = pmy
    Pmz_sim[s] = pmz
    Pmx_wcuts_histo_sim[s] = pmx_h
    Pmy_wcuts_histo_sim[s] = pmy_h
    Pmz_wcuts_histo_sim[s] = pmz_h
    
histos_to_plot = {'Pmx':Pmx_wcuts_histo, 'Pmy':Pmy_wcuts_histo, 
                  'Pmz':Pmz_wcuts_histo}

simhistos_to_plot = {'Pmx':Pmx_wcuts_histo_sim, 'Pmy':Pmy_wcuts_histo_sim, 
                     'Pmz':Pmz_wcuts_histo_sim}

#plotting and fitting the histos

diffs_Pmx = {}
diffs_Pmy = {}
diffs_Pmz = {}
fit_lims = {'Pmx':(-0.025,0.025,-0.025,0.025),'Pmy':(-0.05,0.05,-0.05,0.05),
            'Pmz':(-0.025,0.025,-0.025,0.025)}
for plot in histos_to_plot:
    fig = B.pl.figure(layout='constrained',figsize=(20,10))
    fig.suptitle(f'{plot}', fontsize=16)
    rows = 2
    cols = 4
    i = 1
    for r in RUN.many:      
        s = RUN.setting[r]
        h = histos_to_plot[plot][r]
        hsim = simhistos_to_plot[plot][s]
        
        ax = plt.subplot(rows,cols,i)
        
        hsim.plot_exp(c='#ed6792',marker='+',markersize=8,
                             capsize=0,mew=1,elinewidth=1,label='SIMC')
        hsim.fit(fit_lims[plot][0],fit_lims[plot][1],plot_fit=False)
        hsim.fit(hsim.mean.value-hsim.sigma.value*0.8,hsim.mean.value+hsim.sigma.value*0.8)
        hsim.plot_fit(color='#de425b')
        
        mu_s = hsim.mean.value
        
        h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                             capsize=0,mew=1,elinewidth=1,label='Data') 
        h.fit(fit_lims[plot][2],fit_lims[plot][3],plot_fit=False)
        h.fit(h.mean.value-h.sigma.value*0.8,h.mean.value+h.sigma.value*0.8)
        h.plot_fit(color='#1f77b4')
        
        mu = h.mean.value
    
        # if mu_s > mu:  
        #     diffs = mu_s - mu
        # else:
        #     diffs = mu - mu_s
        
        #DATA-SIMC
        diffs = mu - mu_s
        
        han, lab = fit_legend([h.fit_par,hsim.fit_par])
        
        if plot == 'Pmx':
            diffs_Pmx[r] = diffs
        elif plot == 'Pmy':
            diffs_Pmy[r] = diffs
        elif plot == 'Pmz':
            diffs_Pmz[r] = diffs   
        
        # #print(han,lab)
        
        diff_lab = f'DATA - SIMC = {diffs:.3e}'
        diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
        
        myhan = [han[1],han[2],han[0],han[3],diff_han]
        mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
        B.pl.legend(handles=myhan,labels=mylab,loc='best',fontsize='xx-small')
        
        bot,top = plt.ylim()
    
        # # plot formatting
        # B.pl.title(f'$\Delta th_p$ for {s}')
        # B.pl.xlabel('')
        # B.pl.ylabel('')
        # B.pl.ylim(0,top)
        # #B.pl.legend(handles=han,labels=lab,loc='best')
        # B.pl.xlabel('[rad]')
        B.pl.vlines([hsim.mean.value,h.mean.value], 0,[hsim.A.value,h.A.value],
                    colors=['#de425b','#1f77b4'],linestyles='--')

        ax.set_autoscaley_on(True)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'Run {r} {s}')
        i+=1

#%% pm diff vs delta
delta = np.array([-8,-4,0,4,8,12])
# diffs_sorted = np.array([diffs_Pmx[20841],diffs_Pmx[20846],diffs_Pmx[20851],
#                 diffs_Pmx[20858],diffs_Pmx[20861],diffs_Pmx[20869]])
diffs_sorted = np.array([diffs_Pmy[20841],diffs_Pmy[20846],diffs_Pmy[20851],
                diffs_Pmy[20858],diffs_Pmy[20861],diffs_Pmy[20869]])
# diffs_sorted = np.array([diffs_Pmz[20841],diffs_Pmz[20846],diffs_Pmz[20851],
#                 diffs_Pmz[20858],diffs_Pmz[20861],diffs_Pmz[20869]])

B.pl.figure()
B.pl.scatter(delta, diffs_sorted)
line = B.linefit(delta, diffs_sorted)

plt.title('DATA-SIMC dPf for each $\delta$ setting')
han,lab = fit_legend(line)
B.pl.legend(handles=han,labels=lab)
# B.pl.savefig('STEP1_DATA-SIMC_vs_delta.png')

#%% relative angle offsets

NORM = {}
for r in RUN.many:
    NORM[r] = get_norm(r)

WEIGHTS = {}    
for s in SIMC.many:    
    W = calc_weights(SIMC.Branches[s])
    WEIGHTS[s] = W

pxptar = {}
pyptar = {}
hxptar = {}
hyptar = {}
for r in RUN.many:
    pxptar[r] = RUN.Branches[r]['P.gtr.th']
    pyptar[r] = RUN.Branches[r]['P.gtr.ph']
    hxptar[r] = RUN.Branches[r]['H.gtr.th']
    hyptar[r] = RUN.Branches[r]['H.gtr.ph']

spxptar = {}
spyptar = {}
shxptar = {}
shyptar = {}
for s in SIMC.many:
    spxptar[s] = SIMC.Branches[s]['e_xptar']
    spyptar[s] = SIMC.Branches[s]['e_yptar']
    shxptar[s] = SIMC.Branches[s]['h_xptar']
    shyptar[s] = SIMC.Branches[s]['h_yptar']
#%% cuts for xptar/yptar
W_cut = C.WCUT(0.91,1.0,name='W_cut')
C.HCANA_names.update({'W_cut':'P.kin.primary.W'})
C.SIMC_names.update({'W_cut':'W'})

cuts_list = [W_cut,C.hms_delta]

cuts_to_apply = {}
for r in RUN.many:
    c_list = []
    print(f'Cuts for Run {r}')
    for cut in cuts_list:
        cut.init()
        br = RUN.Branches[r][C.HCANA_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
    
    # add collimator cut
    hxc = RUN.Branches[r]['H.extcor.xsieve']
    hyc = RUN.Branches[r]['H.extcor.ysieve']
    
    hcoll_cut = C.coll_cut(hxc, hyc, spec='HMS')
    c_list.append(hcoll_cut) 
    
    cuts_to_apply[r] = c_list
    
all_cuts = {}
for r in RUN.many:
    all_cuts_arr = cuts_to_apply[r][0]
    for arr in cuts_to_apply[r]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts[r] = all_cuts_arr
    
cuts_to_apply_sim = {}
for s in SIMC.many:
    c_list = []
    print(f'Cuts for SIMC setting {s}')
    for cut in cuts_list:
        cut.init()
        br = SIMC.Branches[s][C.SIMC_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
        
    cuts_to_apply_sim[s] = c_list

all_cuts_sim = {}
for s in SIMC.many:
    all_cuts_arr = cuts_to_apply_sim[s][0]
    for arr in cuts_to_apply_sim[s]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts_sim[s] = all_cuts_arr    

#%% make 1d histograms
pxptar_h = {}
pyptar_h = {}
hxptar_h = {}
hyptar_h = {}  
for r in RUN.many:
    pxptar_h[r] = NORM[r]*B.histo(pxptar[r][all_cuts[r]], 
                                  range=(-0.05,0.05), bins=100)
    pyptar_h[r] = NORM[r]*B.histo(pyptar[r][all_cuts[r]], 
                                  range=(-0.03,0.03), bins=100)
    hxptar_h[r] = NORM[r]*B.histo(hxptar[r][all_cuts[r]], 
                                  range=(-0.1,0.1), bins=100)
    hyptar_h[r] = NORM[r]*B.histo(hyptar[r][all_cuts[r]], 
                                  range=(-0.04,0.04), bins=100)

spxptar_h = {}
spyptar_h = {}
shxptar_h = {}
shyptar_h = {}  
for s in SIMC.many:
    spxptar_h[s] = B.histo(spxptar[s][all_cuts_sim[s]],range=(-0.05,0.05), 
                           bins=100,
                           weights=WEIGHTS[s][all_cuts_sim[s]],calc_w2=True)
    spyptar_h[s] = B.histo(spyptar[s][all_cuts_sim[s]], range=(-0.03,0.03), 
                           bins=100,
                           weights=WEIGHTS[s][all_cuts_sim[s]],calc_w2=True)
    shxptar_h[s] = B.histo(shxptar[s][all_cuts_sim[s]], range=(-0.1,0.1), 
                           bins=100,
                           weights=WEIGHTS[s][all_cuts_sim[s]],calc_w2=True)
    shyptar_h[s] = B.histo(shyptar[s][all_cuts_sim[s]], range=(-0.04,0.04),
                           bins=100,
                           weights=WEIGHTS[s][all_cuts_sim[s]],calc_w2=True)
#%% calculate DATA/SIMC ratios
pxptar_ratio = {}
pyptar_ratio = {}
hxptar_ratio = {}
hyptar_ratio = {}
for r in RUN.many:
    s = RUN.setting[r]
    #set 0 to 1e4
    set_zero(spxptar_h[s].bin_content,1e5)
    set_zero(spyptar_h[s].bin_content,1e5)
    set_zero(shxptar_h[s].bin_content,1e5)
    set_zero(shyptar_h[s].bin_content,1e5)
    
    pxptar_ratio[r] = pxptar_h[r]/spxptar_h[s]
    pyptar_ratio[r] = pyptar_h[r]/spyptar_h[s]
    hxptar_ratio[r] = hxptar_h[r]/shxptar_h[s]
    hyptar_ratio[r] = hyptar_h[r]/shyptar_h[s]
    
    # #set nans to 0
    # set_nan(pxptar_ratio[r].bin_content)
    # set_nan(pyptar_ratio[r].bin_content)
    # set_nan(hxptar_ratio[r].bin_content)
    # set_nan(hyptar_ratio[r].bin_content)
    
#%% plot 1d histos
hist_to_plot = {'SHMS xptar':pxptar_h,'SHMS yptar':pyptar_h,
                'HMS xptar':hxptar_h,'HMS yptar':hyptar_h}
sim_hist_to_plot = {'SHMS xptar':spxptar_h,'SHMS yptar':spyptar_h,
                'HMS xptar':shxptar_h,'HMS yptar':shyptar_h}
hyptar_diffs = []
for plot in hist_to_plot:
    fig = B.pl.figure(layout='constrained',figsize=(20,10))
    fig.suptitle(f'{plot}', fontsize=16)
    rows = 2
    cols = 4
    i = 1
    for r in RUN.many:       
        s = RUN.setting[r]
        ax = plt.subplot(rows,cols,i)
        h = hist_to_plot[plot][r]
        hsim = sim_hist_to_plot[plot][s]
        
        h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                             capsize=0,mew=1,elinewidth=1,label='Data')         
        hsim.plot_exp(c='#ed6792',marker='+',markersize=8,
                             capsize=0,mew=1,elinewidth=1,label='SIMC')

        if plot == 'HMS yptar':
            hsim.fit(plot_fit=False)
            hsim.fit(hsim.mean.value-1.2*hsim.sigma.value,hsim.mean.value+1.2*hsim.sigma.value)
            hsim.plot_fit(color='#de425b')
        
            mu_s = hsim.mean.value
        
            h.fit(plot_fit=False)
            h.fit(h.mean.value-1.2*h.sigma.value,h.mean.value+1.2*h.sigma.value)
            h.plot_fit(color='#1f77b4')
            
            mu = h.mean.value
        
            diffs = mu_s - mu

            hyptar_diffs.append(diffs)
        
            han, lab = fit_legend([h.fit_par,hsim.fit_par])
        
            #print(han,lab)
        
            diff_lab = f'SIMC - DATA = {diffs:.3e}'
            diff_han = mlines.Line2D([], [], color='None', label=diff_lab)
        
            myhan = [han[1],han[2],han[0],han[3],diff_han]
            mylab = [lab[1],lab[2],lab[0],lab[3],diff_lab]
            # B.pl.legend(handles=myhan,labels=mylab,loc='best')
        
        
            B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                    colors=['#de425b','#1f77b4'],linestyles='--')
        
        i+=1
        ax.set_autoscaley_on(True)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'Run {r} {s}')  

#%% plot ratios
hist_to_plot = {'SHMS xptar DATA/SIMC ratio':[pxptar_ratio,(-0.05,0.05)],
                'SHMS yptar DATA/SIMC ratio':[pyptar_ratio,(-0.03,0.03)],
                'HMS xptar DATA/SIMC ratio':[hxptar_ratio,(-0.1,0.1)],
                'HMS yptar DATA/SIMC ratio':[hyptar_ratio,(-0.04,0.04)]}

for plot in hist_to_plot:
    fig = B.pl.figure(layout='constrained',figsize=(20,10))
    fig.suptitle(f'{plot}', fontsize=16)
    rows = 2
    cols = 4
    i = 1
    for r in RUN.many:       
        s = RUN.setting[r]
        ax = plt.subplot(rows,cols,i)
        h = hist_to_plot[plot][0][r]
        
        h.plot_exp(c='#7b8fd4',marker='+',markersize=8,
                             capsize=0,mew=1,elinewidth=1)                 
        B.pl.hlines([1,0.9,0.8],hist_to_plot[plot][1][0],
                    hist_to_plot[plot][1][1],
                    linestyles='--',colors=['g','orange','r'])
        B.pl.text(hist_to_plot[plot][1][0]+0.005,0.91,'90%',color='orange')
        B.pl.text(hist_to_plot[plot][1][0]+0.005,0.81,'80%',color='r')
     
        i+=1
        ax.set_autoscaley_on(True)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'Run {r} {s}')  
