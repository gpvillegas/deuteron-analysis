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
from matplotlib import colormaps as cmp
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

def combine_histos(histos,many=[]):
    loop1 = True 
    for m in many:
        if loop1:
            hsum = histos[m]
            loop1 = False
        else:
            hsum += histos[m]
    
    return hsum    
        

        
    
#%% Declare constants
Eb = 10.542         # GeV
Mp = 0.938272       # GeV

dtr = np.pi/180

#%% Load hcana variables

# kin_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'kin\''))
# gtr_branches = get_list(db.retrieve('deuteron_db.db', 'Branches', 'TTree_Branches',
#                            where = 'Type = \'gtr\''))
# branches_sel = kin_branches + gtr_branches

kin_var = ['P.kin.primary.scat_ang_rad','H.kin.secondary.xangle','P.gtr.p',
           'H.gtr.p','P.kin.primary.x_bj']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

sieve_branches = ['H.extcor.xsieve', 'H.extcor.ysieve',
                  'P.extcor.xsieve', 'P.extcor.ysieve']

br_sel = kin_var + acc_var + calPID_var + sieve_branches

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','theta_e','theta_p',
               'e_pf','h_pf','e_xptar','e_yptar','h_xptar','h_yptar']
#%% load data root files
# make sure files are where they should be!
RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                  select_branches={'T':br_sel})
#%%
# load simc root files
# make sure files are where they should be!
SIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                   select_branches={'SNT':br_sel_SIMC},simc_type='-')     
#%% Assign measured quatities

th_e = {}
xangle = {}
th_p_meas = {}
kf_meas = {}
Pf_meas = {}
for r in RUN.many:
    # electron scattering angle
    th_e[r] = RUN.Branches[r]['P.kin.primary.scat_ang_rad']

    # angle of detected particle (proton) with scattered electron
    xangle[r] = RUN.Branches[r]['H.kin.secondary.xangle']

    # detected proton scattering angle
    th_p_meas[r] = xangle[r] - th_e[r]

    # measured electron final momentum
    kf_meas[r] = RUN.Branches[r]['P.gtr.p']

    # measured proton final momentum
    Pf_meas[r] = RUN.Branches[r]['H.gtr.p']

th_e_sim = {}
th_p_meas_sim = {}
kf_meas_sim = {}
Pf_meas_sim = {}
for s in SIMC.many:
    th_e_sim[s] = SIMC.Branches[s]['theta_e']

    th_p_meas_sim[s] = SIMC.Branches[s]['theta_p']

    kf_meas_sim[s] = SIMC.Branches[s]['e_pf']*(1/1e3)

    Pf_meas_sim[s] = SIMC.Branches[s]['h_pf']*(1/1e3)

# Declare calculated quatities

kf_calc = {}
Pf_x = {}
Pf_y = {}
Pf_z = {}
Pf_calc = {}
th_p_calc = {}

for r in RUN.many:
    # electron final momentum in GeV
    kf_calc[r] = (Mp*Eb)/(Mp + 2.*Eb*np.sin(0.5*th_e[r])*np.sin(0.5*th_e[r]))

    # proton final momentum components in GeV
    Pf_x[r] = kf_calc[r]*np.sin(th_e[r])
    Pf_y[r] = 0.
    Pf_z[r] = Eb - kf_calc[r]*np.cos(th_e[r])

    # proton final momentum in GeV
    Pf_calc[r] = np.sqrt(Pf_x[r]*Pf_x[r] + Pf_y[r]*Pf_y[r] + Pf_z[r]*Pf_z[r])

    # proton recoil angle
    np.seterr(invalid='ignore') # ignore RunTimeWarning due to 0/0 division
            
    th_p_calc[r] = np.arctan(Pf_x[r]/Pf_z[r])

kf_calc_sim = {}
Pf_x_sim = {}
Pf_y_sim = {}
Pf_z_sim = {}
Pf_calc_sim = {}
th_p_calc_sim = {}

for s in SIMC.many:
    kf_calc_sim[s] =\
        (Mp*Eb)/(Mp + 2.*Eb*np.sin(0.5*th_e_sim[s])*np.sin(0.5*th_e_sim[s]))

    Pf_x_sim[s] = kf_calc_sim[s]*np.sin(th_e_sim[s])
    Pf_y_sim[s] = 0.
    Pf_z_sim[s] = Eb - kf_calc_sim[s]*np.cos(th_e_sim[s])

    Pf_calc_sim[s] =\
        np.sqrt(Pf_x_sim[s]*Pf_x_sim[s] + Pf_y_sim[s]*Pf_y_sim[s] +\
                Pf_z_sim[s]*Pf_z_sim[s])

    np.seterr(invalid='ignore') # ignore RunTimeWarning due to 0/0 division
    th_p_calc_sim[s] = np.arctan(Pf_x_sim[s]/Pf_z_sim[s])

# Calculated - Measured quatities

dkf = {}
dPf = {}
dth_p = {}
for r in RUN.many:
    dkf[r] = kf_calc[r] - kf_meas[r]
    dPf[r] = Pf_calc[r] - Pf_meas[r]
    dth_p[r] = th_p_calc[r] - th_p_meas[r]

dkf_sim = {}
dPf_sim = {}
dth_p_sim = {}
for s in SIMC.many:
    dkf_sim[s] = kf_calc_sim[s] - kf_meas_sim[s]
    dPf_sim[s] = Pf_calc_sim[s] - Pf_meas_sim[s]
    dth_p_sim[s] = th_p_calc_sim[s] - th_p_meas_sim[s]

#%%
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
  
#
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

#%% Plot differences with cuts; apply cuts
for r in RUN.many:
    dkf[r] = dkf[r][all_cuts[r]]
    dPf[r] = dPf[r][all_cuts[r]]
    dth_p[r] = dth_p[r][all_cuts[r]]
    
for s in SIMC.many:
    dkf_sim[s] = dkf_sim[s][all_cuts_sim[s]]
    dPf_sim[s] = dPf_sim[s][all_cuts_sim[s]]
    dth_p_sim[s] = dth_p_sim[s][all_cuts_sim[s]]
    
#%% calculate normalization factor for data and locate weigths for SIMC
NORM = {}
for r in RUN.many:
    NORM[r] = get_norm(r)

W = {}    
for s in SIMC.many:    
    WEIGHTS = calc_weights(SIMC.Branches[s])
    W[s] = WEIGHTS[all_cuts_sim[s]]
    
#%% delta theta_p; making the histos
dthp_wcuts_histos = {}
dkf_wcuts_histos = {}
dPf_wcuts_histos = {}

for r in RUN.many:
    dthp_wcuts_histos[r] = NORM[r]*B.histo(dth_p[r],(-0.015,0.015),100)
    dkf_wcuts_histos[r] = NORM[r]*B.histo(dkf[r],(-0.1,0.1),100)
    dPf_wcuts_histos[r] = NORM[r]*B.histo(dPf[r],(-0.1,0.1),100)
    
dthp_wcuts_simhistos = {}
dkf_wcuts_simhistos = {}
dPf_wcuts_simhistos = {}

for s in SIMC.many:
    dthp_wcuts_simhistos[s] = B.histo(dth_p_sim[s],(-0.015,0.015),100,
                                      weights=W[s],calc_w2=True)
    dkf_wcuts_simhistos[s] = B.histo(dkf_sim[s],(-0.1,0.1),100,
                                     weights=W[s],calc_w2=True)
    dPf_wcuts_simhistos[s] = B.histo(dPf_sim[s],(-0.1,0.1),100,
                                     weights=W[s],calc_w2=True)
    
#%% plotting and fitting the histos
hist_to_plot = {'$\Delta \\theta_p$':dthp_wcuts_histos,
                '$\Delta k_f$':dkf_wcuts_histos,
                '$\Delta P_f$':dPf_wcuts_histos}
sim_hist_to_plot = {'$\Delta \\theta_p$':dthp_wcuts_simhistos,
                '$\Delta k_f$':dkf_wcuts_simhistos,
                '$\Delta P_f$':dPf_wcuts_simhistos}
diffs_thp = {}
diffs_kf = {}
diffs_Pf = {}
for plot in hist_to_plot:
    fig = B.pl.figure(layout='constrained',figsize=(20,10))
    fig.suptitle(f'{plot}', fontsize=16)
    rows = 2
    cols = 4
    i = 1
    for r in RUN.many:      
        s = RUN.setting[r]
        h = hist_to_plot[plot][r]
        hsim = sim_hist_to_plot[plot][s]
        
        ax = plt.subplot(rows,cols,i)
        
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
    
        # if mu_s > mu:  
        #     diffs = mu_s - mu
        # else:
        #     diffs = mu - mu_s
        
        #DATA-SIMC
        diffs = mu - mu_s
        
        han, lab = fit_legend([h.fit_par,hsim.fit_par])
        
        if plot == '$\Delta \\theta_p$':
            diffs_thp[r] = diffs
        elif plot == '$\Delta k_f$':
            diffs_kf[r] = diffs
        elif plot == '$\Delta P_f$':
            diffs_Pf[r] = diffs   
        
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
        B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                    colors=['#de425b','#1f77b4'],linestyles='--')

        ax.set_autoscaley_on(True)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'Run {r} {s}')
        i+=1

#%%
        save_where = './step2/'
        fname = save_where + f'{plot}_{r}_{s}_2.png'
        extent = full_extent(ax).transformed(fig.dpi_scale_trans.inverted())
        B.pl.savefig(fname,bbox_inches=extent)
    
#%% plot kf_calc vs delta

kf_calc_cut = {}
kf_meas_cut = {}
delta = {}
xbj = {}
for m in RUN.many:
    
    kf_calc_cut[m] = kf_calc[m][all_cuts[m]]
    kf_meas_cut[m] = kf_meas[m][all_cuts[m]]
    
    d = RUN.Branches[m]['P.gtr.dp']
    delta[m] = d[all_cuts[m]]
    
    xbj[m] = RUN.Branches[m]['P.kin.primary.x_bj'][all_cuts[m]]
                 
# make and store 2d histograms
kf_calc_vs_delta_h = {}
xbj_vs_delta_h = {}
for m in RUN.many:
    kf = kf_calc_cut[m]
    d = delta[m]
    x = xbj[m]
    
    h = B.histo2d(d,100*((kf/8.55)-1)-d,range=[(-10.5,15),(-5,5)],bins=100)
    h2 = B.histo2d(d,x,range=[(-10.5,15),(0.75,1.25)],bins=100)
    
    kf_calc_vs_delta_h[m] = h
    xbj_vs_delta_h[m] = h2
    
    if m == 20840:
        y_fit = 100*((kf/8.55)-1)-d
        x_fit = d
    else:
        y_fit = np.append(y_fit,100*((kf/8.55)-1)-d)
        x_fit = np.append(x_fit,d)
    
# combine all histos

kf_calc_vs_delta_h_all = combine_histos(kf_calc_vs_delta_h,RUN.many) 
xbj_vs_delta_h_all = combine_histos(xbj_vs_delta_h,RUN.many)

range_y = (y_fit >= -0.8) & (y_fit <= 0.1)
range_x = (x_fit <= 13)
range_fit = range_x & range_y

fit = B.polyfit(x_fit[range_fit], y_fit[range_fit],order=1)

B.pl.figure()
kf_calc_vs_delta_h_all.plot(colormap=cmp['viridis'],logz=True)
fit.plot(color='r')
B.pl.title('kf_calc vs $\delta$')
B.pl.xlabel('[%]')
B.pl.ylabel('[%]')

B.pl.figure()
xbj_vs_delta_h_all.plot(colormap=cmp['viridis'],logz=True)
B.pl.title('x_bj vs $\delta$')
B.pl.xlabel('[%]')
B.pl.ylabel('[]')

#%% DATA - SIMC differences

delta = np.array([-8,-4,0,4,8,12])
# diffs_sorted = np.array([diffs_kf[20841],diffs_kf[20846],diffs_kf[20851],
#                 diffs_kf[20858],diffs_kf[20861],diffs_kf[20869]])
# diffs_sorted = np.array([diffs_thp[20841],diffs_thp[20846],diffs_thp[20851],
#                 diffs_thp[20858],diffs_thp[20861],diffs_thp[20869]])
diffs_sorted = np.array([diffs_Pf[20841],diffs_Pf[20846],diffs_Pf[20851],
                diffs_Pf[20858],diffs_Pf[20861],diffs_Pf[20869]])

B.pl.figure()
plt.scatter(delta, diffs_sorted)
plt.title('dpf DATA-SIMC vs $\delta$')

fit1 = B.linefit(delta, diffs_sorted)

B.pl.text(2.5,-0.0065,f'slope = {fit1.slope:.4e}$\pm${fit1.sigma_s:.4e}')
B.pl.text(2.5,-0.0070,f'offset = {fit1.offset:.4e}$\pm${fit1.sigma_o:.4e}')

# B.pl.savefig('STEP1_DATA-SIMC_vs_delta.png') 
#%%
kf_nom = 8.55 #[GeV]
pf_nom = np.array([3.499,3.145,2.783,2.417])
the_nom = np.array([14.145,12.940,11.700,10.430,9.120,7.700])
the_nom_rad = the_nom*np.pi/180

diff_thp = np.array([0.00012,0.00019,0.00026,0.00062,0.00063,0.00553])

dthp_dthe = (kf_nom*kf_nom - Eb*kf_nom*np.cos(the_nom_rad))/\
                Eb*Eb - 2.*Eb*kf_nom*np.cos(the_nom_rad) + kf_nom*kf_nom
 
corr_offset = diff_thp/dthp_dthe

#%%
delta = np.array([-8,-4,0,4,8,12])
kf_corr = 8.5201
the_nom = np.array([14.145,12.940,11.700,10.430,9.120,7.700])
Delta_pf = np.array([diffs_Pf[20841],diffs_Pf[20846],diffs_Pf[20851],
                diffs_Pf[20858],diffs_Pf[20861],diffs_Pf[20840]])
pfx = kf_corr*np.sin(the_nom*dtr)
pfz = Eb - kf_corr*np.cos(the_nom*dtr)
pfcalc = np.sqrt(pfx*pfx + pfz*pfz)

dkf_dthe = -(Mp*Eb*Eb*np.sin(0.5*the_nom*dtr)*np.cos(0.5*the_nom*dtr))/\
    (Mp+2*Eb*np.sin(0.5*the_nom*dtr)*np.sin(0.5*the_nom*dtr))**2

dPfx_dthe = dkf_dthe*np.sin(the_nom*dtr) + kf_corr*np.cos(the_nom*dtr)
dPfz_dthe = -dkf_dthe*np.cos(the_nom*dtr) + kf_corr*np.sin(the_nom*dtr)

dPf_dthe = (pfx*dPfx_dthe + pfz*dPfz_dthe)/pfcalc

the_corr = Delta_pf/dPf_dthe

print(the_corr/dtr)
B.pl.figure()
B.pl.plot(delta,the_corr/dtr,'.')
line = B.linefit(delta,the_corr/dtr,plot_fit=False)
line.plot(color='r')
B.pl.text(1.5,-0.03,f'slope = {line.slope:.4e}$\pm${line.sigma_s:.4e}')
B.pl.text(1.5,-0.033,f'offset = {line.offset:.4e}$\pm${line.sigma_o:.4e}')
#%% momentum component alignment
Pmx = {}
Pmy = {}
Pmz = {}
Pmx_wcuts_histo = {}
Pmy_wcuts_histo = {}
Pmz_wcuts_histo = {}

for r in RUN.many:
    pmx = RUN.Branches[r]['H.kin.secondary.pmiss_x'][all_cuts[r]]
    pmy = RUN.Branches[r]['H.kin.secondary.pmiss_y'][all_cuts[r]]
    pmz = RUN.Branches[r]['H.kin.secondary.pmiss_z'][all_cuts[r]]
    
    pmx_h = NORM[r]*B.histo(pmx, range=(-0.025,0.050), bins=100)
    pmy_h = NORM[r]*B.histo(pmy, range=(-0.025,0.050), bins=100)
    pmz_h = NORM[r]*B.histo(pmz, range=(-0.025,0.050), bins=100)
    
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


    pmx_h = B.histo(pmx, range=(-0.025,0.050), bins=100,
                                  weights=W[s], calc_w2=True)
    pmy_h = B.histo(pmy, range=(-0.025,0.050), bins=100, 
                                  weights=W[s], calc_w2=True)
    pmz_h = B.histo(pmz, range=(-0.025,0.050), bins=100, 
                                  weights=W[s], calc_w2=True)
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

#%% plotting and fitting the histos

diffs_Pmx = {}
diffs_Pmy = {}
diffs_Pmz = {}
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
        B.pl.vlines([hsim.mean.value,h.mean.value], 1,[hsim.A.value,h.A.value],
                    colors=['#de425b','#1f77b4'],linestyles='--')

        ax.set_autoscaley_on(True)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f'Run {r} {s}')
        i+=1