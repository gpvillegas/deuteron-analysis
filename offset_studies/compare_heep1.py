#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:39:09 2024

@author: gvill
"""
import numpy as np
import LT.box as B
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


# import deut_analysis_tools as dat
import database_operations as db
import data_init as D
import cut_handler as C
from matplotlib import colormaps as cmp

'''
Comparison plots btw SIMC and heep data
'''
#%% functions live here
def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(8,7)
    
    
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

###
# function to normalize histograms by efficiencies and getting the total
# charge of the runs given, should I keep the charge per run as well?
###

def normalize_histos(histos,many=[]):
    if many:
        charge = []
        h_enorm = {}
        for m in many:
            h = histos[m]
            enorm = D.get_eff_norm(m)
            charge.append(D.get_charge_norm(m))
            
            h_enorm[m] = enorm*h
            
        charge_tot = np.sum(np.array(charge))
        return charge_tot,h_enorm
    else:
        print('Need to know the run # to get the normalization factors')

###
# function to combine histograms from different runs, will be normalized by
# efficiencies per run, then divided by total charge of combined runs.
###        
        
def combine_histos(histos,many=[]):
    tot_q, enorm_h = normalize_histos(histos,many)
    
    loop1 = True 
    for m in many:
        if loop1:
            hsum = enorm_h[m]
            loop1 = False
        else:
            hsum += enorm_h[m]
     
    qnorm_hsum = hsum*(1/tot_q)
    
    return qnorm_hsum
   
    
#%% load selection of branches for this specific analysis
deut_db = 'deuteron_db.db' 

kin_var = ['P.kin.primary.W','H.kin.secondary.emiss','P.kin.primary.Q2']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

coll_var = ['H.extcor.xsieve','H.extcor.ysieve',
            'P.extcor.xsieve','P.extcor.ysieve']

br_sel = kin_var + acc_var + calPID_var + ['CTime.epCoinTime_ROC2'] + coll_var

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','W','Em','e_xptar',
               'e_yptar','h_xptar','h_yptar']

#%% load root files
root_DIR = "/media/gvill/Gema's T7/ROOTfiles/heepcoin_hyptaroff_hthetacentraloff/"
# root_DIR = "/media/gvill/Gema's T7/ROOTfiles/deep_testing0/"
RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                    select_branches={'T':br_sel})
                    #,ROOTfiles_path=root_DIR)

# root_DIR = "/media/gvill/Gema's T7/ROOTfiles/heepsingles_deltaoptim2/"
# R = [20843,20844,20847,20856,20859,20862,20867]
# RUN = D.DATA_INIT(data_type='deut23_data', run=R, 
#                     select_branches={'T':br_sel},ROOTfiles_path=root_DIR)

# R = [20871,20873,20888,20958]
# R = [20871,20872] # pm_120
# R = [20873,20874,20875,20876,20877,20878,20880,20881,20882,20883,21076,
#      21078,21079,21080,21081,21082,21083,21084,21085,21086,21087,21088,
#      21089,21090,21091,21092,21093,21094,21095,21096,21097,21098,21099,
#      21100,21101,21102] # pm_580
# R = [20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,20896,
#      20897,20898,20899,20900,20901] # pm_800
# R = [20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,20896,
#      20897,20898,20899,20900,20901,20902,20903,20904,20905,20907,20908,
#      20909,20910,20911,20912,20913,20914,20915,20916,20921,20922,20923,
#      20924,20925,20926,20927,20928,20929,20930,20931,20932,20933,20934,
#      20935,20936,20937,20938,20939,20940,20941,20942,20943,20944,20945,
#      20949,20950,20951,20953,20954,20955,20956] # pm_800

# RUN = D.DATA_INIT(data_type='deut23_data', run=R, 
#                     select_branches={'T':br_sel},ROOTfiles_path=root_DIR)

# RUN = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
#                     select_branches={'T':br_sel})

# run = D.DATA_INIT(data_type='deut23_data', run=20851, 
#                     select_branches=br_sel)

#%%
# SIMC files
worksim_DIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_hyptaroff_hthetacentraloff/"
# worksim_DIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/deep_testing0/"

SIMC = D.DATA_INIT(data_type='SIMC', kin_study='heep_coin',
                    select_branches={'SNT':br_sel_SIMC},simc_type='-')
                    #,SIMC_ROOTfiles_path=worksim_DIR)

# Sheep = D.DATA_INIT(data_type='SIMC', kin_study='deep',
#                    simc_type = 'jmlfsi_rad',select_branches={'SNT':br_sel_SIMC})
# SIMC = D.DATA_INIT(data_type='SIMC', kin_study='deep', setting='pm_120',
#                    simc_type = 'jmlfsi_rad',select_branches={'SNT':br_sel_SIMC},
#                    SIMC_ROOTfiles_path=worksim_DIR)

# SIMC = D.DATA_INIT(data_type='SIMC', kin_study='deep', setting='pm_580',
#                    simc_type = 'jmlfsi_rad',select_branches={'SNT':br_sel_SIMC})
# SIMC = D.DATA_INIT(data_type='SIMC', kin_study='deep', setting='pm_800',
#                    simc_type = 'jmlfsi_rad',select_branches={'SNT':br_sel_SIMC})

# SIMC = D.DATA_INIT(data_type='SIMC', kin_study='heep_coin',
#                    simc_type = '-',select_branches={'SNT':br_sel_SIMC})

# worksim_DIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/heepsingles_deltaoptim2/"
# SIMC = D.DATA_INIT(data_type='SIMC', kin_study='heep_singles',
#                    simc_type = '-',select_branches={'SNT':br_sel_SIMC},
#                    SIMC_ROOTfiles_path=worksim_DIR)


#%% declare wanted variables

W = {}
Em = {}
e_delta = {}
h_delta = {}
for r in RUN.many:
    W[r] = RUN.Branches[r]['P.kin.primary.W']
    Em[r] = RUN.Branches[r]['H.kin.secondary.emiss']
    e_delta[r] = RUN.Branches[r]['P.gtr.dp']
    h_delta[r] = RUN.Branches[r]['H.gtr.dp']

SW = {}
SEm = {}
Se_delta = {}
Sh_delta = {}
for s in SIMC.many:
    SW[s] = SIMC.Branches[s]['W']
    SEm[s] = SIMC.Branches[s]['Em']
    Se_delta[s] = SIMC.Branches[s]['e_delta']
    Sh_delta[s] = SIMC.Branches[s]['h_delta']
    
# SW = SIMC.Branches['W']
# SEm = SIMC.Branches['Em']        

#%% Define cuts for data
for m in RUN.many:
    coinTime = RUN.Branches[m]['CTime.epCoinTime_ROC2']
    coinTime_corr = D.get_coinTimeOffset(coinTime)
    
    RUN.Branches[m].update({'CTime.epCoinTime_ROC2_corr':coinTime_corr})

cuts_list = C.acceptance_cuts 
# cuts_list = C.acceptance_cuts  

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
#%%
# Define cuts for SIMC
cuts_list = C.acceptance_cuts

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

#%%
# singular run
cuts_list = C.acceptance_cuts
    
c_list = []
for cut in cuts_list:
    cut.init()
    br = SIMC.Branches[C.SIMC_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    c_list.append(cut_array)
    
all_cuts_sim = c_list[0]
for arr in c_list:    
    all_cuts_sim = all_cuts_sim & arr     

#%% apply cuts to variables
W_c = {}
Em_c = {}
edelta_c = {}
hdelta_c = {}
for r in RUN.many:
    W_c[r] = W[r][all_cuts[r]]
    Em_c[r] = Em[r][all_cuts[r]]
    edelta_c[r] = e_delta[r][all_cuts[r]]
    hdelta_c[r] = h_delta[r][all_cuts[r]]
    
SW_c = {}
SEm_c = {}
Sedelta_c = {}
Shdelta_c = {}
for s in SIMC.many:
    SW_c[s] = SW[s][all_cuts_sim[s]]
    SEm_c[s] = SEm[s][all_cuts_sim[s]]
    Sedelta_c[s] = Se_delta[s][all_cuts_sim[s]]
    Shdelta_c[s] = Sh_delta[s][all_cuts_sim[s]]
    
# SW_c = SW[all_cuts_sim]
# SEm_c = SEm[all_cuts_sim]   
    
#%% calculate normalization factors for data runs and get weights for SIMC
NORM = {}
for r in RUN.many:
    enorm = D.get_eff_norm(r)
    # enorm = D.get_eff_norm(r,run_type='singles')
    charge = D.get_charge_norm(r)
    
    NORM[r] = (1/charge)*enorm
    

WEIGHTS = {}    
for s in SIMC.many:    
    weights = D.calc_weights(SIMC.Branches[s])
    WEIGHTS[s] = weights[all_cuts_sim[s]]
#%% make histograms        
    
W_histos_wcuts = {}
Em_histos_wcuts = {}

ed_W_histos_wcuts = {}
hd_W_histos_wcuts = {}
ed_Em_histos_wcuts = {}
hd_Em_histos_wcuts = {}

edelta_ran = {20840:(9,14),20841:(-10,-4),20846:(-7,0),20851:(-3,4),
              20858:(1,8),20861:(5,12),20868:(9,14),20869:(9,14)}
for r in RUN.many:
    W_histos_wcuts[r] = NORM[r]*B.histo(W_c[r],range = (0.7,1.2),bins = 200)
    Em_histos_wcuts[r] = NORM[r]*B.histo(Em_c[r],range = (-0.1,0.1),bins = 200)
    
    ed_W_histos_wcuts[r] = B.histo2d(W_c[r],edelta_c[r],
                                             range = [(0.7,1.2),edelta_ran[r]],
                                                      bins = 200)
    hd_W_histos_wcuts[r] = B.histo2d(W_c[r],hdelta_c[r],
                                             range = [(0.7,1.2),(-10,10)],
                                                      bins = 200)
    ed_Em_histos_wcuts[r] = B.histo2d(Em_c[r],edelta_c[r],
                                             range = [(-0.1,0.1),edelta_ran[r]],
                                                      bins = 200)
    hd_Em_histos_wcuts[r] = B.histo2d(Em_c[r],hdelta_c[r],
                                             range = [(-0.1,0.1),(-10,10)],
                                                      bins = 200)

# make simc histos

SIMC_W_histos_wcuts = {}
SIMC_Em_histos_wcuts = {}

SIMC_ed_W_histos_wcuts = {}
SIMC_hd_W_histos_wcuts = {}
SIMC_ed_Em_histos_wcuts = {}
SIMC_hd_Em_histos_wcuts = {}
for s in SIMC.many:
    SIMC_W_histos_wcuts[s] = B.histo(SW_c[s],range = (0.7,1.2),bins = 200,
                                     weights=WEIGHTS[s], calc_w2=True)
    SIMC_Em_histos_wcuts[s] = B.histo(SEm_c[s],range = (-0.1,0.1),bins = 200,
                                      weights=WEIGHTS[s], calc_w2=True)
    
    SIMC_ed_W_histos_wcuts[s] = B.histo2d(SW_c[s],Sedelta_c[s],
                                          range=[(0.7,1.2),(-15,15)],
                                          bins=200)
    SIMC_hd_W_histos_wcuts[s] = B.histo2d(SW_c[s],Shdelta_c[s],
                                          range=[(0.7,1.2),(-15,15)],
                                          bins=200)
    SIMC_ed_Em_histos_wcuts[s] = B.histo2d(SEm_c[s],Sedelta_c[s],
                                          range=[(-0.1,0.1),(-15,15)],
                                          bins=200)
    SIMC_hd_Em_histos_wcuts[s] = B.histo2d(SEm_c[s],Shdelta_c[s],
                                          range=[(-0.1,0.1),(-15,15)],
                                          bins=200)

#%% make histograms: deep runs that need to be combined

W_histos_wcuts = {}
Em_histos_wcuts = {}
for r in RUN.many:
    W_histos_wcuts[r] = B.histo(W_c[r],range = (0,2),bins = 25)
    Em_histos_wcuts[r] = B.histo(Em_c[r],range = (-0.1,0.1),bins = 15)

    
weights = D.calc_weights(SIMC.Branches)
WEIGHTS = weights[all_cuts_sim]

SIMC_W_histos_wcuts = B.histo(SW_c,range = (0,2),bins = 25,
                                 weights=WEIGHTS, calc_w2=True)
SIMC_Em_histos_wcuts = B.histo(SEm_c,range = (-0.1,0.1),bins = 15,
                                  weights=WEIGHTS, calc_w2=True)    

# combine histos

W_all = combine_histos(W_histos_wcuts,RUN.many)
Em_all = combine_histos(Em_histos_wcuts,RUN.many)
    
#%% plot histos: W
W_means_data = {}
W_means_SIMC = {}

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('Invariant Mass (with cuts)', fontsize=16)
rows = 2
cols = 4
i = 1
for r in RUN.many:
# for s in SIMC.many:    
    s = RUN.setting[r]
    ax = plt.subplot(rows,cols,i)
    h = W_histos_wcuts[r]
    
    #trying correction factor
    # h.bin_center += f_corr[r]
    
    h.fit(plot_fit=False)
    hmean = h.mean.value
    hsigma = h.sigma.value
    h.fit(hmean-hsigma,hmean+hsigma,plot_fit=False)
    hmean = h.mean.value
    W_means_data[r] = hmean
    
    hS = SIMC_W_histos_wcuts[s]
    hS.fit(plot_fit=False)
    hSmean = hS.mean.value
    hSsigma = hS.sigma.value
    hS.fit(hSmean-hSsigma,hSmean+hSsigma,plot_fit=False)
    hSmean = hS.mean.value
    W_means_SIMC[s] = hSmean
    
    h.plot(hatch='......',facecolor='white',edgecolor='#7b8fd4')
    h.plot_fit(color='#144bcc')
    hS.plot(filled=False, color='#ed6792')
    hS.plot_fit(color='#ed1866')
    
    means_diff = hmean - hSmean
    handles, labels = ax.get_legend_handles_labels()

    diff_lab = f'DATA - SIMC = {means_diff:.3e}'
    diff_han = mlines.Line2D([],[], color='None', label=diff_lab)
    handles.append(diff_han)
    labels.append(diff_lab)

    B.pl.legend(handles=handles,labels=labels,
                loc='upper right',fontsize='medium',handlelength=0)
    B.pl.vlines([hmean,hSmean],ymin=0,ymax=[h.A.value,hS.A.value],
                linestyles='--',colors=['#144bcc','#ed1866'])
    
    i+=1
    ax.set_autoscaley_on(True)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'Run {r} {s}')     
#%%    
plot_2d = True

if plot_2d:
    histos_to_plot = {'SHMS $\delta$ vs W': ed_W_histos_wcuts,
                      'SHMS $\delta$ vs Em': ed_Em_histos_wcuts,
                      'HMS $\delta$ vs W': hd_W_histos_wcuts,
                      'HMS $\delta$ vs Em': hd_Em_histos_wcuts} 
    for plot in histos_to_plot:
        fig = B.pl.figure(layout='constrained',figsize=(20,10))
        fig.suptitle(f'{plot}', fontsize=16)
        rows = 2
        cols = 4
        i = 1
        for r in RUN.many:
            ax = plt.subplot(rows,cols,i)
            h = histos_to_plot[plot][r]
            
            h.plot(colormap=cmp['viridis'])
            
            if plot[0] == 'S' and plot[17] == 'W':
                B.pl.vlines(x=0.944,ymin=edelta_ran[r][0],
                            ymax=edelta_ran[r][1],colors='red',
                            linestyles='--')
                if edelta_ran[r][0] < -7:
                    B.pl.hlines(y=-7,xmin=0.7,xmax=1.2,colors='black',
                                linestyles='--')
                    
            if plot[0] == 'S' and plot[17] == 'E':
                B.pl.vlines(x=0.0,ymin=edelta_ran[r][0],
                            ymax=edelta_ran[r][1],colors='red',
                            linestyles='--')
                if edelta_ran[r][0] < -7:
                    B.pl.hlines(y=-7,xmin=-0.1,xmax=0.1,colors='black',
                                linestyles='--')
                    
            if plot[0] == 'H' and plot[16] == 'W':
                B.pl.vlines(x=0.944,ymin=-10,
                            ymax=10,colors='red',
                            linestyles='--')
            
            if plot[0] == 'H' and plot[16] == 'E':
                B.pl.vlines(x=0.0,ymin=-10,
                            ymax=10,colors='red',
                            linestyles='--')
                    
            i+=1
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title(f'Run {r}')     

        
#%% plot histos: Em
Em_means_data = {}
Em_means_SIMC = {}

fig = B.pl.figure(layout='constrained',figsize=(20,10))
fig.suptitle('Missing Energy (with cuts)', fontsize=16)
rows = 2
cols = 4
i = 1

diffs_em = {}
for r in RUN.many:
    s = RUN.setting[r]
    ax = plt.subplot(rows,cols,i)
    h = Em_histos_wcuts[r]
    h.fit(plot_fit=False)
    hmean = h.mean.value
    hsigma = h.sigma.value
    h.fit(hmean-hsigma,hmean+hsigma,plot_fit=False)
    hmean = h.mean.value
    
    hS = SIMC_Em_histos_wcuts[s]
    hS.fit(-0.01,0.0,plot_fit=False)
    hSmean = hS.mean.value
    hSsigma = hS.sigma.value
    hS.fit(hSmean-hSsigma,hSmean+hSsigma,plot_fit=False)
    hSmean = hS.mean.value
    
    h.plot(hatch='......',facecolor='white',edgecolor='#7b8fd4')
    h.plot_fit(color='#144bcc')
    hS.plot(filled=False, color='#ed6792')
    hS.plot_fit(color='#ed1866')
    
    means_diff = hmean - hSmean
    diffs_em[r] = means_diff
    Em_means_data[r] = hmean
    Em_means_SIMC[s] = hSmean
    
    handles, labels = ax.get_legend_handles_labels()

    diff_lab = f'DATA - SIMC = {means_diff:.3e}'
    diff_han = mlines.Line2D([],[], color='None', label=diff_lab)
    handles.append(diff_han)
    labels.append(diff_lab)

    B.pl.legend(handles=handles,labels=labels,
                loc='upper right',fontsize='medium',handlelength=0)
    B.pl.vlines([hmean,hSmean],ymin=0,ymax=[h.A.value,hS.A.value],
                linestyles='--',colors=['#144bcc','#ed1866'])
    i+=1
    ax.set_autoscaley_on(True)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'Run {r} {s}')    
# B.pl.savefig('./offset_studies/Em_allheep-deut18.png')    

#%% plot histos: W (only 1 plot)
    
B.pl.figure()
h = W_all
# h.fit(plot_fit=False)
# hmean = h.mean.value
# hsigma = h.sigma.value
# h.fit(hmean-1.2*hsigma,hmean+1.2*hsigma,plot_fit=False)
# hmean = h.mean.value

hS = SIMC_W_histos_wcuts
# hS.fit(plot_fit=False)
# hSmean = hS.mean.value
# hSsigma = hS.sigma.value
# hS.fit(hSmean-1.2*hSsigma,hSmean+1.2*hSsigma,plot_fit=False)
# hSmean = hS.mean.value

h.plot(hatch='......',facecolor='white',edgecolor='#7b8fd4')
h.plot_fit(color='#144bcc')
hS.plot(filled=False, color='#ed6792')
hS.plot_fit(color='#ed1866')

# means_diff = hmean - hSmean
ax = B.pl.gca()
handles, labels = ax.get_legend_handles_labels()

diff_lab = f'DATA - SIMC = {means_diff:.3e}'
diff_han = mlines.Line2D([],[], color='None', label=diff_lab)
handles.append(diff_han)
labels.append(diff_lab)

B.pl.legend(handles=handles,labels=labels,
            loc='upper right',fontsize='medium',handlelength=0)
B.pl.vlines([hmean,hSmean],ymin=0,ymax=[h.A.value,hS.A.value],
            linestyles='--',colors=['#144bcc','#ed1866'])

ax.set_autoscaley_on(True)

# pptxify(t = f'Invariant Mass for {SIMC.setting}',
#         x = 'W (GeV)', y = '') 

#%% plot histos: Em (only 1 plot)
B.pl.figure()
h = Em_all
h.fit(plot_fit=False)
hmean = h.mean.value
# hsigma = h.sigma.value
# h.fit(hmean-1.2*hsigma,hmean+1.2*hsigma,plot_fit=False)
# hmean = h.mean.value

hS = SIMC_Em_histos_wcuts
hS.fit(plot_fit=False)
hSmean = hS.mean.value
# hSsigma = hS.sigma.value
# hS.fit(hSmean-1.2*hSsigma,hSmean+1.2*hSsigma,plot_fit=False)
# hSmean = hS.mean.value

h.plot(hatch='......',facecolor='white',edgecolor='#7b8fd4')
h.plot_fit(color='#144bcc')
hS.plot(filled=False, color='#ed6792')
hS.plot_fit(color='#ed1866')

means_diff = hmean - hSmean
ax = B.pl.gca()
handles, labels = ax.get_legend_handles_labels()

diff_lab = f'DATA - SIMC = {means_diff:.3e}'
diff_han = mlines.Line2D([],[], color='None', label=diff_lab)
handles.append(diff_han)
labels.append(diff_lab)

B.pl.legend(handles=handles,labels=labels,
            loc='upper right',fontsize='medium',handlelength=0)
B.pl.vlines([hmean,hSmean],ymin=0,ymax=[h.A.value,hS.A.value],
            linestyles='--',colors=['#144bcc','#ed1866'])

ax.set_autoscaley_on(True)

# pptxify(t = f'Missing Energy for run {SIMC.setting}',
#         x = 'Em (GeV)', y = '')

#%%
means_W_pm120 = np.array([1.48769,1.05351,0.97395])
smeans_W_pm120 = np.array([1.45438,1.00059,0.96586])

means_Em_pm120 = np.array([-0.00266,0.01037,0.01341])
smeans_Em_pm120 = np.array([0.01451,0.01180,0.00418])
delt = np.array([-8,-4,0])

delta = np.array([-4,0,4,8,12])
# diffs_sorted = np.array([diffs_em[20841],diffs_em[20846],diffs_em[20851],
#                 diffs_em[20858],diffs_em[20861],diffs_em[20869]])
# wms = np.array([W_means_SIMC['delta_scan_-8'],W_means_SIMC['delta_scan_-4'],
#                 W_means_SIMC['delta_scan_0'],W_means_SIMC['delta_scan_+4'],
#                 W_means_SIMC['delta_scan_+8'],W_means_SIMC['delta_scan_+12']])
# wmd = np.array([W_means_data[20841],W_means_data[20846],
#                 W_means_data[20851],W_means_data[20858],
#                 W_means_data[20861],W_means_data[20840]])
# emms = np.array([Em_means_SIMC['delta_scan_-8'],Em_means_SIMC['delta_scan_-4'],
#                 Em_means_SIMC['delta_scan_0'],Em_means_SIMC['delta_scan_+4'],
#                 Em_means_SIMC['delta_scan_+8'],Em_means_SIMC['delta_scan_+12']])
# emmd = np.array([Em_means_data[20841],Em_means_data[20846],
#                 Em_means_data[20851],Em_means_data[20858],
#                 Em_means_data[20861],Em_means_data[20840]])
# means_DAT = np.array([means_data[20844],means_data[20847],
#                          means_data[20856],means_data[20859],
#                          means_data[20862],means_data[20867]])

wms = np.array([W_means_SIMC['delta_scan_-4'],
                W_means_SIMC['delta_scan_0'],W_means_SIMC['delta_scan_+4'],
                W_means_SIMC['delta_scan_+8'],W_means_SIMC['delta_scan_+12']])
wmd = np.array([W_means_data[20846],
                W_means_data[20851],W_means_data[20858],
                W_means_data[20861],W_means_data[20840]])
emms = np.array([Em_means_SIMC['delta_scan_-4'],
                Em_means_SIMC['delta_scan_0'],Em_means_SIMC['delta_scan_+4'],
                Em_means_SIMC['delta_scan_+8'],Em_means_SIMC['delta_scan_+12']])
emmd = np.array([Em_means_data[20846],
                Em_means_data[20851],Em_means_data[20858],
                Em_means_data[20861],Em_means_data[20840]])


B.pl.figure()
# plt.scatter(delta, diffs_sorted)
plt.scatter(delta,wmd)
plt.scatter(delta,wms)
B.pl.hlines(0.944,xmin=-5,xmax=13,linestyles='--',color='r')

plt.title('W means vs delta central settings: deltaoptim3')

B.pl.figure()
plt.scatter(delta,emmd)
plt.scatter(delta,emms)
B.pl.hlines(0.0,xmin=-5,xmax=13,linestyles='--',color='r')
plt.title('Em means vs delta central settings: deltaoptim3')

# plt.scatter(delt,means_W_pm120)
# plt.scatter(delt,smeans_W_pm120)

# plt.scatter(delt,means_Em_pm120)
# plt.scatter(delt,smeans_Em_pm120)


# plt.title('W means vs delta central settings: pm120')



#fit0 = B.linefit(np.array([-4,0,4]), np.array([diffs_em[20846],diffs_em[20851],diffs_em[20858]]))
# fit1 = B.linefit(delta, diffs_sorted)

# B.pl.text(-8,0.0015,f'slope = {fit1.slope:.4f}$\pm${fit1.sigma_s:.4f}')
# B.pl.text(-8,0.0010,f'offset = {fit1.offset:.4f}$\pm${fit1.sigma_o:.4f}')
#fit = B.polyfit(delta, diffs_sorted, order=2)

# B.pl.savefig('STEP1_DATA-SIMC_vs_delta.png') 

#%% W delta dep correction

W0 = means_data[20868] # delta +12 setting
W_diff = means_DAT - W0

B.pl.figure()
plt.scatter(delta,W_diff)
fit = B.polyfit(delta, W_diff, order=4)

def corr_term(delta):
    par0 = -0.00252
    par1 = 0.00027
    par2 = -0.00003
    par3 = -0.000001
    par4 = 0.0000003
    
    res = par0 + par1*delta +\
            par2*(delta**2) + par3*(delta**3) +\
                par4*(delta**4)
    return res

f_corr = corr_term(delta)
W_corr = means_DAT - f_corr

B.pl.figure()
plt.scatter(delta, W_corr)
plt.scatter(delta,means_SIM)

#%%
# exclude delta +12
delta = np.array([-8,-4,0,4,8])
means_DAT = np.array([means_data[20841],means_data[20846],
                         means_data[20851],means_data[20858],
                         means_data[20861]])

W0 = means_data[20861] # delta +12 setting
W_diff = means_DAT - W0

B.pl.figure()
plt.scatter(delta,W_diff)
fit = B.polyfit(delta, W_diff, order=2)

def corr_term(delta):
    par0 = fit.par[0]
    par1 = fit.par[1]
    par2 = fit.par[2]
    # par3 = fit.par[3]
    # par4 = 0.0000003
    
    res = par0 + par1*delta +\
            par2*(delta**2) 
                # par3*(delta**3) 
            #+ par4*(delta**4)
    return res

f_corr = corr_term(delta)
W_corr = means_DAT - f_corr

B.pl.figure()
plt.scatter(delta, W_corr)
delta = np.array([-8,-4,0,4,8,12])
plt.scatter(delta,means_SIM)
#%% calculate pf correction from Em difference
Mp = 0.938272
pf_cent = np.array([3.499,3.145,2.783,2.417,2.048,1.664])

dEmdpf = -pf_cent/(np.sqrt(Mp*Mp + pf_cent*pf_cent))

diffem = np.array([-0.0089,-0.0026,-0.0024,-0.0034,-0.0006,0.0007])

pf_corr = diffem/dEmdpf

print(pf_corr,'\n',(1-pf_corr/pf_cent))
