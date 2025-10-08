#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:39:09 2024

@author: gvill
"""
import LT.box as B
import matplotlib.lines as mlines


# import deut_analysis_tools as dat
import database_operations as db
import data_init as D
import cut_handler as C

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
    
#%% load selection of branches for this specific analysis
deut_db = 'deuteron_db.db' 

kin_var = ['P.kin.primary.W','H.kin.secondary.emiss_nuc']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

br_sel = kin_var + acc_var + calPID_var

br_sel_SIMC = ['e_delta','h_delta','Weight','Normfac','W','Em','e_xptar',
               'e_yptar','h_xptar','h_yptar']

#%% load root files
run_num = 20871
sett = 'pm_120'

# heep = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
#                     select_branches=br_sel)

run = D.DATA_INIT(data_type='deut23_data', run=run_num, 
                    select_branches={'T':br_sel})

# SIMC files

# Sheep = D.DATA_INIT(data_type='SIMC', kin_study='heep_coin')

simc = D.DATA_INIT(data_type='SIMC', setting=sett, 
                   select_branches={'SNT':br_sel_SIMC}, simc_type='jmlfsi_rad')

#%% declare wanted variables

W = run.Branches['P.kin.primary.W']
Em = run.Branches['H.kin.secondary.emiss_nuc']

SW = simc.Branches['W']
SEm = simc.Branches['Em']

#%% Define cuts for data

cuts_list = C.acceptance_cuts

cuts_to_apply = []
for cut in cuts_list:
    br = run.Branches[C.HCANA_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply.append(cut_array)

all_cuts = cuts_to_apply[0]
for arr in cuts_to_apply:    
    all_cuts = all_cuts & arr

# Define cuts for SIMC

cuts_to_apply_simc = []
for cut in cuts_list:
    cut.init()
    br = simc.Branches[C.SIMC_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply_simc.append(cut_array)

all_cuts_simc = cuts_to_apply_simc[0]
for arr in cuts_to_apply_simc:    
    all_cuts_simc = all_cuts_simc & arr
#%% apply cuts to variables

W_c = W[all_cuts]
Em_c = Em[all_cuts]
    
SW_c = SW[all_cuts_simc]
SEm_c = SEm[all_cuts_simc]

#%% calculate normalization factors for data runs and get weights for SIMC

ENORM = D.get_eff_norm(run_num)

q = D.get_charge_norm(run_num)
NORM = ENORM*(1/q)
  
weights = D.calc_weights(simc.Branches)
WEIGHTS = weights[all_cuts_simc]

#%% make histograms 

W_histos_wcuts = NORM*B.histo(W_c,range = (0.5,2),bins = 50)
Em_histos_wcuts = NORM*B.histo(Em_c,range = (-0.025,0.1),bins = 50)

SIMC_W_histos_wcuts = B.histo(SW_c,range = (0.5,2),bins = 50,
                                 weights=WEIGHTS, calc_w2=True)
SIMC_Em_histos_wcuts = B.histo(SEm_c,range = (-0.025,0.1),bins = 50,
                                  weights=WEIGHTS, calc_w2=True)

#%% plot histos: W
    
B.pl.figure()
h = W_histos_wcuts
h.fit(plot_fit=False)
hmean = h.mean.value
hsigma = h.sigma.value
h.fit(hmean-hsigma,hmean+hsigma,plot_fit=False)
hmean = h.mean.value

hS = SIMC_W_histos_wcuts
hS.fit(plot_fit=False)
hSmean = hS.mean.value
hSsigma = hS.sigma.value
hS.fit(hSmean-hSsigma,hSmean+hSsigma,plot_fit=False)
hSmean = hS.mean.value

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

#B.pl.title(f'Run {r}')

pptxify(t = f'Invariant Mass for run {run_num} {sett}',
        x = 'W (GeV)', y = '')

# save_where = '/home/gvill/Desktop/plots for pptx/'
# fname = save_where + 'W_20871_pm120.png'
# B.pl.savefig(fname)   
#%% plot histos: Em
B.pl.figure()
h = Em_histos_wcuts
h.fit(plot_fit=False)
hmean = h.mean.value
hsigma = h.sigma.value
h.fit(hmean-hsigma,hmean+hsigma,plot_fit=False)
hmean = h.mean.value

hS = SIMC_Em_histos_wcuts
hS.fit(0,0.025,plot_fit=False)
hSmean = hS.mean.value
hSsigma = hS.sigma.value
hS.fit(hSmean-hSsigma,hSmean+hSsigma,plot_fit=False)
hSmean = hS.mean.value

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

#B.pl.title(f'Run {r}')

pptxify(t = f'Missing Energy for run {run_num} {sett}',
        x = 'Em (GeV)', y = '')

# save_where = '/home/gvill/Desktop/plots for pptx/'
# fname = save_where + 'Em_20871_pm120.png'
# B.pl.savefig(fname)  
