#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:10:46 2024

@author: gvill
"""

import numpy as np
import matplotlib.pyplot as plt
import LT.box as B

import deut_analysis_tools as dat

#%%

b120 = dat.pm120.Branches
r120 = dat.pm120.Runs
bS120 = dat.pm120_SIMC

heep = dat.DEUT_DATA(kin_study='heep_coin',branches_sel=['kins'])
heep_SIMC = dat.DEUT_SIMC(kin_study='heep_coin')

#Pmiss
pmS_pm120 = dat.make_plots_SIMC('Pm', bS120)
pm_pm120 = dat.make_plots('H.kin.secondary.pmiss', b120, r120)

pm_pm120.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,SIMC=pmS_pm120)

#Emiss
pmS_pm120 = dat.make_plots_SIMC('Em', bS120)
pm_pm120 = dat.make_plots('H.kin.secondary.emiss_nuc', b120, r120)

pm_pm120.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,SIMC=pmS_pm120)

#BCM4A_current
# curr = {}
# for run in r120:
#     dat.pm120.Trees[run].select_tree('TSP')
#     curr[run] = dat.pm120.Trees[run].get_branches('P.BCM4A.scalerCurrent')
    
#     b120[run].update(curr[run])

curr = {}
for run in r120:
    curr[run] = dat.get_list(\
        dat.db.retrieve(dat.deutDB_name,'BCM4A_current','RUN_LIST',
                        where = f"run = \'{run}\'"))[0]
        
heep_curr = {}
for run in heep.Runs:
    heep_curr[run] = dat.get_list(\
        dat.db.retrieve(dat.deutDB_name,'BCM4A_current','RUN_LIST',
                        where = f"run = \'{run}\'"))[0]        

# #BCM4A_charge
# BCM4A_q = {}
# for run in r120:
#     dat.pm120.Trees[run].select_tree('TSP')
#     q = dat.pm120.Trees[run].get_branches('P.BCM4A.scalerCharge')['P.BCM4A.scalerCharge']
#     histo_q = B.histo(q,(-0.1,2500),100)
#     BCM4A_q[run] = histo_q.sum()
    
#     b120[run].update(BCM4A_q[run])

       
chargedb = {}
for run in r120:
    chargedb[run] = dat.get_list(\
        dat.db.retrieve(dat.deutDB_name,'BCM4A_charge','RUN_LIST',
                        where = f"run = \'{run}\'"))[0]  

heep_charge = {}
for run in heep.Runs:
    heep_charge[run] = dat.get_list(\
        dat.db.retrieve(dat.deutDB_name,'BCM4A_charge','RUN_LIST',
                        where = f"run = \'{run}\'"))[0]          

### heep yield

# W vs Emiss
heepSIMC_emW= {}
heepSIMC_histos = {}
for sett in heep_SIMC.setting:
    heepSIMC_emW[sett] = dat.make_plots_SIMC(['Em','W'], 
                                             heep_SIMC.Branches[sett])
    
    heepSIMC_histos[sett] = heepSIMC_emW[sett].make_2Dhistos(
                                [(-0.1,0.6),100,(0.8,1.0),100],
                                 fig_title=f'W(Emiss) SIMC {sett}',
                                 apply_weights=True,plot=True)
    
heep_emW = dat.make_plots(['H.kin.secondary.emiss','P.kin.primary.W'], 
                          heep.Branches, heep.Runs)

heep_histos = heep_emW.make_2Dhistos([(-0.1,0.6),100,(0.8,1.0),100],
                                     fig_title='W(Emiss)',
                                     apply_norm=False,plot=True)



#Emiss vs Pmiss
empmS_pm120 = dat.make_plots_SIMC(['Pm','Em'],bS120)
empm_pm120 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.emiss_nuc'],b120,20871)

h = empm_pm120.make_2Dhistos([(-0.1,1.0),200,(-0.1,1.0),200],
                             fig_title='Emiss vs. Pmiss',
                             apply_norm=False,plot=True)

hS = empmS_pm120.make_2Dhistos([(-0.1,1.0),100,(-0.1,1.0),100],
                               fig_title='Emiss vs. Pmiss SIMC',
                               apply_weights=True,plot=False)


#%%
# xbin number
bin_range = np.arange(10,28)


Em_min = -0.05
Em_max_a = np.array([0.05, 0.07, 0.09, 0.11])
Em_max =  0.09

B.pl.close('all')
plot_all = True
plot_separate = False

#%% perform projections and integrate
def Em_integration(bin_range, Emin = Em_min, Emax = Em_max, RUNS=r120,
                   plot_all = False):
    counts_simc = []
    for nx in bin_range[:]:
        hsimc = hS.project_y(bins = [nx])
        
        counts_simc.append(hsimc.sum(Emin, Emax))
    
    counts_simc = np.array(counts_simc)           
    
    counts_exp_RUNS = {}
    for run in RUNS:
        counts_exp = []
        for nx in bin_range[:]:
            hexp = h[run].project_y(bins = [nx])
            
            counts_exp.append(hexp.sum(Emin, Emax))
            
            if plot_all:
                B.pl.figure()
                hexp.plot(filled = False)
                hsimc.plot(filled = False)           
        
        counts_exp = np.array(counts_exp)
        counts_exp_RUNS[run] = counts_exp
    
    pm_values = {}
    R = {}
    dR = {}
    for run in RUNS:                
        # calculate ratios    
        pm_values[run] = h[run].x_bin_center[bin_range]
        
        R[run] = counts_exp_RUNS[run][:,0]/counts_simc[:,0]
        
        # assuming the errors in the exp. histo are correct
        dR[run] = counts_exp_RUNS[run][:,1]/counts_simc[:,0]
    
    return pm_values, R, dR, counts_exp_RUNS

#%%

save = False

yield_RUNs = {}
dyield_RUNs = {}

for run in [20871]:
    B.pl.figure()
    for Emax in Em_max_a:
        pm_values, R, dR, counts =\
            Em_integration(bin_range, Em_min, Emax, RUNS=[20871], 
                           plot_all = False)

        # B.plot_exp(pm_values[run], R[run], dR[run], 
        #            label = f'Em_max = {Emax:.2f} (GeV)')
        
        YIELD = counts[run][:,0]
        dY = counts[run][:,1]
        yield_RUNs[run] = max(YIELD)
        dyield_RUNs[run] = max(dY)
        B.plot_exp(pm_values[run], YIELD, dY, 
                   label = f'Em_max = {Emax:.2f} (GeV)')
    
    B.pl.legend()
    B.pl.xlabel('Pm (GeV/c)')
    B.pl.ylabel('Y_exp/Y_simc')
    B.pl.title(f'Yield Ratios for {run}')

    if save:
        B.pl.savefig(f'Yield_ratios_{run}.pdf')
        
#%%

yield_RUNs = []
dyield_RUNs = []

for run in r120:
    pm_values, R, dR, counts =\
        Em_integration(bin_range, Em_min, Em_max, RUNS=r120, 
                       plot_all = False)
    
    YIELD = counts[run][:,0]
    dY = counts[run][:,1]
    yield_RUNs.append(max(YIELD))
    dyield_RUNs.append(max(dY))
    
c = np.array(list(curr.values()))

yi = np.array([17.76,33.11,61.21,50.55,68.97,65.62,74.01,76.50])
yi_rel = yi/min(yi)

dyi = np.array([0.6,0.53,0.5,0.71,0.75,0.73,0.78,0.79])
dyi_rel = dyi/yi + 0.53/17.76

c = np.sort(np.array([35.59,9.92,19.61,29.29,38.98,47.42,56.34,65.36]))

B.plot_exp(c, yi_rel, dyi_rel, label = '$LD_2$ target')

B.plot_line(c, yi_rel)

B.pl.legend()
B.pl.xlabel('Avg Current ($\mu$A)')
B.pl.ylabel('Relative Y_exp/mC')
B.pl.title('Exp. Yield vs. Current')

if save:
    B.pl.savefig(f'Yield_ratios_{run}.pdf')
    
#%%
x1 = 0.01
x2 = 0.15
y1 = -0.01
y2 = 0.03


YIELD = {}
for run in r120:
    htemp = h[run]
    rc = htemp.rect_cut(x1,x2,y1,y2)
    
    YIELD[run] = htemp.sum(rc)
    
qy = []
for run in r120:
    qy.append([chargedb[run],YIELD[run][0],YIELD[run][1]]) 
    
qy = np.array(qy)

q = qy[:,0]
qnorm_yield = qy[:,1]/qy[:,0]
dyield = qy[:,2]/qy[:,0]

B.plot_exp(q, qnorm_yield, dyield, label = '$LD_2$ target')

B.pl.legend()
B.pl.xlabel('Accumulated Charge (mC)')
B.pl.ylabel('Relative Y_exp/mC')
B.pl.title('Exp. Yield vs. Charge')
 
    
ay = []
for run in r120:
    ay.append([curr[run],YIELD[run][0]])

ay = np.array(ay)

#%%
x1 = -0.01
x2 = 0.03
y1 = 0.9
y2 = 1.0

YIELD = {}
for run in heep.Runs:
    htemp = heep_histos[run]
    rc = htemp.rect_cut(x1,x2,y1,y2)
        
    YIELD[run] = htemp.sum(rc)
        
    
qy = []
for run in heep.Runs:
    qy.append([heep_charge[run],YIELD[run][0],YIELD[run][1]]) 
    
qy = np.array(qy)

q = qy[:,0]
qnorm_yield = qy[:,1]/qy[:,0]
dyield = qy[:,2]/qy[:,0]

B.plot_exp(q, qnorm_yield, dyield, label = '$LD_2$ target')

B.pl.legend()
B.pl.xlabel('Accumulated Charge (mC)')
B.pl.ylabel('Relative Y_exp/mC')
B.pl.title('Exp. Yield vs. Charge')

## all
ay = []
for run in heep.Runs:
    ay.append([heep_curr[run],YIELD[run][0],YIELD[run][1]])

ay = np.array(ay)
        
I = ay[:,0]
qnorm_yield = ay[:,1]/ay[:,0]
dyield = ay[:,2]/ay[:,0]
    
B.plot_exp(I, qnorm_yield/14095, dyield/14095, label = '$LH_2$ target')
  
B.pl.legend()
B.pl.xlabel('Avg Current ($\mu$A)')
B.pl.ylabel('$Y_{Norm}(I) / Y_{Norm}(I_{min})$')
B.pl.title('Rel.Exp. Yield vs. Current for heep_coin') 
 
## delta +12    
ay = []
for run in [20840,20868,20869]:
    ay.append([heep_curr[run],YIELD[run][0],YIELD[run][1]])

ay = np.array(ay)
        
I = ay[:,0]
qnorm_yield = ay[:,1]/ay[:,0]
dyield = ay[:,2]/ay[:,0]
    
B.plot_exp(I, qnorm_yield/14095, dyield/14095, label = '$LH_2$ target')
  
B.pl.legend()
B.pl.xlabel('Avg Current ($\mu$A)')
B.pl.ylabel('$Y_{Norm}(I) / Y_{Norm}(I_{min})$')
B.pl.title('Rel.Exp. Yield vs. Current for heep_coin delta_scan_+12') 

## delta +8/-8    
ay = []
for run in [20841,20861]:
    ay.append([heep_curr[run],YIELD[run][0],YIELD[run][1]])

ay = np.array(ay)
        
I = ay[:,0]
qnorm_yield = ay[:,1]/ay[:,0]
dyield = ay[:,2]/ay[:,0]
    
B.plot_exp(I, qnorm_yield/3245, dyield/3245, label = '$LH_2$ target')
  
B.pl.legend()
B.pl.xlabel('Avg Current ($\mu$A)')
B.pl.ylabel('$Y_{Norm}(I) / Y_{Norm}(I_{min})$')
B.pl.title('Rel.Exp. Yield vs. Current for heep_coin delta_scan_$\pm$8')

## delta +4/-4    
ay = []
for run in [20846,20858]:
    ay.append([heep_curr[run],YIELD[run][0],YIELD[run][1]])

ay = np.array(ay)
        
I = ay[:,0]
qnorm_yield = ay[:,1]/ay[:,0]
dyield = ay[:,2]/ay[:,0]
    
B.plot_exp(I, qnorm_yield/2083, dyield/2083, label = '$LH_2$ target')
  
B.pl.legend()
B.pl.xlabel('Avg Current ($\mu$A)')
B.pl.ylabel('$Y_{Norm}(I) / Y_{Norm}(I_{min})$')
B.pl.title('Rel.Exp. Yield vs. Current for heep_coin delta_scan_$\pm$4')        
        
            