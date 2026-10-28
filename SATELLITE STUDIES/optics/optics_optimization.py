#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:19:39 2024

@author: gvill
"""

import numpy as np
import matplotlib.pyplot as plt
import LT.box as B

import deut_analysis_tools as dat
import data_init as D

#%%
"""
xptar, yptar plots
"""

br = ['P.gtr.th','P.gtr.ph','H.gtr.th','H.gtr.ph','P.kin.primary.W']
brSIMC = ['e_xptar','e_yptar','h_xptar','h_yptar','W','Weight','Normfac']

heep_coin = D.DATA_INIT(data_type='deut23_data',kin_study='heep_coin',
                        select_branches=br)
heepSIMC = D.DATA_INIT(data_type='SIMC',kin_study='heep_coin',
                       select_branches=brSIMC)
# heep_coin = dat.DEUT_DATA(kin_study='heep_coin',branches_sel=br)
# heepSIMC = dat.DEUT_SIMC(kin_study='heep_coin',branches_sel=brSIMC)
#heep_singles = dat.DEUT_DATA(kin_study='heep_singles',branches_sel=br)
#%%
HMSang_acc = dat.make_plots(['H.gtr.th','H.gtr.ph'],heep_coin.Branches,
                            heep_coin.many,apply_cuts=False)
SIMC_HMSang_acc = dat.make_plots_SIMC(['h_xptar','h_yptar'],
                                      heepSIMC.Branches['delta_scan_+12'],
                                      apply_cuts=False)

SHMSang_acc = dat.make_plots(['P.gtr.th','P.gtr.ph'],heep_coin.Branches,
                            heep_coin.many,apply_cuts=False)

HMSWxptar = dat.make_plots(['H.gtr.th','P.kin.primary.W'],heep_coin.Branches,
                            heep_coin.many,apply_cuts=False)

HMSWyptar = dat.make_plots(['H.gtr.ph','P.kin.primary.W'],heep_coin.Branches,
                            heep_coin.many,apply_cuts=False)

SHMSWxptar = dat.make_plots(['P.gtr.th','P.kin.primary.W'],heep_coin.Branches,
                            heep_coin.many,apply_cuts=False)

SHMSWyptar = dat.make_plots(['P.gtr.ph','P.kin.primary.W'],heep_coin.Branches,
                            heep_coin.many,apply_cuts=False)
#%%

HMSang_acc.make_1Dhistos([(-0.1,0.1),1000],fig_title='HMS')

SHMSang_acc.make_1Dhistos([(-0.1,0.1),1000],fig_title='SHMS')

HMSang_acc.make_2Dhistos([(-0.1,0.1),100,(-0.06,0.06),100], 
                         fig_title='HMS Xptar vs. Yptar')

SHMSang_acc.make_2Dhistos([(-0.05,0.05),100,(-0.02,0.02),100], 
                          fig_title='SHMS Xptar vs. Yptar')

HMSWxptar.make_2Dhistos([(-0.1,0.1),100,(0.5,1.5),200],
                        fig_title='HMS W vs. Xptar')

HMSWyptar.make_2Dhistos([(-0.06,0.06),100,(0.5,1.5),200],
                        fig_title='HMS W vs. Yptar')

SHMSWxptar.make_2Dhistos([(-0.05,0.05),100,(0.5,1.5),200],
                         fig_title='SHMS W vs. Xptar')

SHMSWyptar.make_2Dhistos([(-0.02,0.02),100,(0.5,1.5),200],
                         fig_title='SHMS W vs. Yptar')
#%%
"""
focal plane plots x_fp, y_fp, xp_fp, yp_fp
"""

br = dat.branches_types('optics')

#brSIMC = ['e_xptar','e_yptar','h_xptar','h_yptar','W','Weight','Normfac']

heep_coin = dat.DEUT_DATA(kin_study='heep_coin',branches_sel=br)
#heepSIMC = dat.DEUT_SIMC(kin_study='heep_coin',branches_sel=brSIMC)
#heep_singles = dat.DEUT_DATA(kin_study='heep_singles',branches_sel=br)

p = dat.make_plots(['P.dc.y_fp','P.dc.x_fp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.dc.y_fp','H.dc.x_fp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-35,35),200,(-40,25),200],fig_title='SHMS x_fp:y_fp',
                save_as='SHMS_XfpYfp.pdf')
h.make_2Dhistos([(-25,25),200,(-40,35),200],fig_title='HMS x_fp:y_fp',
                save_as='HMS_XfpYfp.pdf')

p = dat.make_plots(['P.dc.x_fp','P.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.dc.x_fp','H.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-40,25),200,(-15,25),200],fig_title='SHMS dp:x_fp')
h.make_2Dhistos([(-40,35),200,(-15,15),200],fig_title='HMS dp:x_fp')

p = dat.make_plots(['P.dc.y_fp','P.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.dc.y_fp','H.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-35,35),200,(-15,25),200],fig_title='SHMS dp:y_fp')
h.make_2Dhistos([(-25,25),200,(-15,15),200],fig_title='HMS dp:y_fp')

p = dat.make_plots(['P.dc.xp_fp','P.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.dc.xp_fp','H.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-0.1,0.1),200,(-15,25),200],fig_title='SHMS dp:xp_fp')
h.make_2Dhistos([(-0.1,0.1),200,(-15,15),200],fig_title='HMS dp:xp_fp')

p = dat.make_plots(['P.dc.yp_fp','P.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.dc.yp_fp','H.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-0.05,0.05),200,(-15,25),200],fig_title='SHMS dp:yp_fp')
h.make_2Dhistos([(-0.05,0.05),200,(-15,15),200],fig_title='HMS dp:yp_fp')

p = dat.make_plots(['P.gtr.th','P.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.gtr.th','H.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-0.05,0.05),200,(-15,25),200],fig_title='SHMS dp:xptar')
h.make_2Dhistos([(-0.1,0.1),200,(-15,15),200],fig_title='HMS dp:xptar')

p = dat.make_plots(['P.gtr.ph','P.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.gtr.ph','H.gtr.dp'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-0.05,0.05),200,(-15,25),200],fig_title='SHMS dp:yptar')
h.make_2Dhistos([(-0.05,0.05),200,(-15,15),200],fig_title='HMS dp:yptar')

p = dat.make_plots(['P.extcor.xsieve','P.extcor.ysieve'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)
h = dat.make_plots(['H.extcor.xsieve','H.extcor.ysieve'], heep_coin.Branches, 
                   heep_coin.Runs, apply_cuts=False)

# p.make_1Dhistos()
p.make_2Dhistos([(-10,10),200,(-15,25),200],fig_title='SHMS xsieve:ysieve')
h.make_2Dhistos([(-10,10),200,(-15,15),200],fig_title='HMS xsieve:ysieve')

#%%

# Declare constants
Mp = 0.938272       # Proton mass
Eb = 10.542         # Beam energy
Th_p = 0.           # Measured proton angle

cosTh_p = np.cos(Th_p)
sinTh_p = np.sin(Th_p)

# Calculated momentum P(Eb,Th_p) 
Pcalc = (2.*Mp*Eb*(Eb+Mp)*cosTh_p)/\
            (Mp*Mp + 2.*Mp*Eb + Eb*Eb*sinTh_p*sinTh_p)  
            
#%%

run = dat.DEUT_DATA(Runs = 20851, branches_sel=dat.branches_types('kins'))          

Q2 = run.Branches[20851]['P.kin.primary.Q2']
th_e = run.Branches[20851]['P.kin.primary.scat_ang_rad']
E1 = 10.542
E2 = 8.55
Q2_calc = 4.*E1*E2*np.sin(th_e/2.)*np.sin(th_e/2.)

run.Branches[20851].update({'Q2_calc':Q2_calc})

Q2_hist = dat.make_plots(['P.kin.primary.Q2','Q2_calc'],run.Branches,
                         run.Runs)
Q2_hist_nocuts = dat.make_plots(['P.kin.primary.Q2','Q2_calc'],run.Branches,
                                run.Runs,apply_cuts=False)

Q2_hist.make_1Dhistos([(2.,5.),200])
Q2_hist_nocuts.make_1Dhistos([(2.,5.),200])

B.pl.figure()
diff_hist_nocuts = B.histo(Q2_calc-Q2,bins=200,range=(-1,1))
diff_hist_nocuts.plot()

Q2_cut = Q2_hist.v['P.kin.primary.Q2'][20851]
Q2_calc_cut = Q2_hist.v['Q2_calc'][20851]

diff_hist = B.histo(Q2_calc_cut-Q2_cut,bins=200,range=(-1,1))
diff_hist.plot(color='r',filled=False)
