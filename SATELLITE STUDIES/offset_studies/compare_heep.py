#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:39:09 2024

@author: gvill
"""
import numpy as np
import LT.box as B
import matplotlib.pyplot as plt

import deut_analysis_tools as dat

'''
Comparison plots
'''
#heep delta 0
heep0 = dat.DEUT_DATA(branches_sel=dat.branches_types('kins'),
                      Runs=20851)
heep0_SIMC = dat.DEUT_SIMC(kin_study='heep_coin',setting='delta_scan_0')
heep0_SIMC = heep0_SIMC.Branches['delta_scan_0']

pmS = dat.make_plots_SIMC('Pm', heep0_SIMC)
pm = dat.make_plots('H.kin.secondary.pmiss', heep0.Branches, heep0.Runs,)

pm.make_1Dhistos([(-0.1,2.5),200],apply_norm=True,SIMC=pmS,
                 fig_title='Pmiss')

EmS = dat.make_plots_SIMC('Em', heep0_SIMC)
Em = dat.make_plots('H.kin.secondary.emiss', heep0.Branches, heep0.Runs)
                    

Em.make_1Dhistos([(-0.1,1.0),200],apply_norm=True,SIMC=EmS,
                 fig_title='Emiss')

WS = dat.make_plots_SIMC('W', heep0_SIMC)
W = dat.make_plots('P.kin.primary.W', heep0.Branches, heep0.Runs)
                   

W.make_1Dhistos([(0.5,2.0),200],apply_norm=True,SIMC=WS,
                fig_title='W')
'''
'''
#heep delta +8
heep8 = dat.DEUT_DATA(branches_sel=dat.branches_types('kins'),
                      Runs=20861)
heep8_SIMC = dat.DEUT_SIMC(kin_study='heep_coin',setting='delta_scan_+8')
heep8_SIMC = heep8_SIMC.Branches['delta_scan_+8']

pmS = dat.make_plots_SIMC('Pm', heep8_SIMC)
pm = dat.make_plots('H.kin.secondary.pmiss', heep8.Branches, heep8.Runs)
                    

pm.make_1Dhistos([(-0.1,2.5),200],apply_norm=True,SIMC=pmS,
                 fig_title='Pmiss')

EmS = dat.make_plots_SIMC('Em', heep8_SIMC)
Em = dat.make_plots('H.kin.secondary.emiss', heep8.Branches, heep8.Runs)
                   

Em.make_1Dhistos([(-0.1,1.0),200],apply_norm=True,SIMC=EmS,
                 fig_title='Emiss')

WS = dat.make_plots_SIMC('W', heep8_SIMC)
W = dat.make_plots('P.kin.primary.W', heep8.Branches, heep8.Runs)

W.make_1Dhistos([(0.5,2.0),200],apply_norm=True,SIMC=WS,
                fig_title='W')

#heep delta -8
heepm8 = dat.DEUT_DATA(branches_sel=dat.branches_types('kins'),
                      Runs=20841)
heepm8_SIMC = dat.DEUT_SIMC(kin_study='heep_coin',setting='delta_scan_-8')
heepm8_SIMC = heepm8_SIMC.Branches['delta_scan_-8']

pmS = dat.make_plots_SIMC('Pm', heepm8_SIMC)
pm = dat.make_plots('H.kin.secondary.pmiss', heepm8.Branches, heepm8.Runs)
                    

pm.make_1Dhistos([(-0.1,2.5),200],apply_norm=True,SIMC=pmS,
                 fig_title='Pmiss')

EmS = dat.make_plots_SIMC('Em', heepm8_SIMC)
Em = dat.make_plots('H.kin.secondary.emiss', heepm8.Branches, heepm8.Runs)
                    

Em.make_1Dhistos([(-0.1,1.0),200],apply_norm=True,SIMC=EmS,
                 fig_title='Emiss')

WS = dat.make_plots_SIMC('W', heepm8_SIMC)
W = dat.make_plots('P.kin.primary.W', heepm8.Branches, heepm8.Runs)
                   

W.make_1Dhistos([(0.5,2.0),200],apply_norm=True,SIMC=WS,
                fig_title='W')

