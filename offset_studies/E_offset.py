#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:38:18 2024

@author: gvill
"""

import numpy as np
import LT.box as B
from tabulate import tabulate

import deut_analysis_tools as dat

#%%

fname = 'Ebeam_12GeV_meas.data'
f = B.get_file(fname)

E_beam = f['E_beam']
dE_beam = f['dE_beam']
sync_rad = f['sync_rad']
hallcp = f['HALLC:p']/1e03      # in GeV

corr = ((E_beam+dE_beam) - 0.5*sync_rad)/hallcp
corr1 = ((E_beam) - 0.5*sync_rad)/hallcp

tab = {'Date':['04-30-2018','10-04-2018','03-13-2019','02-20-2020',
                '08-19-2020','06-12-2022','09-19-2023'],
        'E_beam':E_beam,'dE_beam':dE_beam,'sync_rad':sync_rad,'HALLC:p':hallcp,
        'corr':corr,'corr1':corr1}

print(tabulate(tab,headers=tab.keys(),tablefmt=''))
print('')
print('corr = ((E_beam+dE_beam) - 0.5*sync_rad)/hallcp')
print('corr mean = ', corr.mean())
print('')
print('corr1 = ((E_beam) - 0.5*sync_rad)/hallcp')
print('corr1 mean = ', corr1.mean())

#%%
# getting HALLC:p value from E tree

e_runs = [20845,20847,20851,20856,20859,20862,
               20867,20872,20878,20922,20956,
                   21005,21048,21054,21061,21070,
                       21090]

e_data = dat.DEUT_DATA(Runs=e_runs,branches_sel=['HALLC_p'],trees_sel=['E'])

hallcp_means = []
for r in e_runs:
    m = e_data.Branches[r]['HALLC_p'].mean()
    hallcp_means.append(m)

hallcp_means = np.array(hallcp_means)
hallcp_value = hallcp_means.mean()

print('HALLC_p = ', hallcp_value)
