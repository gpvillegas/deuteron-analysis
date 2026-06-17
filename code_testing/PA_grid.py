#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 12:14:54 2026

@author: gvill
"""

import data_init as D
import numpy as np

T_branches = ['H.dc.x_fp','H.dc.xp_fp','H.dc.y_fp','H.dc.yp_fp',
              'H.gtr.dp']

R = 20851

ROOTfiles_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

RUN = D.DATA_INIT(data_type='deut23_data',run=R,
                         select_branches={'T':T_branches},
                         select_trees=['T'],
                         ROOTfiles_path=ROOTfiles_DIR) 

xfp = RUN.Branches['H.dc.x_fp']
yfp = RUN.Branches['H.dc.y_fp']
xpfp = RUN.Branches['H.dc.xp_fp']
ypfp = RUN.Branches['H.dc.yp_fp']
dp = RUN.Branches['H.gtr.dp']
#%%

with open('20851_hmsKIN.txt','w') as f:
    for i in range(xfp.size):
        if dp[i] > 10 or dp[i] < -10 or dp[i] != dp[i] or\
            xpfp[i] > 0.08 or xpfp[i] < -0.08 or ypfp[i] > 0.04 or\
                ypfp[i] < -0.04: continue #get rid of nan
        f.write(f'{xfp[i]:0.6f} {yfp[i]:0.6f} {xpfp[i]:0.6f} {ypfp[i]:0.6f} {dp[i]:0.6f}\n')
        
    
