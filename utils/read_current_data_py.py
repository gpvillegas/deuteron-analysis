#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 15:08:15 2026

read a asci simc output file

@author: boeglinw
"""

import numpy as np
import LT.box as B

d = './d2_pm800_Paris_fsi_rad.data'


l = open(d).readlines()



#%%
n_col = int(l[2].strip()) 

names = {}
for i in range(n_col+1):
    names[l[2+i].strip()] = i-1
#%%

i_start = 2+n_col+1
data = np.loadtxt(l[i_start:])

