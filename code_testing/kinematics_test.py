#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 11:36:11 2024

@author: boeglinw
"""

import numpy as np
import LT.box as B
import root_util as RU


#%% masses

MP = 938.27208*1e-3
MN = 939.565421*1e-3
MD = 1875.612928*1e-3


#%%
file = '/home/gvill/deuteron/ROOTfiles/prod/deut_replay_prod_20871_-1.root'

rf = RU.root_tree(file, trees = ['T'])
rf.load_trees()

# print the list of branches
rf.list_branches('T')

#%% load the interesting branches
branch_list = ['P.kin.primary.W','P.gtr.p', 'P.kin.primary.nu', 'P.kin.primary.Q2',
               'H.gtr.p','H.kin.secondary.pmiss']

rf.get_branches(branch_list)

#%%
P_f = rf.tree_data['H.gtr.p']
sel_P_f = P_f<5.0

Pm = rf.tree_data['H.kin.secondary.pmiss']
#Em = rf.tree_data['Em']
sel_Pm = Pm < 1.0

nu = rf.tree_data['P.kin.primary.nu']
sel_nu = nu < 11.0

Q2 = rf.tree_data['P.kin.primary.Q2']
sel_Q2 = Q2 < 10.0

is_ok = sel_P_f & sel_Pm & sel_nu & sel_Q2
#missing momentum cut
Pm > 0.0
Pm < 0.5
set1= Pm > 0.0 & Pm < 0.5
set2= Pm > 0.5 & Pm < 1.0
set3= Pm > 1.0 & Pm < 1.5
set4= Pm > 1.5 & Pm < 2.0

cuts = [set1, set2, set3, set4]

#%% calc kin as in averaged kin
E_tot = nu[is_ok] + MD  # total energy at targe

# total proton energy
Ep_f = np.sqrt(P_f[is_ok]**2 + MP**2)


# total energy of recoiling particle 
Er_f = E_tot - Ep_f

# invariant mass of recoiling particle
M_r2 = Er_f**2 - Pm[is_ok]**2
sel_Mr2 = M_r2 > 0.0
M_r = np.sqrt(M_r2[sel_Mr2])

h_Mr = B.histo(M_r, range = (0.5,1.5), bins = 200, 
               xlabel = 'M_recoil', 
               ylabel = 'Counts', 
               title = 'Inv. mass of recoiling particle')


h_Mr.plot(filled = False)

for sets in cuts:
    E_tot = nu[is_ok & sets] + MD  # total energy at targe

    # total proton energy
    Ep_f = np.sqrt(P_f[is_ok]**2 + MP**2)


    # total energy of recoiling particle 
    Er_f = E_tot - Ep_f

    # invariant mass of recoiling particle
    M_r2 = Er_f**2 - Pm[is_ok]**2
    sel_Mr2 = M_r2 > 0.0
    M_r = np.sqrt(M_r2[sel_Mr2])

    h_Mr = B.histo(M_r, range = (0.5,1.5), bins = 200, 
                   xlabel = 'M_recoil', 
                   ylabel = 'Counts', 
                   title = 'Inv. mass of recoiling particle')


    h_Mr.plot(filled = False)

