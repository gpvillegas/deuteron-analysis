#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 11:13:22 2025

@author: gvill
"""

import data_init as D
import cut_handler as C
import numpy as np
import LT.box as B
from matplotlib import colormaps as cmp
import database_operations as db


#%%
TSP_branches = ['evNumber','P.BCM1.scalerCurrent','P.BCM4A.scalerCurrent',
                'P.BCM4C.scalerCurrent','P.BCM1.scalerChargeCut',
                'P.BCM4A.scalerChargeCut','P.BCM4C.scalerChargeCut',
                'P.BCM4A.scalerCharge',
                'P.pTRIG6.scaler']
TSH_branches = ['H.BCM4A.scalerCharge']

R = 20871

ROOTfiles_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

RUN = D.DATA_INIT(data_type='deut23_data',run=R,
                         select_branches={'TSP':TSP_branches,'TSH':TSH_branches},
                         select_trees=['TSP','TSH'],
                         ROOTfiles_path=ROOTfiles_DIR) 

#%% testing scaler cut
Acurr_cut = C.current_cut(RUN,current='BCM4A',cut_lim='tgt_boil')
Acurr_cut_arr = Acurr_cut()

charge = RUN.Branches['P.BCM4A.scalerCharge']
n = Acurr_cut.current.size
pos = np.arange(0,n,1)
inter = np.vstack((Acurr_cut.cut_array,np.roll(Acurr_cut.cut_array,-1))) 

j=0
temp = []
t = [-1,-1]
for i in inter.T:
    if i[0] == True and i[1] == False or i[0] == False and i[1] == True:
        print(i)
        print('Found interval at pos ', pos[j])
        if i[0] == True:
            t[0] = int(pos[j])
        elif i[0] == False:
            t[1] = int(pos[j]) +1
            
    if t[0]!= -1 and t[1] != -1:
        print(t)
        temp.append(t)
        t = [-1,-1]
    j+=1

bad = 0
for i in temp:
    try:
        bad += charge[i[1]] - charge[i[0]]
    except IndexError:
        bad += charge[i[1]-1] - charge[i[0]]     

#%%


M1curr_cut = C.current_cut(RUN,cut_lim='tgt_boil',current='BCM1')
M1curr_cut_arr = M1curr_cut()

Ccurr_cut = C.current_cut(RUN,cut_lim='tgt_boil',current='BCM4C')
Ccurr_cut_arr = Ccurr_cut()

Acharge = RUN.Branches['P.BCM4A.scalerChargeCut']
M1charge = RUN.Branches['P.BCM1.scalerChargeCut']
Ccharge = RUN.Branches['P.BCM4C.scalerChargeCut']

T6scl = RUN.Branches['P.pTRIG6.scaler']


print('T6 scl tot - BCM4A cut = ', T6scl[Acurr_cut.cut_array][-1])
print('T6 scl tot - BCM1 cut = ', T6scl[M1curr_cut.cut_array][-1])
print('T6 scl tot - BCM4C cut = ', T6scl[Ccurr_cut.cut_array][-1])

# bad_A_scl = T6scl[~Acurr_cut.cut_array]
# bad_A_tot = 0.

# for i in range(bad_A_scl.size-1):
#     bad_A_tot = bad_A_tot + (bad_A_scl[i+1] - bad_A_scl[i])
    
# print('bad_A_tot = ', bad_A_tot)

good_A_scl = T6scl[Acurr_cut.cut_array]
good_A_tot = 0.

for i in range(good_A_scl.size-1):
    good_A_tot = good_A_tot + (good_A_scl[i+1] - good_A_scl[i])
    
print('good_A_tot = ', good_A_tot)

good_M1_scl = T6scl[M1curr_cut.cut_array]
good_M1_tot = 0.

for i in range(good_M1_scl.size-1):
    good_M1_tot = good_M1_tot + (good_M1_scl[i+1] - good_M1_scl[i])
    
print('good_M1_tot = ', good_M1_tot)

good_C_scl = T6scl[Ccurr_cut.cut_array]
good_C_tot = 0.

for i in range(good_C_scl.size-1):
    good_C_tot = good_C_tot + (good_C_scl[i+1] - good_C_scl[i])
    
print('good_C_tot = ', good_C_tot)

#%% comparing bcm calib 
TSP_branches = ['evcount','evNumber','P.BCM1.scalerCurrent','P.BCM2.scalerCurrent',
                'P.BCM4A.scalerCurrent','P.BCM4B.scalerCurrent',
                'P.BCM4C.scalerCurrent','P.BCM1.scalerCharge',
                'P.BCM2.scalerCharge','P.BCM4A.scalerCharge',
                'P.BCM4B.scalerCharge','P.BCM4C.scalerCharge',
                'P.1MHz.scalerTime']

# R = [17756,18531,18723,19783]
R = 21059

BCMcalib_current_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

XEM2_BCMcalib_DIR = "/media/gvill/Gema's T7/ROOTfiles/xem2BCMcalib/"

RUN_old = D.DATA_INIT(data_type='deut23_data',run=R,
                  select_branches={'TSP':TSP_branches},
                  select_trees=['TSP'],
                  ROOTfiles_path=BCMcalib_current_DIR)

RUN_xem2 = D.DATA_INIT(data_type='deut23_data',run=R,
                  select_branches={'TSP':TSP_branches},
                  select_trees=['TSP'],
                  ROOTfiles_path=XEM2_BCMcalib_DIR)

#%%
# initialize current cuts
BCM1_cut_old = C.current_cut(RUN_old,current='BCM1',many=True)
BCM2_cut_old = C.current_cut(RUN_old,current='BCM2',many=True)
BCM4A_cut_old = C.current_cut(RUN_old,current='BCM4A',many=True)
BCM4B_cut_old = C.current_cut(RUN_old,current='BCM4B',many=True)
BCM4C_cut_old = C.current_cut(RUN_old,current='BCM4C',many=True)

BCM1_cut_old()
BCM2_cut_old()
BCM4A_cut_old()
BCM4B_cut_old()
BCM4C_cut_old()

charge_BCM1_cut_old = {}
charge_BCM2_cut_old = {}
charge_BCM4A_cut_old = {}
charge_BCM4B_cut_old = {}
charge_BCM4C_cut_old = {}
sclTime_BCM1_cut_old = {}
sclTime_BCM2_cut_old = {}
sclTime_BCM4A_cut_old = {}
sclTime_BCM4B_cut_old = {}
sclTime_BCM4C_cut_old = {}
charge_old = {}
avg_current_old = {}

for r in R:

    # get charge with current cut
    charge_BCM1_cut_old[r] = BCM1_cut_old.scaler(RUN_old.Branches[r]['P.BCM1.scalerCharge'],BCM1_cut_old.cut_array[r])
    charge_BCM2_cut_old[r] = BCM2_cut_old.scaler(RUN_old.Branches[r]['P.BCM2.scalerCharge'],BCM2_cut_old.cut_array[r])
    charge_BCM4A_cut_old[r] = BCM4A_cut_old.scaler(RUN_old.Branches[r]['P.BCM4A.scalerCharge'],BCM4A_cut_old.cut_array[r])
    charge_BCM4B_cut_old[r] = BCM4B_cut_old.scaler(RUN_old.Branches[r]['P.BCM4B.scalerCharge'],BCM4B_cut_old.cut_array[r])
    charge_BCM4C_cut_old[r] = BCM4C_cut_old.scaler(RUN_old.Branches[r]['P.BCM4C.scalerCharge'],BCM4C_cut_old.cut_array[r])
    
    # get scl time with current cut
    sclTime_BCM1_cut_old[r] = BCM1_cut_old.scaler(RUN_old.Branches[r]['P.1MHz.scalerTime'],BCM1_cut_old.cut_array[r])
    sclTime_BCM2_cut_old[r] = BCM2_cut_old.scaler(RUN_old.Branches[r]['P.1MHz.scalerTime'],BCM2_cut_old.cut_array[r])
    sclTime_BCM4A_cut_old[r] = BCM4A_cut_old.scaler(RUN_old.Branches[r]['P.1MHz.scalerTime'],BCM4A_cut_old.cut_array[r])
    sclTime_BCM4B_cut_old[r] = BCM4B_cut_old.scaler(RUN_old.Branches[r]['P.1MHz.scalerTime'],BCM4B_cut_old.cut_array[r])
    sclTime_BCM4C_cut_old[r] = BCM4C_cut_old.scaler(RUN_old.Branches[r]['P.1MHz.scalerTime'],BCM4C_cut_old.cut_array[r])
    
    charge_old[r] = [charge_BCM1_cut_old[r],charge_BCM2_cut_old[r],charge_BCM4A_cut_old[r],
                  charge_BCM4B_cut_old[r],charge_BCM4C_cut_old[r]]
    
    avg_current_old[r] = [charge_BCM1_cut_old[r]/sclTime_BCM1_cut_old[r],
                       charge_BCM2_cut_old[r]/sclTime_BCM2_cut_old[r],
                       charge_BCM4A_cut_old[r]/sclTime_BCM4A_cut_old[r],
                       charge_BCM4B_cut_old[r]/sclTime_BCM4B_cut_old[r],
                       charge_BCM4C_cut_old[r]/sclTime_BCM4C_cut_old[r]]

# initialize current cuts
BCM1_cut_xem2 = C.current_cut(RUN_xem2,current='BCM1',many=True)
BCM2_cut_xem2 = C.current_cut(RUN_xem2,current='BCM2',many=True)
BCM4A_cut_xem2 = C.current_cut(RUN_xem2,current='BCM4A',many=True)
BCM4B_cut_xem2 = C.current_cut(RUN_xem2,current='BCM4B',many=True)
BCM4C_cut_xem2 = C.current_cut(RUN_xem2,current='BCM4C',many=True)

BCM1_cut_xem2()
BCM2_cut_xem2()
BCM4A_cut_xem2()
BCM4B_cut_xem2()
BCM4C_cut_xem2()

charge_BCM1_cut_xem2 = {}
charge_BCM2_cut_xem2 = {}
charge_BCM4A_cut_xem2 = {}
charge_BCM4B_cut_xem2 = {}
charge_BCM4C_cut_xem2 = {}
sclTime_BCM1_cut_xem2 = {}
sclTime_BCM2_cut_xem2 = {}
sclTime_BCM4A_cut_xem2 = {}
sclTime_BCM4B_cut_xem2 = {}
sclTime_BCM4C_cut_xem2 = {}
charge_xem2 = {}
avg_current_xem2 = {}

for r in R:

    # get charge with current cut
    charge_BCM1_cut_xem2[r] = BCM1_cut_xem2.scaler(RUN_xem2.Branches[r]['P.BCM1.scalerCharge'],BCM1_cut_xem2.cut_array[r])
    charge_BCM2_cut_xem2[r] = BCM2_cut_xem2.scaler(RUN_xem2.Branches[r]['P.BCM2.scalerCharge'],BCM2_cut_xem2.cut_array[r])
    charge_BCM4A_cut_xem2[r] = BCM4A_cut_xem2.scaler(RUN_xem2.Branches[r]['P.BCM4A.scalerCharge'],BCM4A_cut_xem2.cut_array[r])
    charge_BCM4B_cut_xem2[r] = BCM4B_cut_xem2.scaler(RUN_xem2.Branches[r]['P.BCM4B.scalerCharge'],BCM4B_cut_xem2.cut_array[r])
    charge_BCM4C_cut_xem2[r] = BCM4C_cut_xem2.scaler(RUN_xem2.Branches[r]['P.BCM4C.scalerCharge'],BCM4C_cut_xem2.cut_array[r])
    
    # get scl time with current cut
    sclTime_BCM1_cut_xem2[r] = BCM1_cut_xem2.scaler(RUN_xem2.Branches[r]['P.1MHz.scalerTime'],BCM1_cut_xem2.cut_array[r])
    sclTime_BCM2_cut_xem2[r] = BCM2_cut_xem2.scaler(RUN_xem2.Branches[r]['P.1MHz.scalerTime'],BCM2_cut_xem2.cut_array[r])
    sclTime_BCM4A_cut_xem2[r] = BCM4A_cut_xem2.scaler(RUN_xem2.Branches[r]['P.1MHz.scalerTime'],BCM4A_cut_xem2.cut_array[r])
    sclTime_BCM4B_cut_xem2[r] = BCM4B_cut_xem2.scaler(RUN_xem2.Branches[r]['P.1MHz.scalerTime'],BCM4B_cut_xem2.cut_array[r])
    sclTime_BCM4C_cut_xem2[r] = BCM4C_cut_xem2.scaler(RUN_xem2.Branches[r]['P.1MHz.scalerTime'],BCM4C_cut_xem2.cut_array[r])
    
    charge_xem2[r] = [charge_BCM1_cut_xem2[r],charge_BCM2_cut_xem2[r],charge_BCM4A_cut_xem2[r],
                  charge_BCM4B_cut_xem2[r],charge_BCM4C_cut_xem2[r]]
    
    avg_current_xem2[r] = [charge_BCM1_cut_xem2[r]/sclTime_BCM1_cut_xem2[r],
                       charge_BCM2_cut_xem2[r]/sclTime_BCM2_cut_xem2[r],
                       charge_BCM4A_cut_xem2[r]/sclTime_BCM4A_cut_xem2[r],
                       charge_BCM4B_cut_xem2[r]/sclTime_BCM4B_cut_xem2[r],
                       charge_BCM4C_cut_xem2[r]/sclTime_BCM4C_cut_xem2[r]]

#%%

# initialize current cuts
BCM1_cut_old = C.current_cut(RUN_old,current='BCM1')
BCM2_cut_old = C.current_cut(RUN_old,current='BCM2')
BCM4A_cut_old = C.current_cut(RUN_old,current='BCM4A')
BCM4B_cut_old = C.current_cut(RUN_old,current='BCM4B')
BCM4C_cut_old = C.current_cut(RUN_old,current='BCM4C')

BCM1_cut_old()
BCM2_cut_old()
BCM4A_cut_old()
BCM4B_cut_old()
BCM4C_cut_old()

# get charge with current cut
charge_BCM1_cut_old = BCM1_cut_old.scaler(RUN_old.Branches['P.BCM1.scalerCharge'],BCM1_cut_old.cut_array)
charge_BCM2_cut_old = BCM2_cut_old.scaler(RUN_old.Branches['P.BCM2.scalerCharge'],BCM2_cut_old.cut_array)
charge_BCM4A_cut_old = BCM4A_cut_old.scaler(RUN_old.Branches['P.BCM4A.scalerCharge'],BCM4A_cut_old.cut_array)
charge_BCM4B_cut_old = BCM4B_cut_old.scaler(RUN_old.Branches['P.BCM4B.scalerCharge'],BCM4B_cut_old.cut_array)
charge_BCM4C_cut_old = BCM4C_cut_old.scaler(RUN_old.Branches['P.BCM4C.scalerCharge'],BCM4C_cut_old.cut_array)

# get scl time with current cut
sclTime_BCM1_cut_old = BCM1_cut_old.scaler(RUN_old.Branches['P.1MHz.scalerTime'],BCM1_cut_old.cut_array)
sclTime_BCM2_cut_old = BCM2_cut_old.scaler(RUN_old.Branches['P.1MHz.scalerTime'],BCM2_cut_old.cut_array)
sclTime_BCM4A_cut_old = BCM4A_cut_old.scaler(RUN_old.Branches['P.1MHz.scalerTime'],BCM4A_cut_old.cut_array)
sclTime_BCM4B_cut_old = BCM4B_cut_old.scaler(RUN_old.Branches['P.1MHz.scalerTime'],BCM4B_cut_old.cut_array)
sclTime_BCM4C_cut_old = BCM4C_cut_old.scaler(RUN_old.Branches['P.1MHz.scalerTime'],BCM4C_cut_old.cut_array)

charge_old = [charge_BCM1_cut_old,charge_BCM2_cut_old,charge_BCM4A_cut_old,
              charge_BCM4B_cut_old,charge_BCM4C_cut_old]

avg_current_old = [charge_BCM1_cut_old/sclTime_BCM1_cut_old,
                   charge_BCM2_cut_old/sclTime_BCM2_cut_old,
                   charge_BCM4A_cut_old/sclTime_BCM4A_cut_old,
                   charge_BCM4B_cut_old/sclTime_BCM4B_cut_old,
                   charge_BCM4C_cut_old/sclTime_BCM4C_cut_old]

# initialize current cuts
BCM1_cut_xem2 = C.current_cut(RUN_xem2,current='BCM1')
BCM2_cut_xem2 = C.current_cut(RUN_xem2,current='BCM2')
BCM4A_cut_xem2 = C.current_cut(RUN_xem2,current='BCM4A')
BCM4B_cut_xem2 = C.current_cut(RUN_xem2,current='BCM4B')
BCM4C_cut_xem2 = C.current_cut(RUN_xem2,current='BCM4C')

BCM1_cut_xem2()
BCM2_cut_xem2()
BCM4A_cut_xem2()
BCM4B_cut_xem2()
BCM4C_cut_xem2()

# get charge with current cut
charge_BCM1_cut_xem2 = BCM1_cut_xem2.scaler(RUN_xem2.Branches['P.BCM1.scalerCharge'],BCM1_cut_xem2.cut_array)
charge_BCM2_cut_xem2 = BCM2_cut_xem2.scaler(RUN_xem2.Branches['P.BCM2.scalerCharge'],BCM2_cut_xem2.cut_array)
charge_BCM4A_cut_xem2 = BCM4A_cut_xem2.scaler(RUN_xem2.Branches['P.BCM4A.scalerCharge'],BCM4A_cut_xem2.cut_array)
charge_BCM4B_cut_xem2 = BCM4B_cut_xem2.scaler(RUN_xem2.Branches['P.BCM4B.scalerCharge'],BCM4B_cut_xem2.cut_array)
charge_BCM4C_cut_xem2 = BCM4C_cut_xem2.scaler(RUN_xem2.Branches['P.BCM4C.scalerCharge'],BCM4C_cut_xem2.cut_array)

# get scl time with current cut
sclTime_BCM1_cut_xem2 = BCM1_cut_xem2.scaler(RUN_xem2.Branches['P.1MHz.scalerTime'],BCM1_cut_xem2.cut_array)
sclTime_BCM2_cut_xem2 = BCM2_cut_xem2.scaler(RUN_xem2.Branches['P.1MHz.scalerTime'],BCM2_cut_xem2.cut_array)
sclTime_BCM4A_cut_xem2 = BCM4A_cut_xem2.scaler(RUN_xem2.Branches['P.1MHz.scalerTime'],BCM4A_cut_xem2.cut_array)
sclTime_BCM4B_cut_xem2 = BCM4B_cut_xem2.scaler(RUN_xem2.Branches['P.1MHz.scalerTime'],BCM4B_cut_xem2.cut_array)
sclTime_BCM4C_cut_xem2 = BCM4C_cut_xem2.scaler(RUN_xem2.Branches['P.1MHz.scalerTime'],BCM4C_cut_xem2.cut_array)

charge_xem2 = [charge_BCM1_cut_xem2,charge_BCM2_cut_xem2,charge_BCM4A_cut_xem2,
              charge_BCM4B_cut_xem2,charge_BCM4C_cut_xem2]

avg_current_xem2 = [charge_BCM1_cut_xem2/sclTime_BCM1_cut_xem2,
                   charge_BCM2_cut_xem2/sclTime_BCM2_cut_xem2,
                   charge_BCM4A_cut_xem2/sclTime_BCM4A_cut_xem2,
                   charge_BCM4B_cut_xem2/sclTime_BCM4B_cut_xem2,
                   charge_BCM4C_cut_xem2/sclTime_BCM4C_cut_xem2]

#%% plot currents: define histograms

current_BCM1_old = RUN_old.Branches['P.BCM1.scalerCurrent']
current_BCM2_old = RUN_old.Branches['P.BCM2.scalerCurrent']
current_BCM4A_old = RUN_old.Branches['P.BCM4A.scalerCurrent']
current_BCM4B_old = RUN_old.Branches['P.BCM4B.scalerCurrent']
current_BCM4C_old = RUN_old.Branches['P.BCM4C.scalerCurrent']

current_BCM1_xem2 = RUN_xem2.Branches['P.BCM1.scalerCurrent']
current_BCM2_xem2 = RUN_xem2.Branches['P.BCM2.scalerCurrent']
current_BCM4A_xem2 = RUN_xem2.Branches['P.BCM4A.scalerCurrent']
current_BCM4B_xem2 = RUN_xem2.Branches['P.BCM4B.scalerCurrent']
current_BCM4C_xem2 = RUN_xem2.Branches['P.BCM4C.scalerCurrent']

h_curr_BCM1_old = B.histo(current_BCM1_old,range=(-20,70),bins=100)
h_curr_BCM2_old = B.histo(current_BCM2_old,range=(-20,70),bins=100)
h_curr_BCM4A_old = B.histo(current_BCM4A_old,range=(-20,70),bins=100)
h_curr_BCM4B_old = B.histo(current_BCM4B_old,range=(-20,70),bins=100)
h_curr_BCM4C_old = B.histo(current_BCM4C_old,range=(-20,70),bins=100)

h_curr_BCM1_xem2 = B.histo(current_BCM1_xem2,range=(-20,70),bins=100)
h_curr_BCM2_xem2 = B.histo(current_BCM2_xem2,range=(-20,70),bins=100)
h_curr_BCM4A_xem2 = B.histo(current_BCM4A_xem2,range=(-20,70),bins=100)
h_curr_BCM4B_xem2 = B.histo(current_BCM4B_xem2,range=(-20,70),bins=100)
h_curr_BCM4C_xem2 = B.histo(current_BCM4C_xem2,range=(-20,70),bins=100)

current_histos = {'BCM1':[h_curr_BCM1_old,h_curr_BCM1_xem2],
                  'BCM2':[h_curr_BCM2_old,h_curr_BCM2_xem2],
                  'BCM4A':[h_curr_BCM4A_old,h_curr_BCM4A_xem2],
                  'BCM4B':[h_curr_BCM4B_old,h_curr_BCM4B_xem2],
                  'BCM4C':[h_curr_BCM4C_old,h_curr_BCM4C_xem2]}

#%% plot currents
B.pl.figure(figsize=(20,5),layout='constrained')
B.pl.suptitle(R)
i=1

for h in current_histos:
    B.pl.subplot(1,5,i)
    
    current_histos[h][0].plot(color=cmp['Paired'].colors[5],
                              label='old')
    current_histos[h][1].plot(filled=False,color='black',
                              label='xem2')
    B.pl.title(h)
    B.pl.xlabel('')
    B.pl.ylabel('')
    B.pl.legend()
    i+=1
    
#%% plot charge vs currents
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14
}
B.pl.rc('font', **font)
B.pl.figure()

current_labels = ['BCM1','BCM2','BCM4A','BCM4B','BCM4C']

B.plot_exp(avg_current_old, np.array(charge_old)/1e3, label='old',
           markersize=8)
B.plot_exp(avg_current_xem2, np.array(charge_xem2)/1e3, label='xem2',
           markersize=8)

for i, label in enumerate(current_labels):
    B.pl.annotate(label, (avg_current_old[i], np.array(charge_old)[i]/1e3),
                  textcoords="offset points", xytext=(5,-10), ha='left')
    
for i, label in enumerate(current_labels):
    B.pl.annotate(label, (avg_current_xem2[i], np.array(charge_xem2)[i]/1e3),
                  textcoords="offset points", xytext=(5,10), ha='right')

B.pl.title(R)
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.ylabel('Charge [mC]')
B.pl.legend()

#%% write results to a file
with open('BCMcalib_comparison_results.txt','a') as f:
    header = '#! run[i,0]/calib[s,1]/BCM1_charge[f,2]/BCM1_avg_current[f,3]/'+\
        'BCM2_charge[f,4]/BCM2_avg_current[f,5]/BCM4A_charge[f,6]/BCM4A_avg_current[f,7]/'+\
        'BCM4B_charge[f,8]/BCM4B_avg_current[f,9]/BCM4C_charge[f,10]/BCM4C_avg_current[f,11]/'+\
        'BCM1_ratio[f,12]/BCM2_ratio[f,13]/BCM4A_ratio[f,14]/BCM4B_ratio[f,15]/'+\
        'BCM4C_ratio[f,16]/\n'
    res_old = f'{R} old {charge_BCM1_cut_old/1e3:0.4f} {charge_BCM1_cut_old/sclTime_BCM1_cut_old:0.4f}'+\
        f' {charge_BCM2_cut_old/1e3:0.4f} {charge_BCM2_cut_old/sclTime_BCM2_cut_old:0.4f}'+\
        f' {charge_BCM4A_cut_old/1e3:0.4f} {charge_BCM4A_cut_old/sclTime_BCM4A_cut_old:0.4f}'+\
        f' {charge_BCM4B_cut_old/1e3:0.4f} {charge_BCM4B_cut_old/sclTime_BCM4B_cut_old:0.4f}'+\
        f' {charge_BCM4C_cut_old/1e3:0.4f} {charge_BCM4C_cut_old/sclTime_BCM4C_cut_old:0.4f}'+\
        f' {charge_BCM1_cut_old/charge_BCM1_cut_xem2:0.4f}'+\
        f' {charge_BCM2_cut_old/charge_BCM2_cut_xem2:0.4f}'+\
        f' {charge_BCM4A_cut_old/charge_BCM4A_cut_xem2:0.4f}'+\
        f' {charge_BCM4B_cut_old/charge_BCM4B_cut_xem2:0.4f}'+\
        f' {charge_BCM4C_cut_old/charge_BCM4C_cut_xem2:0.4f}\n'
    res_xem2 = f'{R} xem2 {charge_BCM1_cut_xem2/1e3:0.4f} {charge_BCM1_cut_xem2/sclTime_BCM1_cut_xem2:0.4f}'+\
        f' {charge_BCM2_cut_xem2/1e3:0.4f} {charge_BCM2_cut_xem2/sclTime_BCM2_cut_xem2:0.4f}'+\
        f' {charge_BCM4A_cut_xem2/1e3:0.4f} {charge_BCM4A_cut_xem2/sclTime_BCM4A_cut_xem2:0.4f}'+\
        f' {charge_BCM4B_cut_xem2/1e3:0.4f} {charge_BCM4B_cut_xem2/sclTime_BCM4B_cut_xem2:0.4f}'+\
        f' {charge_BCM4C_cut_xem2/1e3:0.4f} {charge_BCM4C_cut_xem2/sclTime_BCM4C_cut_xem2:0.4f}'+\
        f' {charge_BCM1_cut_old/charge_BCM1_cut_xem2:0.4f}'+\
        f' {charge_BCM2_cut_old/charge_BCM2_cut_xem2:0.4f}'+\
        f' {charge_BCM4A_cut_old/charge_BCM4A_cut_xem2:0.4f}'+\
        f' {charge_BCM4B_cut_old/charge_BCM4B_cut_xem2:0.4f}'+\
        f' {charge_BCM4C_cut_old/charge_BCM4C_cut_xem2:0.4f}\n'
        
    # f.write(header)
    f.write(res_old)
    f.write(res_xem2)
    
#%% plot ratios

fname = 'BCMcalib_comparison_results.txt'
f = B.get_file(fname)

runs = f['run']
BCM1ratio = f['BCM1_ratio']
BCM2ratio = f['BCM2_ratio']
BCM4Aratio = f['BCM4A_ratio']
BCM4Bratio = f['BCM4B_ratio']
BCM4Cratio = f['BCM4C_ratio']

def rm_odd(arr):
    new_arr = []
    for i in range(arr.size):
        if i == 0 or i%2 == 0:
            new_arr.append(arr[i])
    return np.array(new_arr)

runs = rm_odd(runs)
BCM1ratio = rm_odd(BCM1ratio)
BCM2ratio = rm_odd(BCM2ratio)
BCM4Aratio = rm_odd(BCM4Aratio)
BCM4Bratio = rm_odd(BCM4Bratio)
BCM4Cratio = rm_odd(BCM4Cratio)

B.pl.figure()
B.plot_line(runs,BCM1ratio,color=cmp['Paired'].colors[0])
B.plot_line(runs,BCM2ratio,color=cmp['Paired'].colors[2])
B.plot_line(runs,BCM4Aratio,color=cmp['Paired'].colors[4])
B.plot_line(runs,BCM4Bratio,color=cmp['Paired'].colors[6])
B.plot_line(runs,BCM4Cratio,color=cmp['Paired'].colors[8])

B.plot_exp(runs,BCM1ratio,label='BCM1',markersize=8,color=cmp['Paired'].colors[1])
B.plot_exp(runs,BCM2ratio,label='BCM2',markersize=8,color=cmp['Paired'].colors[3])
B.plot_exp(runs,BCM4Aratio,label='BCM4A',markersize=8,color=cmp['Paired'].colors[5])
B.plot_exp(runs,BCM4Bratio,label='BCM4B',markersize=8,color=cmp['Paired'].colors[7])
B.plot_exp(runs,BCM4Cratio,label='BCM4C',markersize=8,color=cmp['Paired'].colors[9])

ax= B.pl.gca()
ax.set_xticks(runs)

B.pl.xlabel('Run #')
B.pl.ylabel('Charge Ratio old/xem2')

B.pl.legend()

#%% plot current vs evcount
BCM1_current_old = {}
BCM2_current_old = {}
BCM4A_current_old = {}
BCM4C_current_old = {}
BCM1_current_xem2 = {}
BCM2_current_xem2 = {}
BCM4A_current_xem2 = {}
BCM4C_current_xem2 = {}
evCount = {}

for r in R:
    BCM1_current_old[r] = RUN_old.Branches[r]['P.BCM1.scalerCurrent']
    BCM2_current_old[r] = RUN_old.Branches[r]['P.BCM2.scalerCurrent']
    BCM4A_current_old[r] = RUN_old.Branches[r]['P.BCM4A.scalerCurrent']
    BCM4C_current_old[r] = RUN_old.Branches[r]['P.BCM4C.scalerCurrent']
    BCM1_current_xem2[r] = RUN_xem2.Branches[r]['P.BCM1.scalerCurrent']
    BCM2_current_xem2[r] = RUN_xem2.Branches[r]['P.BCM2.scalerCurrent']
    BCM4A_current_xem2[r] = RUN_xem2.Branches[r]['P.BCM4A.scalerCurrent']
    BCM4C_current_xem2[r] = RUN_xem2.Branches[r]['P.BCM4C.scalerCurrent']
    evCount[r] = RUN_old.Branches[r]['evcount']
    


for r in R:
    B.pl.figure(figsize=(10,20),layout='constrained')
    i=1
    B.pl.suptitle(f'Run {r}')
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM1_current_old[r],color=cmp['Paired'].colors[1],
               label='old')
    B.pl.title('BCM1')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM1_current_xem2[r],color=cmp['Paired'].colors[5],
           label='xem2')
    B.pl.title('BCM1')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM2_current_old[r],color=cmp['Paired'].colors[1],
               label='old')
    B.pl.title('BCM2')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM2_current_xem2[r],color=cmp['Paired'].colors[5],
               label='xem2')
    B.pl.title('BCM2')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM4A_current_old[r],color=cmp['Paired'].colors[1],
               label='old')
    B.pl.title('BCM4A')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM4A_current_xem2[r],color=cmp['Paired'].colors[5],
               label='xem2')
    B.pl.title('BCM4A')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM4C_current_old[r],color=cmp['Paired'].colors[1],
               label='old')
    B.pl.title('BCM4C')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()
    i+=1
    B.pl.subplot(4,2,i)
    B.plot_exp(evCount[r],BCM4C_current_xem2[r],color=cmp['Paired'].colors[5],
               label='xem2')
    B.pl.title('BCM4C')
    B.pl.xlabel('Scaler Event Count')
    B.pl.ylabel('Scaler Current')
    B.pl.legend()

#%%plot current ratios
for r in R:
    B.pl.figure(figsize=(5,20),layout='constrained')
    i=1
    B.pl.suptitle(f'Run {r}')
    B.pl.subplot(4,1,i)
    B.plot_exp(evCount[r],BCM1_current_old[r]/BCM1_current_xem2[r],color=cmp['Paired'].colors[1])
    B.pl.title('BCM1')
    B.pl.ylim(0.95,1.15)
    B.pl.hlines(1.0,-80,2200,color='black')


    i+=1
    B.pl.subplot(4,1,i)
    B.plot_exp(evCount[r],BCM2_current_old[r]/BCM2_current_xem2[r],color=cmp['Paired'].colors[1])
    B.pl.title('BCM2')
    B.pl.ylim(0.95,1.15) 
    B.pl.hlines(1.0,-80,2200,color='black')


    i+=1
    B.pl.subplot(4,1,i)
    B.plot_exp(evCount[r],BCM4A_current_old[r]/BCM4A_current_xem2[r],color=cmp['Paired'].colors[1])
    B.pl.title('BCM4A')
    B.pl.ylim(-1,2)
    B.pl.ylabel('Scaler Current Ratio old/xem2')

    i+=1
    B.pl.subplot(4,1,i)
    B.plot_exp(evCount[r],BCM4C_current_old[r]/BCM4C_current_xem2[r],color=cmp['Paired'].colors[1])
    B.pl.title('BCM4C')
    B.pl.ylim(-1,5)
    B.pl.xlabel('Scaler Event Count')


#%% plot some charge ratios for deuteron runs 
# To get rid of unhelpful groupings (lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0]

db_info = db.retrieve('deuteron_db.db','run, BCM1_charge, BCM2_charge, BCM4A_charge, BCM1_current',
                      'RUN_LIST_UPDATED')

RUN_LIST = []
BCM1_charge = []
BCM2_charge = []
BCM4A_charge = []
BCM1_current = []

for d in db_info:
    run, bcm1, bcm2, bcm4a, bcm1c = d
    RUN_LIST.append(run)
    BCM1_charge.append(bcm1)
    BCM2_charge.append(bcm2)
    BCM4A_charge.append(bcm4a)
    BCM1_current.append(bcm1c)

BCM1_charge = np.array(BCM1_charge)
BCM2_charge = np.array(BCM2_charge)
BCM4A_charge = np.array(BCM4A_charge)

#%%    
B.pl.figure()
B.plot_exp(RUN_LIST, BCM1_charge/BCM4A_charge)
B.pl.title('old BCM Calib')
B.pl.xlabel('Run #')
B.pl.ylabel('BCM1_charge/BCM4A_charge')

B.pl.figure()
B.plot_exp(RUN_LIST, BCM2_charge/BCM4A_charge)
B.pl.title('old BCM Calib')
B.pl.xlabel('Run #')
B.pl.ylabel('BCM2_charge/BCM4A_charge')

B.pl.figure()
B.plot_exp(RUN_LIST, BCM1_charge/BCM2_charge)
B.pl.title('old BCM Calib')
B.pl.xlabel('Run #')
B.pl.ylabel('BCM1_charge/BCM2_charge')

B.pl.figure()
B.plot_exp(BCM1_current, BCM1_charge/BCM2_charge)
B.pl.title('old BCM Calib')
B.pl.xlabel('BCM1_avg_current')
B.pl.ylabel('BCM1_charge/BCM2_charge')
