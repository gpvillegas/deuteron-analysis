#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 13:04:40 2025

@author: gvill

Plot SHMS/HMS delta acceptance for all settings
"""
import numpy as np
import LT.box as B
import data_init as D
import cut_handler as C

#%%
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
    #tot_q, enorm_h = normalize_histos(histos,many)
    tot_q, enorm_h = 1,histos
    loop1 = True 
    for m in many:
        if loop1:
            hsum = enorm_h[m]
            loop1 = False
        else:
            hsum += enorm_h[m]
     
    qnorm_hsum = hsum*(1/tot_q)
    
    return qnorm_hsum

#%% directories for each type of run
# hc = heep_coin; d = deep

hc_root_DIR = "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim3_allfp/"
d_root_DIR = "/media/gvill/Gema's T7/ROOTfiles/deep_testing0/"

hc_SIMC_DIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim3/"
d_SIMC_DIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/deep_testing0/"

#%% select branches

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

br_sel = acc_var + ['H.kin.secondary.emiss_nuc','H.kin.secondary.emiss']

br_sel_SIMC = ['e_delta','h_delta','Em']

#%% load data and SIMC

hc = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin', 
                    select_branches={'T':br_sel},ROOTfiles_path=hc_root_DIR)

# hc_SIMC = D.DATA_INIT(data_type='SIMC', kin_study='heep_coin',
#                     select_branches={'SNT':br_sel_SIMC},simc_type='-',
#                     SIMC_ROOTfiles_path=d_SIMC_DIR)

#%%
R_pm120 = [20871,20872]
R_pm580 = [20873,20874,20875,20876,20877,20878,20880,20881,20882,20883,21076,
           21078,21079,21080,21081,21082,21083,21084,21085,21086,21087,21088,
           21089,21090,21091,21092,21093,21094,21095,21096,21097,21098,21099,
           21100,21101,21102]
R_pm800 = [20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,20896,
          20897,20898,20899,20900,20901,20902,20903,20905,20904,20907,20908,
          20909,20910,20911,20912,20913,20914,20915,20916,20921,20922,20923,
          20924,20925,20926,20927,20928,20929,20930,20931,20932,20933,20934,
          20935,20936,20937,20938,20939,20940,20941,20942,20943,20944,20945,
          20949,20950,20951,20953,20954,20955,20956]
R_pm900 = [20958,20959,20960,20961,20962,20963,20965,20966,20969,20970,20971,
           20972,20973,20974,20975,20976,20977,20978,20979,20980,20981,20982,
           20983,20984,20985,20986,20987,20988,20989,20990,20991,20992,20993,
           20994,20995,20996,20997,20998,20999,21000,21001,21002,21003,21004,
           21005,21006,21007,21008,21009,21011,21012,21013,21014,21015,21016,
           21017,21018,21019,21020,21021,21022,21023,21024,21025,21026,21027,
           21028,21029,21030,21031,21032,21033,21034,21036,21037,21038,
           21039,21040,21041,21042,21043,21044,21045,21046,21047,21048]

#%%

d_pm120 = D.DATA_INIT(data_type='deut23_data', run=R_pm120, 
                    select_branches={'T':br_sel},ROOTfiles_path=d_root_DIR)
#%%
d_pm580 = D.DATA_INIT(data_type='deut23_data', run=R_pm580, 
                    select_branches={'T':br_sel},ROOTfiles_path=d_root_DIR)
#%%
d_pm800 = D.DATA_INIT(data_type='deut23_data', run=R_pm800, 
                    select_branches={'T':br_sel},ROOTfiles_path=d_root_DIR)
#%%
d_pm900 = D.DATA_INIT(data_type='deut23_data', run=R_pm900, 
                    select_branches={'T':br_sel},ROOTfiles_path=d_root_DIR)

#%% Define cuts for data

# heep_coin cuts
hc_cuts = {}
for r in hc.many:
    hc_cuts[r] = C.Em_cut_hc(hc.Branches[r]['H.kin.secondary.emiss'])

# pm120 cuts
pm120_cuts = {}
for r in d_pm120.many:
    C.Em_cut.init()
    pm120_cuts[r] = C.Em_cut(d_pm120.Branches[r]['H.kin.secondary.emiss_nuc'])

# pm580 cuts
pm580_cuts = {}
for r in d_pm580.many:
    C.Em_cut.init()
    pm580_cuts[r] = C.Em_cut(d_pm580.Branches[r]['H.kin.secondary.emiss_nuc'])

# pm800 cuts
pm800_cuts = {}
for r in d_pm800.many:
    C.Em_cut.init()
    pm800_cuts[r] = C.Em_cut(d_pm800.Branches[r]['H.kin.secondary.emiss_nuc'])

# pm900 cuts
pm900_cuts = {}
for r in d_pm900.many:
    C.Em_cut.init()
    pm900_cuts[r] = C.Em_cut(d_pm900.Branches[r]['H.kin.secondary.emiss_nuc'])
    
    # add collimator cut
    # hxc = RUN.Branches[r]['H.extcor.xsieve']
    # hyc = RUN.Branches[r]['H.extcor.ysieve']
    
    # hcoll_cut = C.coll_cut(hxc, hyc, spec='HMS')
    # c_list.append(hcoll_cut) 
    
    # cuts_to_apply[r] = c_list

#%% apply cuts
pd_hc = {}
hd_hc = {}
for r in hc.many:
    pd_hc[r] = hc.Branches[r]['P.gtr.dp'][hc_cuts[r]]
    hd_hc[r] = hc.Branches[r]['H.gtr.dp'][hc_cuts[r]]
 
pd_pm120 = {}
hd_pm120 = {}  
for r in d_pm120.many:
    pd_pm120[r] = d_pm120.Branches[r]['P.gtr.dp'][pm120_cuts[r]]
    hd_pm120[r] = d_pm120.Branches[r]['H.gtr.dp'][pm120_cuts[r]]
    
pd_pm580 ={}
hd_pm580 = {}    
for r in d_pm580.many:
    pd_pm580[r] = d_pm580.Branches[r]['P.gtr.dp'][pm580_cuts[r]] 
    hd_pm580[r] = d_pm580.Branches[r]['H.gtr.dp'][pm580_cuts[r]] 
    
pd_pm800 ={} 
hd_pm800 = {}   
for r in d_pm800.many:
    pd_pm800[r] = d_pm800.Branches[r]['P.gtr.dp'][pm800_cuts[r]]  
    hd_pm800[r] = d_pm800.Branches[r]['H.gtr.dp'][pm800_cuts[r]]
    
pd_pm900 ={} 
hd_pm900 = {}   
for r in d_pm900.many:
    pd_pm900[r] = d_pm900.Branches[r]['P.gtr.dp'][pm900_cuts[r]] 
    hd_pm900[r] = d_pm900.Branches[r]['H.gtr.dp'][pm900_cuts[r]]

#%%
pdelta_hc_histo = []
hdelta_hc_histo = []
for r in hc.many:
    if r != 20869:
        h = B.histo(pd_hc[r], range=(-15,15),bins=50)
        pdelta_hc_histo.append(h)
    
        h = B.histo(hd_hc[r], range=(-15,15),bins=50)
        hdelta_hc_histo.append(h)
        
pdelta_pm120_histo = {}
hdelta_pm120_histo = {}
for r in d_pm120.many:
    h = B.histo(pd_pm120[r], range=(-15,15),bins=50)
    pdelta_pm120_histo[r] = h 
    
    h = B.histo(hd_pm120[r], range=(-15,15),bins=50)
    hdelta_pm120_histo[r] = h  

pdelta_pm580_histo = {}
hdelta_pm580_histo = {}
for r in d_pm580.many:
    h = B.histo(pd_pm580[r], range=(-15,15),bins=50)
    pdelta_pm580_histo[r] = h
    
    h = B.histo(hd_pm580[r], range=(-15,15),bins=50)
    hdelta_pm580_histo[r] = h

pdelta_pm800_histo = {}
hdelta_pm800_histo = {}
for r in d_pm800.many:
    h = B.histo(pd_pm800[r], range=(-15,15),bins=50)
    pdelta_pm800_histo[r] = h 
    
    h = B.histo(hd_pm800[r], range=(-15,15),bins=50)
    hdelta_pm800_histo[r] = h 

pdelta_pm900_histo = {}
hdelta_pm900_histo = {}
for r in d_pm900.many:
    h = B.histo(pd_pm900[r], range=(-15,15),bins=50)
    pdelta_pm900_histo[r] = h 

    h = B.histo(hd_pm900[r], range=(-15,15),bins=50)
    hdelta_pm900_histo[r] = h           
         

#%%
pdelta_pm120_histo_all = combine_histos(pdelta_pm120_histo,d_pm120.many)
pdelta_pm580_histo_all = combine_histos(pdelta_pm580_histo,d_pm580.many)
pdelta_pm800_histo_all = combine_histos(pdelta_pm800_histo,d_pm800.many)
pdelta_pm900_histo_all = combine_histos(pdelta_pm900_histo,d_pm900.many)

hdelta_pm120_histo_all = combine_histos(hdelta_pm120_histo,d_pm120.many)
hdelta_pm580_histo_all = combine_histos(hdelta_pm580_histo,d_pm580.many)
hdelta_pm800_histo_all = combine_histos(hdelta_pm800_histo,d_pm800.many)
hdelta_pm900_histo_all = combine_histos(hdelta_pm900_histo,d_pm900.many)

#%%SHMS delta
colors = ['#b7e6a5','#7ccba2','#46aea0','#089099','#00718b','#045275','#003147']
linecolors = ['#165682','#2c6d97','#4085ac','#559dc1','#6cb6d6','#83cfea','#9ce9ff']
B.pl.figure()


# pdelta_pm800_histo_all.plot(color='#AF58BA',alpha=0.4,label='pm_800')
# pdelta_pm580_histo_all.plot(color='#009ADE',alpha=0.3,label='pm_580') 
pdelta_pm900_histo_all.plot(color='#FFC61E',alpha=0.5,label='pm_900')
# pdelta_pm120_histo_all.plot(color='#FF1F5B',alpha=0.5,label='pm_120')

# for h in range(len(pdelta_hc_histo)):
#     pdelta_hc_histo[h].plot(color=colors[h],alpha=0.3)

B.pl.vlines(-7,0,2000,linestyles='--',colors='black')
ax = B.pl.gca()
ax.set_autoscaley_on(True)

B.pl.title('SHMS $\delta$ w/ Em cut')
B.pl.xlabel('')
B.pl.ylabel('')
B.pl.legend()
 
# B.pl.figure()
# for h in hdelta_h:
#     h.plot(filled=False)
#%%HMS delta
B.pl.figure()
hdelta_pm800_histo_all.plot(color='#AF58BA',alpha=0.4,label='pm_800')
hdelta_pm580_histo_all.plot(color='#009ADE',alpha=0.3,label='pm_580') 
hdelta_pm900_histo_all.plot(color='#FFC61E',alpha=0.5,label='pm_900')
hdelta_pm120_histo_all.plot(color='#FF1F5B',alpha=0.5,label='pm_120')


for h in range(len(hdelta_hc_histo)):
    hdelta_hc_histo[h].plot(color=colors[h],alpha=0.3)
    
B.pl.title('HMS $\delta$')
B.pl.xlabel('')
B.pl.ylabel('')    
B.pl.legend()