#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 10:42:34 2025

@author: gvill

Investigating starttime issue
"""

import LT.box as B
import data_init as D
import cut_handler as C

#%%
br_sel = ['P.hod.goodstarttime','H.hod.goodstarttime',
          'P.hod.starttime','H.hod.starttime',
          'T.coin.pEDTM_tdcTimeRaw',
          'P.hod.TimeHist_StartTime_Sigma','P.hod.TimeHist_StartTime_Hits',
          'P.hod.TimeHist_FpTime_Sigma','P.hod.TimeHist_FpTime_Hits',
          'H.hod.TimeHist_StartTime_Sigma','H.hod.TimeHist_StartTime_Hits',
          'H.hod.TimeHist_FpTime_Sigma','H.hod.TimeHist_FpTime_Hits',
          'P.cal.etottracknorm']

RUNS = [20871,20956,21048,21102]  # pm120,pm800,pm900,pm580

ptof10_DIR = "/media/gvill/Gema's T7/ROOTfiles/starttime_issue_testing/"
ptof2_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

ptof10_RUN = D.DATA_INIT(data_type='deut23_data', run=RUNS, 
                    select_branches={'T':br_sel},ROOTfiles_path=ptof10_DIR)

ptof2_RUN = D.DATA_INIT(data_type='deut23_data', run=RUNS, 
                    select_branches={'T':br_sel},ROOTfiles_path=ptof2_DIR)

data_obj = [ptof10_RUN,ptof2_RUN]
#%% cuts
cuts_list = [C.noEDTM]

for d in data_obj:
    cuts_to_apply = {}
    for r in d.many:
        c_list = []
        print(f'Cuts for Run {r}')
        for cut in cuts_list:
            cut.init()
            br = d.Branches[r][C.HCANA_names[cut.name]]
            cut_array = cut(br)
            cut.stats()
            c_list.append(cut_array)
        
        cuts_to_apply[r] = c_list    
        
    all_cuts = {}
    for r in d.many:
        all_cuts_arr = cuts_to_apply[r][0]
        for arr in cuts_to_apply[r]:    
            all_cuts_arr = all_cuts_arr & arr    
        
        all_cuts[r] = all_cuts_arr
    
    for i in br_sel:
        for r in d.many:
            v = d.Branches[r][i]
            v_cut = v[all_cuts[r]]
            
            d.Branches[r].update({f'{i}_cut':v_cut})

#%% plot
i=1
for d in data_obj:
    if i==1:
        ptof='ptof_tolerance = 10'
    elif i==2:
        ptof='ptof_tolerance = 2'
    i+=1
    for r in d.many:
        # SHMS starttime var
        pgoodstarttime = d.Branches[r]['P.hod.goodstarttime']
        pstarttime = d.Branches[r]['P.hod.starttime']
        pstarttime_sigma = d.Branches[r]['P.hod.TimeHist_StartTime_Sigma']
        pstarttime_hits = d.Branches[r]['P.hod.TimeHist_StartTime_Hits']
        pfptime_sigma = d.Branches[r]['P.hod.TimeHist_FpTime_Sigma']
        pfptime_hits = d.Branches[r]['P.hod.TimeHist_FpTime_Hits']
        
        pgoodstarttime_h = B.histo(pgoodstarttime,range=(-0.5,1.5),bins=100,
                                   title='pgoodstarttime',xlabel='',ylabel='')
        pstarttime_h = B.histo(pstarttime,range=(-1100,300),bins=100,
                               title='pstarttime',xlabel='',ylabel='')
        pstarttime_sigma_h = B.histo(pstarttime_sigma,range=(0,5),bins=100,
                                     title='pstarttime_Sigma',xlabel='',ylabel='')
        pstarttime_hits_h = B.histo(pstarttime_hits,range=(0,20),bins=100,
                                    title='pstarttime_Hits',xlabel='',ylabel='')
        pfptime_sigma_h = B.histo(pfptime_sigma,range=(0,5),bins=100,
                                  title='pfptime_Sigma',xlabel='',ylabel='')
        pfptime_hits_h = B.histo(pfptime_hits,range=(0,20),bins=100,
                                 title='pfptime_Hits',xlabel='',ylabel='')
        
        fig, ax = B.pl.subplots(2,3,figsize=(15,8),
                                constrained_layout=True)
        B.pl.suptitle(f'Run {r} SHMS StartTime no Cuts\n{ptof}',fontsize=14)
        ax= B.pl.subplot(2,3,1)
        pgoodstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,2)
        pstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,3)
        pstarttime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,4)
        pstarttime_hits_h.plot()
        
        ax= B.pl.subplot(2,3,5)
        pfptime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,6)
        pfptime_hits_h.plot()
        
        print(f'Stats for Run {r}')
        print('pgoodstarttime total events = ', pgoodstarttime.size)
        print('pgoodstarttime == 0 events = ', (pgoodstarttime == 0).sum(),' (',
              ((pgoodstarttime == 0).sum()/pgoodstarttime.size)*100,'%)')
        print('pgoodstarttime == 1 events = ', (pgoodstarttime == 1).sum(),' (',
                ((pgoodstarttime == 1).sum()/pgoodstarttime.size)*100,'%)')
#%%
i=1
for d in data_obj:
    if i==1:
        ptof='ptof_tolerance = 10'
    elif i==2:
        ptof='ptof_tolerance = 2'
    i+=1
    for r in d.many:
        # SHMS starttime var
        pgoodstarttime = d.Branches[r]['P.hod.goodstarttime_cut']
        pstarttime = d.Branches[r]['P.hod.starttime_cut']
        pstarttime_sigma = d.Branches[r]['P.hod.TimeHist_StartTime_Sigma_cut']
        pstarttime_hits = d.Branches[r]['P.hod.TimeHist_StartTime_Hits_cut']
        pfptime_sigma = d.Branches[r]['P.hod.TimeHist_FpTime_Sigma_cut']
        pfptime_hits = d.Branches[r]['P.hod.TimeHist_FpTime_Hits_cut']
        
        pgoodstarttime_h = B.histo(pgoodstarttime,range=(-0.5,1.5),bins=100,
                                   title='pgoodstarttime',xlabel='',ylabel='')
        pstarttime_h = B.histo(pstarttime,range=(-1100,300),bins=100,
                               title='pstarttime',xlabel='',ylabel='')
        pstarttime_sigma_h = B.histo(pstarttime_sigma,range=(0,5),bins=100,
                                     title='pstarttime_Sigma',xlabel='',ylabel='')
        pstarttime_hits_h = B.histo(pstarttime_hits,range=(0,20),bins=100,
                                    title='pstarttime_Hits',xlabel='',ylabel='')
        pfptime_sigma_h = B.histo(pfptime_sigma,range=(0,5),bins=100,
                                  title='pfptime_Sigma',xlabel='',ylabel='')
        pfptime_hits_h = B.histo(pfptime_hits,range=(0,20),bins=100,
                                 title='pfptime_Hits',xlabel='',ylabel='')
        
        fig, ax = B.pl.subplots(2,3,figsize=(15,8),
                                constrained_layout=True)
        B.pl.suptitle(f'Run {r} SHMS StartTime noEDTM cut\n{ptof}',fontsize=14)
        ax= B.pl.subplot(2,3,1)
        pgoodstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,2)
        pstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,3)
        pstarttime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,4)
        pstarttime_hits_h.plot()
        
        ax= B.pl.subplot(2,3,5)
        pfptime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,6)
        pfptime_hits_h.plot()
        
        print(f'Stats for Run {r}')
        print('pgoodstarttime total events = ', pgoodstarttime.size)
        print('pgoodstarttime == 0 events = ', (pgoodstarttime == 0).sum(),' (',
              ((pgoodstarttime == 0).sum()/pgoodstarttime.size)*100,'%)')
        print('pgoodstarttime == 1 events = ', (pgoodstarttime == 1).sum(),' (',
                ((pgoodstarttime == 1).sum()/pgoodstarttime.size)*100,'%)')
        
        
#%% plot
i=1
for d in data_obj:
    if i==1:
        ptof='ptof_tolerance = 10'
    elif i==2:
        ptof='ptof_tolerance = 2'
    i+=1
    for r in d.many:
        # HMS starttime var
        hgoodstarttime = d.Branches[r]['H.hod.goodstarttime']
        hstarttime = d.Branches[r]['H.hod.starttime']
        hstarttime_sigma = d.Branches[r]['H.hod.TimeHist_StartTime_Sigma']
        hstarttime_hits = d.Branches[r]['H.hod.TimeHist_StartTime_Hits']
        hfptime_sigma = d.Branches[r]['H.hod.TimeHist_FpTime_Sigma']
        hfptime_hits = d.Branches[r]['H.hod.TimeHist_FpTime_Hits']
        
        hgoodstarttime_h = B.histo(hgoodstarttime,range=(-0.5,1.5),bins=100,
                                   title='hgoodstarttime',xlabel='',ylabel='')
        hstarttime_h = B.histo(hstarttime,range=(-1100,300),bins=100,
                               title='hstarttime',xlabel='',ylabel='')
        hstarttime_sigma_h = B.histo(hstarttime_sigma,range=(0,5),bins=100,
                                     title='hstarttime_Sigma',xlabel='',ylabel='')
        hstarttime_hits_h = B.histo(hstarttime_hits,range=(0,20),bins=100,
                                    title='hstarttime_Hits',xlabel='',ylabel='')
        hfptime_sigma_h = B.histo(hfptime_sigma,range=(0,5),bins=100,
                                  title='hfptime_Sigma',xlabel='',ylabel='')
        hfptime_hits_h = B.histo(hfptime_hits,range=(0,20),bins=100,
                                 title='hfptime_Hits',xlabel='',ylabel='')
        
        fig, ax = B.pl.subplots(2,3,figsize=(15,8),
                                constrained_layout=True)
        B.pl.suptitle(f'Run {r} HMS StartTime no Cuts\n{ptof}',fontsize=14)
        ax= B.pl.subplot(2,3,1)
        hgoodstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,2)
        hstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,3)
        hstarttime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,4)
        hstarttime_hits_h.plot()
        
        ax= B.pl.subplot(2,3,5)
        hfptime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,6)
        hfptime_hits_h.plot()
        
        print(f'Stats for Run {r}')
        print('hgoodstarttime total events = ', hgoodstarttime.size)
        print('hgoodstarttime == 0 events = ', (hgoodstarttime == 0).sum(),' (',
              ((hgoodstarttime == 0).sum()/hgoodstarttime.size)*100,'%)')
        print('hgoodstarttime == 1 events = ', (hgoodstarttime == 1).sum(),' (',
                ((hgoodstarttime == 1).sum()/hgoodstarttime.size)*100,'%)')
#%%
i=1
for d in data_obj:
    if i==1:
        ptof='ptof_tolerance = 10'
    elif i==2:
        ptof='ptof_tolerance = 2'
    i+=1
    for r in d.many:
        # SHMS starttime var
        hgoodstarttime = d.Branches[r]['H.hod.goodstarttime_cut']
        hstarttime = d.Branches[r]['H.hod.starttime_cut']
        hstarttime_sigma = d.Branches[r]['H.hod.TimeHist_StartTime_Sigma_cut']
        hstarttime_hits = d.Branches[r]['H.hod.TimeHist_StartTime_Hits_cut']
        hfptime_sigma = d.Branches[r]['H.hod.TimeHist_FpTime_Sigma_cut']
        hfptime_hits = d.Branches[r]['H.hod.TimeHist_FpTime_Hits_cut']
        
        hgoodstarttime_h = B.histo(hgoodstarttime,range=(-0.5,1.5),bins=100,
                                   title='hgoodstarttime',xlabel='',ylabel='')
        hstarttime_h = B.histo(hstarttime,range=(-1100,300),bins=100,
                               title='hstarttime',xlabel='',ylabel='')
        hstarttime_sigma_h = B.histo(hstarttime_sigma,range=(0,5),bins=100,
                                     title='hstarttime_Sigma',xlabel='',ylabel='')
        hstarttime_hits_h = B.histo(hstarttime_hits,range=(0,20),bins=100,
                                    title='hstarttime_Hits',xlabel='',ylabel='')
        hfptime_sigma_h = B.histo(hfptime_sigma,range=(0,5),bins=100,
                                  title='hfptime_Sigma',xlabel='',ylabel='')
        hfptime_hits_h = B.histo(hfptime_hits,range=(0,20),bins=100,
                                 title='hfptime_Hits',xlabel='',ylabel='')
        
        fig, ax = B.pl.subplots(2,3,figsize=(15,8),
                                constrained_layout=True)
        B.pl.suptitle(f'Run {r} HMS StartTime noEDTM cut\n{ptof}',fontsize=14)
        ax= B.pl.subplot(2,3,1)
        hgoodstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,2)
        hstarttime_h.plot()
        
        ax= B.pl.subplot(2,3,3)
        hstarttime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,4)
        hstarttime_hits_h.plot()
        
        ax= B.pl.subplot(2,3,5)
        hfptime_sigma_h.plot()
        
        ax= B.pl.subplot(2,3,6)
        hfptime_hits_h.plot()
        
        print(f'Stats for Run {r}')
        print('hgoodstarttime total events = ', hgoodstarttime.size)
        print('hgoodstarttime == 0 events = ', (hgoodstarttime == 0).sum(),' (',
              ((hgoodstarttime == 0).sum()/hgoodstarttime.size)*100,'%)')
        print('hgoodstarttime == 1 events = ', (hgoodstarttime == 1).sum(),' (',
                ((hgoodstarttime == 1).sum()/hgoodstarttime.size)*100,'%)')



