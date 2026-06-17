#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 10:37:26 2025

@author: gvill

Extract heep yield ratios
"""

import data_init as D
import LT.box as B
import cut_handler as C
import database_operations as db

import numpy as np
import pickle as P

#%% helper functions
def set_value(inpt,value=np.nan,set_value_to=0.):
    is_value = np.where(inpt == value)
    for i in is_value[0]:
        inpt[i] = set_value_to
    print(f'Set {is_value[0].size} {value} to {set_value_to}\n')
    return

# To get rid of unhelpful groupings (lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0]

#%% analysis functions

def apply_evt_sel_cuts(data_obj,is_SIMC=False,show_stats=False,save_report=False):
    if is_SIMC:
        try:
            if type(data_obj.many) == list:
                W = {}
                Q2 = {}
                WEIGHTS = {}
                WEIGHTS_PS = {}
                for m in data_obj.many:
                    W[m] = data_obj.Branches[m]['W']
                    Q2[m] = data_obj.Branches[m]['Q2']
                    WEIGHTS[m] = D.calc_weights(data_obj.Branches[m])
                    WEIGHTS_PS[m] = D.calc_weights_PS(data_obj.Branches[m])
                
                cuts_applied_list = []
                
                cuts_list = [C.hms_delta,C.shms_delta,C.Em_cut_hc]
                cuts_applied_list = cuts_list.copy()
                
                hcoll_cut = C.coll_cut(data_obj,spec='HMS',is_SIMC=True,many=True)
                hcoll_cut_arrays = hcoll_cut()
                cuts_applied_list.append(hcoll_cut)
                
                zt_cut = C.ztar_cut(data_obj,is_SIMC=True,many=True)
                zt_cut_arrays = zt_cut()
                cuts_applied_list.append(zt_cut)
    
                cuts_to_apply = {}
                for m in data_obj.many:
                   cuts_to_apply_list = []
                   for cut in cuts_list:
                       br = data_obj.Branches[m][C.SIMC_names[cut.name]]
                       cut_array = cut(br)
                       # cut.stats()
                       
                       cuts_to_apply_list.append(cut_array)
                   # add special cut arrays to clist
                   cuts_to_apply_list.append(hcoll_cut_arrays[m])
                   cuts_to_apply_list.append(zt_cut_arrays[m])
                   
                   cuts_to_apply[m] = cuts_to_apply_list
                
                all_cuts_sim = {}
                for m in data_obj.many:
                    all_cuts_sim[m] = cuts_to_apply[m][0]
                    for arr in cuts_to_apply[m]:
                        all_cuts_sim[m] = all_cuts_sim[m] & arr 
               
                for m in data_obj.many:
                   W_cut = W[m][all_cuts_sim[m]]
                   Q2_cut = Q2[m][all_cuts_sim[m]]
                   WEIGHTS_cut = WEIGHTS[m][all_cuts_sim[m]]
                   WEIGHTS_PS_cut = WEIGHTS_PS[m][all_cuts_sim[m]]
           
                   data_obj.Branches[m].update({'W_cut':W_cut,
                                                'Q2_cut':Q2_cut,
                                                'WEIGHTS_cut':WEIGHTS_cut,
                                                'WEIGHTS_PS_cut':WEIGHTS_PS_cut})
                
                print('SIMC data object updated with the branches:',
                      'W_cut, Q2_cut,\n',
                      'WEIGHTS_cut, WEIGHTS_PS_cut\n')
                print('Applied the following cuts:\n',
                      [c for c in cuts_applied_list])
            
                if show_stats:
                    print('Cuts STATS:\n',
                          [c.stats() for c in cuts_applied_list])
                    
                if save_report:
                    filename = 'yield_n_xsec/current_cuts_SIMC.txt'
                    with open(filename,'w') as f:
                        f.write('CUTS REPORT\n')
                    a = [c.report(filename) for c in cuts_applied_list]
            else:
                raise(TypeError)
                
        except TypeError:
            W = data_obj.Branches['W']
            Q2 = data_obj.Branches['Q2']
            WEIGHTS = D.calc_weights(data_obj.Branches)
            WEIGHTS_PS = D.calc_weights_PS(data_obj.Branches)
            
            cuts_applied_list = []
            
            cuts_list = [C.hms_delta,C.shms_delta,C.Em_cut_hc,C.W_cut]
            cuts_applied_list = cuts_list.copy()
            
            hcoll_cut = C.coll_cut(data_obj,spec='HMS',is_SIMC=True)
            hcoll_cut_arrays = hcoll_cut()
            cuts_applied_list.append(hcoll_cut)
            
            zt_cut = C.ztar_cut(data_obj,is_SIMC=True)
            zt_cut_arrays = zt_cut()
            cuts_applied_list.append(zt_cut)
            
            cuts_to_apply_list = []
            for cut in cuts_list:
                br = data_obj.Branches[C.SIMC_names[cut.name]]
                cut_array = cut(br)
                # cut.stats()
                
                cuts_to_apply_list.append(cut_array)
            # add special cut arrays to clist
            cuts_to_apply_list.append(hcoll_cut_arrays)
            cuts_to_apply_list.append(zt_cut_arrays)
            
            all_cuts_sim = cuts_to_apply_list[0]
            for arr in cuts_to_apply_list:    
                all_cuts_sim = all_cuts_sim & arr  

            W_cut = W[all_cuts_sim]
            Q2_cut = Q2[all_cuts_sim]
            WEIGHTS_cut = WEIGHTS[all_cuts_sim]
            WEIGHTS_PS_cut = WEIGHTS_PS[all_cuts_sim]

            data_obj.Branches.update({'W_cut':W_cut,
                                      'Q2_cut':Q2_cut,
                                      'WEIGHTS_cut':WEIGHTS_cut,
                                      'WEIGHTS_PS_cut':WEIGHTS_PS_cut})
           
            print('SIMC data object updated with the branches:',
                  'W_cut, Q2_cut,\n',
                  'WEIGHTS_cut, WEIGHTS_PS_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_applied_list])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_applied_list])
                
            if save_report:
                filename = 'yield_n_xsec/current_cuts_SIMC.txt'
                with open(filename,'w') as f:
                    f.write('CUTS REPORT\n')
                a = [c.report(filename) for c in cuts_applied_list]

    else:
        try:
            W = {}
            Q2 = {}

            for m in data_obj.many:
                W[m] = data_obj.Branches[m]['P.kin.primary.W']
                Q2[m] = data_obj.Branches[m]['P.kin.primary.Q2']    

            ## making the cut lists
            # first select the cuts to apply, defined already in cut_handler
            cuts_applied_list = []
            
            cuts_list = [C.hms_delta, C.shms_delta, C.shms_calPID, 
                         C.Em_cut_hc, C.noEDTM]
            cuts_applied_list = cuts_list.copy()
            
            # initialize special cut classes: coll_cut, current_cut, ztar_cut,
            # CTime_cut 
            hcoll_cut = C.coll_cut(data_obj,spec='HMS',many=True)
            hcoll_cut_arrays = hcoll_cut()
            cuts_applied_list.append(hcoll_cut)
            
            curr_cut = C.current_cut(data_obj,current='BCM4A',many=True)
            curr_cut_arrays = curr_cut()
            cuts_applied_list.append(curr_cut)
            
            zt_cut = C.ztar_cut(data_obj,many=True)
            zt_cut_arrays = zt_cut()
            cuts_applied_list.append(zt_cut)
            
            ct_cut = C.CTime_cut(data_obj,many=True)
            ct_cut_arrays = ct_cut()
            cuts_applied_list.append(ct_cut)
            
            #then make the cut arrays by getting the desired array from the DATA_INIT 
            # object and print stats for said cut, the result is a list of boolean arrays
            # for each cut, e.g. a Q2 cut will make a cut on the 'P.kin.primary.Q2' 
            # variable and the result will be a boolean array
            # note you will have a list of boolean arrays PER RUN, this is necessary to 
            # preserve the length of the arrays, since each run has a different # of array
            # elements.

            cuts_to_apply = {}
            for m in data_obj.many:

                clist = []
                for cut in cuts_list:
                    br = data_obj.Branches[m][C.HCANA_names[cut.name]]
                    cut_array = cut(br)
                    # cut.stats()
                    
                    clist.append(cut_array)
                
                # add special cut arrays to clist
                clist.append(hcoll_cut_arrays[m])
                clist.append(curr_cut_arrays[m])
                clist.append(zt_cut_arrays[m])
                clist.append(ct_cut_arrays[m])
                
                # each run has a clist
                cuts_to_apply[m] = clist

            # to apply multiple cuts efficiently I combine all boolean arrays into one
            # again there will be a total array for each run to preserve array length
            # WARNING: if you try to make a cut by array slicing 
            #  e.g. cut_arr = arr[bool_arr] 
            # you will get a complaint if the boolean array is of a different size than 
            # the original array

            all_cuts = {}
            for m in data_obj.many:
                all_cuts[m] = cuts_to_apply[m][0]
                for arr in cuts_to_apply[m]:    
                    all_cuts[m] = all_cuts[m] & arr  

            # now apply the cuts and store the cut arrays back in the DATA_INIT obj
            for m in data_obj.many:
                W_cut = W[m][all_cuts[m]]
                Q2_cut = Q2[m][all_cuts[m]]
                
                data_obj.Branches[m].update({'P.kin.primary.W_cut':W_cut,
                                               'P.kin.primary.Q2_cut':Q2_cut})
            
            print('Data object updated with the branches:\n',
                  'P.kin.primary.W_cut,',
                  'P.kin.primary.Q2_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_applied_list])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_applied_list])
                
            if save_report:
                filename = 'yield_n_xsec/current_cuts_SIMC.txt'
                with open(filename,'w') as f:
                    f.write('CUTS REPORT\n')
                a = [c.report(filename) for c in cuts_applied_list]
                
        except TypeError:
            
            W = data_obj.Branches['P.kin.primary.W']
            Q2 = data_obj.Branches['P.kin.primary.Q2']
    

            ## making the cut lists
            # first select the cuts to apply, defined already in cut_handler
            cuts_applied_list = []
            
            cuts_list = [C.hms_delta, C.shms_delta, C.shms_calPID, 
                         C.Em_cut_hc, C.noEDTM, C.W_cut]
            cuts_applied_list = cuts_list.copy()
            # initialize special cut classes: coll_cut, current_cut, ztar_cut,
            # CTime_cut 
            hcoll_cut = C.coll_cut(data_obj,spec='HMS')
            hcoll_cut_arrays = hcoll_cut()
            cuts_applied_list.append(hcoll_cut)
            
            curr_cut = C.current_cut(data_obj,current='BCM4A')
            curr_cut_arrays = curr_cut()
            cuts_applied_list.append(curr_cut)
            
            zt_cut = C.ztar_cut(data_obj)
            zt_cut_arrays = zt_cut()
            cuts_applied_list.append(zt_cut)
            
            ct_cut = C.CTime_cut(data_obj)
            ct_cut_arrays = ct_cut() 
            cuts_applied_list.append(ct_cut)
             
            #then make the cut arrays by getting the desired array from the DATA_INIT 
            # object and print stats for said cut, the result is a list of boolean arrays
            # for each cut, e.g. a Q2 cut will make a cut on the 'P.kin.primary.Q2' 
            # variable and the result will be a boolean array
            # note you will have a list of boolean arrays PER RUN, this is necessary to 
            # preserve the length of the arrays, since each run has a different # of array
            # elements.

            cuts_to_apply_list = []
            for cut in cuts_list:
                br = data_obj.Branches[C.HCANA_names[cut.name]]
                cut_array = cut(br)
                # cut.stats()
                
                cuts_to_apply_list.append(cut_array)
                
            # add special cut arrays to cuts_to_apply_list
            cuts_to_apply_list.append(hcoll_cut_arrays)
            cuts_to_apply_list.append(curr_cut_arrays)
            cuts_to_apply_list.append(zt_cut_arrays)
            cuts_to_apply_list.append(ct_cut_arrays)
                
            # to apply multiple cuts efficiently I combine all boolean arrays into one
            # again there will be a total array for each run to preserve array length
            # WARNING: if you try to make a cut by array slicing 
            #  e.g. cut_arr = arr[bool_arr] 
            # you will get a complaint if the boolean array is of a different size than 
            # the original array

            all_cuts = cuts_to_apply_list[0]
            for arr in cuts_to_apply_list:    
                all_cuts = all_cuts & arr  

            # now apply the cuts and store the cut arrays back in the DATA_INIT obj
  
            W_cut = W[all_cuts]
            Q2_cut = Q2[all_cuts]
            
            data_obj.Branches.update({'P.kin.primary.W_cut':W_cut,
                                           'P.kin.primary.Q2_cut':Q2_cut})
        
            print('Data object updated with the branches:\n',
                  'P.kin.primary.W_cut, P.kin.primary.Q2_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_applied_list])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_applied_list])
                
            if save_report:
                filename = 'yield_n_xsec/current_cuts_SIMC.txt'
                with open(filename,'w') as f:
                    f.write('CUTS REPORT\n')
                a = [c.report(filename) for c in cuts_applied_list]

def get_yield(h2,norm,bin_range,x_min,x_max):
    counts = []
    for nx in bin_range[:]:
        h = norm*h2.project_x(bins = [nx])
        
        counts.append(h.sum(x_min,x_max))
    
    counts = np.array(counts)
       
    y_values = h2.y_bin_center[bin_range]
    yield_values = counts[:,0]
    yield_errors = counts[:,1]
    
    return (y_values,yield_values,yield_errors)

#%% analysis classes

class yield_plots:
    def __init__(self, data_obj, use_cuts = False, data_type = ''):
        d = data_obj.Branches
        if type(data_obj.many) == list:
            self.many = data_obj.many
        else:
            self.many = False
        self.histos_1D = {'W':None, 'Q2':None}
        self.histos_2D = {'W_vs_Q2':None}
        
        self.set_hist_lim()
        
        if data_type == 'deut23_data':
            try:
                self.W = {}
                self.Q2 = {}
                for m in self.many:
                    if use_cuts:
                        self.W[m] = d[m]['P.kin.primary.W_cut']
                        self.Q2[m] = d[m]['P.kin.primary.Q2_cut']
                    else:    
                        self.W[m] = d[m]['P.kin.primary.W']
                        self.Q2[m] = d[m]['P.kin.primary.Q2']

            except TypeError:
                if use_cuts:
                    self.W = d['P.kin.primary.W_cut']
                    self.Q2 = d['P.kin.primary.Q2_cut']
                else:
                    self.W = d['P.kin.primary.W']
                    self.Q2 = d['P.kin.primary.Q2']
        
        elif data_type == 'SIMC':
            try:
                self.W = {}
                self.Q2 = {}
                self.WEIGHTS = {}
                self.WEIGHTS_PS = {}
                for m in self.many:
                    if use_cuts:
                        self.W[m] = d[m]['W_cut']
                        self.Q2[m] = d[m]['Q2_cut']
                        self.WEIGHTS[m] = d[m]['WEIGHTS_cut']
                        self.WEIGHTS_PS[m] = d[m]['WEIGHTS_PS_cut']
                        
                    else:    
                        self.W[m] = d[m]['W']
                        self.Q2[m] = d[m]['Q2']
                        self.WEIGHTS[m] = D.calc_weights(d[m])
                        self.WEIGHTS_PS[m] = D.calc_weights_PS(d[m])
                    
            except TypeError:               
                if use_cuts:
                    self.W = d['W_cut']
                    self.Q2 = d['Q2_cut']
                    self.WEIGHTS = d['WEIGHTS_cut']
                    self.WEIGHTS_PS = d['WEIGHTS_PS_cut']
                    
                else:    
                    self.W = d['W']
                    self.Q2 = d['Q2']
                    self.WEIGHTS = D.calc_weights(d)
                    self.WEIGHTS_PS = D.calc_weights_PS(d)
                    
        else:
            print('No data type chosen.\n',
                  'Available types: "deut23_data", "SIMC"')
    ###
    # method to create desired histograms, note for SIMC histograms will be 
    # weighted if the flag weights is set to True.
    ###    
    def make_histos(self, weights = False):
        histos_to_plot = {'W':[self.W,self.W_range,self.W_bins],
                          'Q2':[self.Q2,self.Q2_range,self.Q2_bins]}
        
        histos2D_to_plot = {'W_vs_Q2':[self.Q2,self.W,
                                      [self.Q2_range,self.Q2_bins],
                                      [self.W_range,self.W_bins]]}
                            

        # make 1D histograms and store them in self.histos_1D
        try:
            for v in histos_to_plot:
                h = {}
                if type(self.many) is list:
                    for m in self.many:
                        histo = histos_to_plot[v][0][m]
                        ran = histos_to_plot[v][1]
                        nbins = histos_to_plot[v][2]
                        
                        if weights:
                            h[m] = B.histo(histo,range = ran, bins = nbins,
                                           weights = self.WEIGHTS[m],
                                           calc_w2=True)
                        else:
                            h[m] = B.histo(histo,range = ran, bins = nbins)                
                    self.histos_1D[v] = h
                else:
                    raise(TypeError)
                
        except TypeError:
            for v in histos_to_plot:
                histo = histos_to_plot[v][0]
                ran = histos_to_plot[v][1]
                nbins = histos_to_plot[v][2]
                
                if weights:
                    h = B.histo(histo,range = ran, bins = nbins,
                                   weights = self.WEIGHTS,
                                   calc_w2=True)
                else:    
                    h = B.histo(histo,range = ran, bins = nbins)
                    
                self.histos_1D[v] = h
            
                
        # make 2D histograms and store them in self.histos_2D
        try:
            for v in histos2D_to_plot:
                h = {}
                if type(self.many) is list:

                    for m in self.many:
                        hx = histos2D_to_plot[v][0][m]
                        hy = histos2D_to_plot[v][1][m]
                        ranx = histos2D_to_plot[v][2][0]
                        rany = histos2D_to_plot[v][3][0]
                        nxbins = histos2D_to_plot[v][2][1]
                        nybins = histos2D_to_plot[v][3][1]
                        
                        if weights:
                            h[m] = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins],
                                             weights = self.WEIGHTS[m],
                                             calc_w2=True)
                        else:
                            h[m] = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins])               
                    self.histos_2D[v] = h
                else:
                    raise(TypeError)
                
        except TypeError:
            for v in histos2D_to_plot:
                hx = histos2D_to_plot[v][0]
                hy = histos2D_to_plot[v][1]
                ranx = histos2D_to_plot[v][2][0]
                rany = histos2D_to_plot[v][3][0]
                nxbins = histos2D_to_plot[v][2][1]
                nybins = histos2D_to_plot[v][3][1]
                
                if weights:
                    h = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins],
                                     weights = self.WEIGHTS,
                                     calc_w2=True)
                else:
                    h = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins])
                
                self.histos_2D[v] = h
    
    def make_PS_histos(self):
        histos_to_plot = {'W_PS':[self.W,self.W_range,self.W_bins],
                          'Q2_PS':[self.Q2,self.Q2_range,self.Q2_bins]}
        
        histos2D_to_plot = {'W_vs_Q2_PS':[self.Q2,self.W,
                                          [self.Q2_range,self.Q2_bins],
                                          [self.W_range,self.W_bins]]}
        
        try:
            h = {}
            for m in self.many:
                # make 1D phase space histograms
                for v in histos_to_plot:
                    histo = histos_to_plot[v][0][m]
                    ran = histos_to_plot[v][1]
                    nbins = histos_to_plot[v][2]
                    
                    h[m] = B.histo(histo,range = ran, bins = nbins,
                                   weights = self.WEIGHTS_PS[m],
                                   calc_w2=True)
                   
                self.histos_1D[v] = h
            
            h = {}
            for m in self.many:   
                # make 2D phase space histograms
                for v in histos2D_to_plot:
                    hx = histos2D_to_plot[v][0][m]
                    hy = histos2D_to_plot[v][1][m]
                    ranx = histos2D_to_plot[v][2][0]
                    rany = histos2D_to_plot[v][3][0]
                    nxbins = histos2D_to_plot[v][2][1]
                    nybins = histos2D_to_plot[v][3][1]
                    
                    h[m] = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins],
                                     weights = self.WEIGHTS_PS[m],
                                     calc_w2=True)
                     
                self.histos_2D[v] = h
            
        except TypeError:
            # make 1D phase space histograms
            for v in histos_to_plot:
                histo = histos_to_plot[v][0]
                ran = histos_to_plot[v][1]
                nbins = histos_to_plot[v][2]
                
                h = B.histo(histo,range = ran, bins = nbins,
                               weights = self.WEIGHTS_PS,
                               calc_w2=True)
               
                self.histos_1D[v] = h
    
            # make 2D phase space histograms
            for v in histos2D_to_plot:
                hx = histos2D_to_plot[v][0]
                hy = histos2D_to_plot[v][1]
                ranx = histos2D_to_plot[v][2][0]
                rany = histos2D_to_plot[v][3][0]
                nxbins = histos2D_to_plot[v][2][1]
                nybins = histos2D_to_plot[v][3][1]
                
                h = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins],
                                 weights = self.WEIGHTS_PS,
                                 calc_w2=True)
                 
                self.histos_2D[v] = h
    
    def set_hist_lim(self,Q2_range=(3.2,4.4),Q2_bins=50,W_range=(0.8,1.1),
                     W_bins=100):
        self.Q2_range = Q2_range
        self.Q2_bins = Q2_bins
        self.W_range = W_range
        self.W_bins = W_bins
    
    def plot_histos(self,histo=[]):
        all_histos = {}
        all_histos.update(self.histos_1D)
        all_histos.update(self.histos_2D)
        
        if histo:
            try:
                for h in histo:
                    # print('Plotting ', h)
                    for m in self.many:
                        B.pl.figure()
                        p = all_histos[h][m]
                        p.plot()
                        B.pl.title(f'{h} {m}')
                        B.pl.xlabel('')
                        B.pl.ylabel('')
                                     
            except TypeError:
                # print('Plotting ', h)  
                B.pl.figure()
                p = all_histos[h]
                p.plot()
                B.pl.title(h)
                B.pl.xlabel('')
                B.pl.ylabel('')
                
        else:
            try: 
                for h in all_histos:
                    # print('Plotting ', h)
                    for m in self.many:  
                        B.pl.figure()
                        p = all_histos[h][m]
                        p.plot()
                        B.pl.title(f'{h} {m}')
                        B.pl.xlabel('')
                        B.pl.ylabel('')
                        
            except TypeError:    
                for h in all_histos:
                    # print('Plotting ', h)  
                    B.pl.figure()
                    p = all_histos[h]
                    p.plot()
                    B.pl.title(h)
                    B.pl.xlabel('')
                    B.pl.ylabel('')


#%% select branches
T_sel = ['P.kin.primary.W','P.kin.primary.Q2','H.kin.secondary.emiss',
         'H.extcor.xsieve','H.extcor.ysieve','P.extcor.xsieve',
         'P.extcor.ysieve','CTime.epCoinTime_ROC2','H.gtr.dp',
         'P.gtr.dp','P.cal.etottracknorm','H.react.z',
         'P.react.z','P.cal.etotnorm','H.cal.etottracknorm','H.cal.etotnorm',
         'T.coin.pEDTM_tdcTimeRaw','P.kin.primary.x_bj']

br_sel_sim = ['Em','W','Q2','e_delta','h_delta','Weight','Normfac',
              'e_xptar','e_yptar','h_xptar','h_yptar','tar_x',
              'h_zv','e_zv','h_ytar','e_ytar','sig','sigcc','xcoll','ycoll',
              'Jacobian_corr','probabs']             
    
TSP_sel = ['evNumber','P.BCM4A.scalerCharge','P.BCM4A.scalerChargeCut',
           'P.BCM4A.scalerCurrent','P.1MHz.scalerTime']

br_sel = {'T':T_sel,'TSP':TSP_sel}
t_sel = ['T','TSP']

#%% load runs with selected branches
DATA_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"
sett = input('Select heep delta setting:\n')

if sett == '-8':
    RUN = 20841
elif sett == '-4':
    RUN = 20846
elif sett == '0':
    RUN = 20851
elif sett == '+4':
    RUN = 20858
elif sett == '+8':
    RUN = 20861
elif sett == '+12':
    RUN = 20840

heep_dp = D.DATA_INIT(data_type='deut23_data',run=RUN,
                         select_branches=br_sel,
                         select_trees=t_sel,
                         ROOTfiles_path= DATA_DIR)

# load SIMC files
# directory paths for different SIMC files
SIMC_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/worksim" +\
                "/heep_alloffsets/"
SETTING = f'delta_scan_{sett}'

SIMC_rad = D.DATA_INIT(data_type='SIMC',setting = SETTING,
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='-')
SIMC_norad = D.DATA_INIT(data_type='SIMC',setting = SETTING,
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='norad')

# apply event selection cuts to data
apply_evt_sel_cuts(heep_dp)
apply_evt_sel_cuts(SIMC_norad,is_SIMC=True)
apply_evt_sel_cuts(SIMC_rad,is_SIMC=True)

# create necessary plots using helper class
data_plots_nocuts = yield_plots(heep_dp,use_cuts=False,data_type='deut23_data')
# data_plots = yield_plots(heep_dp,use_cuts=True,data_type='deut23_data')
data_plots_nocuts.set_hist_lim(Q2_range=(1,6),Q2_bins=100)
data_plots_nocuts.make_histos()

data_plots = yield_plots(heep_dp,use_cuts=True,data_type='deut23_data')
# data_plots = yield_plots(heep_dp,use_cuts=True,data_type='deut23_data')
data_plots.set_hist_lim(Q2_range=(1,6),Q2_bins=100)
data_plots.make_histos()

simc_norad_plots = yield_plots(SIMC_norad,use_cuts=True,data_type='SIMC')
simc_norad_plots.set_hist_lim(Q2_range=(1,6),Q2_bins=100)
simc_norad_plots.make_histos(weights=True)
simc_norad_plots.make_PS_histos()

simc_rad_plots = yield_plots(SIMC_rad,use_cuts=True,data_type='SIMC')
simc_rad_plots.set_hist_lim(Q2_range=(1,6),Q2_bins=100)
simc_rad_plots.make_histos(weights=True)
simc_rad_plots.make_PS_histos()

#%%
# data normalization factor
data_norm = 1/(D.get_eff_norm(RUN)*D.get_charge_norm(RUN))

# simc radiative correction factor
radcorr_fac = simc_norad_plots.histos_1D['W']/\
                simc_rad_plots.histos_1D['W']

# raw counts (no counts)
raw_counts_h = data_plots_nocuts.histos_1D['W']

# raw yield
yield_h = data_plots.histos_1D['W']

# normalized yield
yield_norm_h = data_norm*data_plots.histos_1D['W']

# SIMC phase space
yield_PS_h = simc_norad_plots.histos_1D['W_PS']

# data cross section
yield_radcorr_h = radcorr_fac*yield_norm_h

Xsec_h = yield_radcorr_h/yield_PS_h

# SIMC Normalized Yield w/ rad
SIMC_yield_norm_rad_h = simc_rad_plots.histos_1D['W']

# SIMC Normalized Yield no rad
SIMC_yield_norm_norad_h = simc_norad_plots.histos_1D['W']

# SIMC Cross Section
SIMC_Xsec_h = SIMC_yield_norm_norad_h/yield_PS_h

# yield ratio
Ratio_h = yield_radcorr_h/SIMC_yield_norm_norad_h

#%%
# save histograms to a pickle file 
# define dict of histograms to save
hist = {'Raw_counts': raw_counts_h,
        'heep_yield': yield_h,
        'heep_yield_norm': yield_norm_h,
        'heep_radcorr_factor': radcorr_fac,
        'heep_yield_PS': yield_PS_h,
        'heep_yield_radcorr': yield_radcorr_h,
        'heep_xsec': Xsec_h,
        'heep_SIMC_yield_norm_rad': SIMC_yield_norm_rad_h,
        'heep_SIMC_yield_norm_norad': SIMC_yield_norm_norad_h,
        'heep_SIMC_xsec': SIMC_Xsec_h,
        'heep_ratio': Ratio_h}

filename = f'yield_n_xsec/pickle/{SETTING}_histos.pcl'

with open(filename, 'wb') as f:
    P.dump(hist,f)
print(f'Saved {filename}\n')

#%%heep yield ratio
data_norm = 1/(D.get_eff_norm(RUN)*D.get_charge_norm(RUN))

data_yield_h = data_norm*data_plots.histos_2D['W_vs_Q2']

simc_yield_h = simc_rad_plots.histos_2D['W_vs_Q2']

data_yi_proj = data_yield_h.project_x(range=(0.85,1.05))

data_yield = data_yi_proj.sum()

simc_yi_proj = simc_yield_h.project_x(range=(0.85,1.05))

simc_yield = simc_yi_proj.sum()

yield_ratio = data_yield[0]/simc_yield[0]

yield_ratio_err = data_yield[1]/data_yield[0] + simc_yield[1]/simc_yield[0]

print('yield ratio = ', yield_ratio,'+/-',yield_ratio_err)

#%%
write_header = False
with open('./yield_n_xsec/heep_yield.txt','a') as f:
    if write_header:
        f.write('#! Run[i,0]/ Setting[s,1]/ Yield[f,2]/ Yield Error[f,3]/'+\
                ' SIMC Yield[f,4]/ SIMC Yield Error[f,5]/ Yield Ratio[f,6]/'+\
                    ' Yield Ratio Error[f,7]/\n')
    f.write(f'{RUN} {SETTING} {data_yield[0]:0.4f} {data_yield[1]:0.4f}'+\
            f' {simc_yield[0]:0.4f} {simc_yield[1]:0.4f} {yield_ratio:0.4f}'+\
                f' {yield_ratio_err:0.4f}\n')
#%% plot yield ratios

fname = './yield_n_xsec/heep_yield.txt'
file = B.get_file(fname)

runs = file['Run']
delta = [-8,-4,0,4,8,12]
yr = file['Yield Ratio']
yr_err = file['Yield Ratio Error']

B.pl.figure()
B.plot_exp(delta, yr, yr_err)
B.pl.ylim(0.5,1.5)
#%%
data_yield = {}
simc_yield = {}
yield_ratio = {}
yield_ratio_err = {}
for m in heep_dp_all.many:
    s = heep_dp.setting[m]
    
    data_norm = 1/(D.get_eff_norm(m)*D.get_charge_norm(m))
    # print('\n',data_norm,'\n')

    data_yield_h = data_norm*data_plots_all.histos_2D['W_vs_Q2'][m]

    simc_yield_h = simc_rad_plots_all.histos_2D['W_vs_Q2'][s]

    data_yi_proj = data_yield_h.project_x(range=(0.85,1.05))

    data_yield[m] = data_yi_proj.sum()

    simc_yi_proj = simc_yield_h.project_x(range=(0.85,1.05))

    simc_yield[s] = simc_yi_proj.sum()

    yield_ratio[m] = data_yield[m][0]/simc_yield[s][0]

    yield_ratio_err[m] = data_yield[m][1]/data_yield[m][0] + simc_yield[s][1]/simc_yield[s][0]

    print(f'yield ratio run {m} = ', yield_ratio[m],'+/-',yield_ratio_err[m],'\n')
#%%
runs = [20840,20841,20846,20851,20858,20861,20868,20869]

ps1x = []
q2_means = []
for r in runs:
    ps1x.append(get_list(db.retrieve('deuteron_db.db', 
            ' pS1X_scl_rates','RUN_LIST_UPDATED', where = f"run=\'{r}\'")))
    
    h = data_plots.histos_1D['Q2'][r]
    h.fit()
    q2_means.append(h.mean.value)
    
yr = [yield_ratio[r] for r in yield_ratio]
yr_e = [yield_ratio_err[r] for r in yield_ratio_err]

#%%
B.pl.figure()
B.plot_exp(runs, yr, yr_e)
xmin, xmax = B.pl.xlim()
B.pl.hlines([1,1.1,0.9,1.05,0.95],xmin,xmax,linestyle=['-','--','--','--','--'],
            color=['lightgray','r','r','g','g'])
B.pl.ylim(0.6,1.2)
B.pl.xlabel('Run #')
B.pl.ylabel('$Y_{exp}/Y_{SIMC}$')

B.pl.figure()
B.plot_exp(ps1x, yr, yr_e)
xmin, xmax = B.pl.xlim()
B.pl.hlines([1,1.1,0.9,1.05,0.95],xmin,xmax,linestyle=['-','--','--','--','--'],
            color=['lightgray','r','r','g','g'])
B.pl.ylim(0.6,1.2)
B.pl.xlabel('SHMS S1X Rate')
B.pl.ylabel('$Y_{exp}/Y_{SIMC}$')
# ax = B.pl.gca()
# twin1 = ax.twiny()
B.pl.figure()
B.plot_exp(q2_means,yr,yr_e)
xmin, xmax = B.pl.xlim()
B.pl.hlines([1,1.1,0.9,1.05,0.95],xmin,xmax,linestyle=['-','--','--','--','--'],
            color=['lightgray','r','r','g','g'])
B.pl.ylim(0.6,1.2)
B.pl.xlabel('$Q^2$')
B.pl.ylabel('$Y_{exp}/Y_{SIMC}$')
# twin1.plot_exp()

#%%heep yield ratio 1d histos
data_norm = 1/(D.get_eff_norm(20851)*D.get_charge_norm(20851))
radcorr = simc_norad_plots.histos_1D['W']/simc_rad_plots.histos_1D['W']

W_norm = data_norm*data_plots.histos_1D['W']*radcorr
W_stats_err = np.sqrt(data_plots.histos_1D['W'].bin_content)

simc_W = simc_norad_plots.histos_1D['W']

W_data_yield = W_norm.sum(0.92,0.96)
W_data_yield_err = np.sqrt(W_data_yield[0])

W_simc_yield = simc_W.sum(0.92,0.96)
W_simc_yield_err = np.sqrt(W_simc_yield[0])

W_ratio = W_data_yield[0]/W_simc_yield[0]

print('W_data_yield = ', float(W_data_yield[0]),'+/-',W_data_yield_err)
print('W_simc_yield = ', float(W_simc_yield[0]),'+/-',W_simc_yield_err)
print('W_ratio = ', W_ratio)

set_value(simc_rad_plots.histos_1D['Q2'].bin_content,
          value=0.,set_value_to=1.)
Q2radcorr = simc_norad_plots.histos_1D['Q2']/simc_rad_plots.histos_1D['Q2']

# print(Q2radcorr.bin_content)
Q2_norm = data_norm*data_plots.histos_1D['Q2']*Q2radcorr
Q2_stats_err = np.sqrt(data_plots.histos_1D['Q2'].bin_content)

simc_Q2 = simc_norad_plots.histos_1D['Q2']
q2_ranges = np.arange(3.25,4.25,0.25)

for r in range(q2_ranges.size-1):
    Q2_data_yield = Q2_norm.sum(q2_ranges[r],q2_ranges[r+1])
    Q2_data_yield_err = np.sqrt(Q2_data_yield[0])
    
    Q2_simc_yield = simc_Q2.sum(q2_ranges[r],q2_ranges[r+1])
    Q2_simc_yield_err = np.sqrt(Q2_simc_yield[0])
    
    Q2_ratio = Q2_data_yield[0]/Q2_simc_yield[0]
    
    print('Q2_data_yield = ', float(Q2_data_yield[0]),'+/-',Q2_data_yield_err)
    print('Q2_simc_yield = ', float(Q2_simc_yield[0]),'+/-',Q2_simc_yield_err)
    print('Q2_ratio = ', Q2_ratio)
    
Q2_data_yield = Q2_norm.sum(3.25,4.25)
Q2_data_yield_err = np.sqrt(Q2_data_yield[0])

Q2_simc_yield = simc_Q2.sum(3.25,4.25)
Q2_simc_yield_err = np.sqrt(Q2_simc_yield[0])

Q2_ratio = Q2_data_yield[0]/Q2_simc_yield[0]

print('Q2_data_yield = ', float(Q2_data_yield[0]),'+/-',Q2_data_yield_err)
print('Q2_simc_yield = ', float(Q2_simc_yield[0]),'+/-',Q2_simc_yield_err)
print('Q2_ratio = ', Q2_ratio)


#%% heep cross section
data_norm = 1/(D.get_eff_norm(20851)*D.get_charge_norm(20851))

Wradcorr = simc_norad_plots.histos_1D['W']/simc_rad_plots.histos_1D['W']
Q2radcorr = simc_norad_plots.histos_1D['Q2']/simc_rad_plots.histos_1D['Q2']

W_norm_radcorr = data_norm*data_plots.histos_1D['W']*Wradcorr
Q2_norm_radcorr = data_norm*data_plots.histos_1D['Q2']*Q2radcorr

W_xsec = W_norm_radcorr/simc_norad_plots.histos_1D['W_PS']
Q2_xsec = Q2_norm_radcorr/simc_norad_plots.histos_1D['Q2_PS']

# W_stats_err = np.sqrt(data_plots.histos_1D['W'].bin_content)/\
#                             simc_rad_plots.histos_1D['W_PS'].bin_content
W_stats_err = data_plots.histos_1D['W'].bin_error/\
                            (simc_rad_plots.histos_1D['W_PS'].bin_content*\
                                 data_plots.histos_1D['W'].bin_content)                           
# Q2_stats_err = np.sqrt(data_plots.histos_1D['Q2'].bin_content)/\
#                             simc_rad_plots.histos_1D['Q2_PS'].bin_content 
Q2_stats_err = data_plots.histos_1D['Q2'].bin_error/\
                            (simc_rad_plots.histos_1D['Q2_PS'].bin_content*\
                                data_plots.histos_1D['Q2'].bin_content)

simc_W_xsec = simc_norad_plots.histos_1D['W']/\
                                    simc_norad_plots.histos_1D['W_PS']
simc_Q2_xsec = simc_norad_plots.histos_1D['Q2']/\
                                    simc_norad_plots.histos_1D['Q2_PS']

simc_W_stats_err = np.sqrt(simc_norad_plots.histos_1D['W'].bin_content)/\
     (simc_norad_plots.histos_1D['W_PS'].bin_content*\
        simc_norad_plots.histos_1D['W'].bin_content)
simc_Q2_stats_err = np.sqrt(simc_norad_plots.histos_1D['Q2'].bin_content)/\
    (simc_norad_plots.histos_1D['Q2_PS'].bin_content*\
     simc_norad_plots.histos_1D['Q2'].bin_content)
                                    

# plot xsecs
plot_range = (W_xsec.bin_center<0.96) & (W_xsec.bin_center>0.92)

data_x = W_xsec.bin_center[plot_range]
data_y = W_xsec.bin_content[plot_range]
data_yerr = W_stats_err[plot_range]

simc_x = simc_W_xsec.bin_center[plot_range]
simc_y = simc_W_xsec.bin_content[plot_range]
simc_yerr = simc_W_stats_err[plot_range]

B.pl.figure()
B.plot_exp(data_x,data_y,data_yerr,label='data')
B.plot_exp(simc_x,simc_y,simc_yerr,label='SIMC')
B.pl.xlim(0.9,1.0)
B.pl.ylim(0.,0.002)
B.pl.xlabel('W [GeV]')
B.pl.ylabel('${d\sigma}/{d\Omega}$')
B.pl.legend()


plot_range = (Q2_xsec.bin_center<4.1) & (Q2_xsec.bin_center>3.4)

data_x = Q2_xsec.bin_center[plot_range]
data_y = Q2_xsec.bin_content[plot_range]
data_yerr = Q2_stats_err[plot_range]

simc_x = simc_Q2_xsec.bin_center[plot_range]
simc_y = simc_Q2_xsec.bin_content[plot_range]
simc_yerr = simc_Q2_stats_err[plot_range]

B.pl.figure()
B.plot_exp(data_x,data_y,data_yerr,label='data')
B.plot_exp(simc_x,simc_y,simc_yerr,label='SIMC')
B.pl.xlim(3.2,4.4)
B.pl.ylim(-0.002,0.004)
B.pl.xlabel('$Q^2$ [$GeV^2$]')
B.pl.ylabel('${d\sigma}/{d\Omega}$')
B.pl.legend()
#%%
# 2d histo

radcorr2d = simc_norad_plots.histos_2D['W_vs_Q2']/simc_rad_plots.histos_2D['W_vs_Q2']

W_vs_Q2_norm_radcorr = data_norm*data_plots.histos_2D['W_vs_Q2']*(radcorr2d)

W_vs_Q2_xsec_2d = W_vs_Q2_norm_radcorr/simc_norad_plots.histos_2D['W_vs_Q2_PS']

W_vs_Q2_xsec = W_vs_Q2_xsec_2d.project_x(range=(0.935,0.955))

W_vs_Q2_yield_2d = data_plots.histos_2D['W_vs_Q2']*radcorr2d
W_vs_Q2_yield = W_vs_Q2_yield_2d.project_x(range=(0.935,0.955))

W_vs_Q2_stats_err = np.sqrt(W_vs_Q2_yield.bin_content)/\
                        W_vs_Q2_yield.bin_content
       

W_vs_Q2_xsec_2d_simc = simc_norad_plots.histos_2D['W_vs_Q2']/\
                        simc_norad_plots.histos_2D['W_vs_Q2_PS']
                        
W_vs_Q2_xsec_simc = W_vs_Q2_xsec_2d_simc.project_x(range=(0.935,0.955))

W_vs_Q2_xsec.plot_exp()
W_vs_Q2_xsec_simc.plot_exp()

W_vs_Q2_xsec = W_vs_Q2_xsec_2d.project_y(range=(3.4,4.1))

W_vs_Q2_xsec_simc = W_vs_Q2_xsec_2d_simc.project_y(range=(3.4,4.1))

B.pl.figure()
W_vs_Q2_xsec.plot_exp() 

#%%
norm = D.get_eff_norm(20851)*D.get_charge_norm(20851)

radcorr_fac = simc_norad_plots.histos_1D['W']/simc_rad_plots.histos_1D['W']

radcorr_fac_2d = simc_norad_plots.histos_2D['W_vs_Q2']/simc_rad_plots.histos_2D['W_vs_Q2']

W_radcorr = data_plots.histos_1D['W']*radcorr_fac 

W_vs_Q2 = (1/norm)*data_plots.histos_2D['W_vs_Q2']

W_vs_Q2_sim = simc_norad_plots.histos_2D['W_vs_Q2']

brange_cut = (W_vs_Q2.y_bin_center<0.98) & (W_vs_Q2.y_bin_center>0.9)
brange = W_vs_Q2.y_bin_center[brange_cut]

counts = []
counts_sim = []
for b in brange[:]:
    h = W_vs_Q2.project_x(range=(b,b))
    hs = W_vs_Q2_sim.project_x(range=(b,b))
    
    counts.append(h.sum(3.4,4.25))
    counts_sim.append(hs.sum(3.4,4.25))
    
xbrange_cut = (W_vs_Q2.x_bin_center<4.25) & (W_vs_Q2.x_bin_center>3.4)
xbrange = W_vs_Q2.x_bin_center[brange_cut]

counts = np.array(counts)
counts_sim = np.array(counts_sim)

counts_sim_err = np.sqrt(counts_sim[:,0])
counts_err = np.sqrt(counts[:,0])

#%%
WEIGHT = D.calc_weights(SIMC_rad.Branches)

cuts_list = C.heep_event_selection_SIMC

hcoll_cut = C.coll_cut(SIMC_rad,spec='HMS',is_SIMC=True)
hcoll_cut_arrays = hcoll_cut()

zt_cut = C.ztar_cut(SIMC_rad,is_SIMC=True)
zt_cut_arrays = zt_cut()

cuts_to_apply_list = []
for cut in cuts_list:
    br = SIMC_rad.Branches[C.SIMC_names[cut.name]]
    cut_array = cut(br)
    # cut.stats()
    
    cuts_to_apply_list.append(cut_array)
# add special cut arrays to clist
cuts_to_apply_list.append(hcoll_cut_arrays)
cuts_to_apply_list.append(zt_cut_arrays)

all_cuts_sim = cuts_to_apply_list[0]
for arr in cuts_to_apply_list:    
    all_cuts_sim = all_cuts_sim & arr
    
sigcc = SIMC_rad.Branches['sigcc'][all_cuts_sim]

h_sigcc = B.histo(sigcc,range=(0,0.0025),bins=50)

