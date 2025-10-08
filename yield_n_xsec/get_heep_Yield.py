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

import numpy as np

#%% analysis functions

def apply_evt_sel_cuts(data_obj,is_SIMC=False,show_stats = False):
    if is_SIMC:
        W = data_obj.Branches['W']
        Q2 = data_obj.Branches['Q2']
        WEIGHTS = D.calc_weights(data_obj.Branches)
        WEIGHTS_PS = D.calc_weights_PS(data_obj.Branches)

        cuts_list = C.heep_event_selection_SIMC
        
        hcoll_cut = C.coll_cut(data_obj,spec='HMS',is_SIMC=True)
        hcoll_cut_arrays = hcoll_cut()
        
        zt_cut = C.ztar_cut(data_obj,is_SIMC=True)
        zt_cut_arrays = zt_cut()

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
              [c for c in cuts_list + [hcoll_cut,zt_cut]])
        
        if show_stats:
            print('Cuts STATS:\n',
                  [c.stats() for c in cuts_list])

    else:
        try:
            W = {}
            Q2 = {}

            for m in data_obj.many:
                W[m] = data_obj.Branches[m]['P.kin.primary.W']
                Q2[m] = data_obj.Branches[m]['P.kin.primary.Q2']    

            ## making the cut lists
            # first select the cuts to apply, defined already in cut_handler
            
            cuts_list = C.heep_event_selection_cuts # contains hms_delta, shms_delta,
                                                # shms_calPID,Em_cut,Q2_cut
            
            # initialize special cut classes: coll_cut, current_cut, ztar_cut,
            # CTime_cut 
            hcoll_cut = C.coll_cut(data_obj,spec='HMS',many=True)
            hcoll_cut_arrays = hcoll_cut()
            
            curr_cut = C.current_cut(data_obj,current='BCM4A',many=True)
            curr_cut_arrays = curr_cut()
            
            zt_cut = C.ztar_cut(data_obj,many=True)
            zt_cut_arrays = zt_cut()
            
            ct_cut = C.CTime_cut(data_obj,many=True)
            ct_cut_arrays = ct_cut() 
            
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
                  [c for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])
                
        except TypeError:
            
            W = data_obj.Branches['P.kin.primary.W']
            Q2 = data_obj.Branches['P.kin.primary.Q2']
    

            ## making the cut lists
            # first select the cuts to apply, defined already in cut_handler
            cuts_list = C.heep_event_selection_cuts
            
            # initialize special cut classes: coll_cut, current_cut, ztar_cut,
            # CTime_cut 
            hcoll_cut = C.coll_cut(data_obj,spec='HMS')
            hcoll_cut_arrays = hcoll_cut()
            
            curr_cut = C.current_cut(data_obj,current='BCM4A')
            curr_cut_arrays = curr_cut()
            
            zt_cut = C.ztar_cut(data_obj)
            zt_cut_arrays = zt_cut()
            
            ct_cut = C.CTime_cut(data_obj)
            ct_cut_arrays = ct_cut() 
             
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
                  [c for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])

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
        self.many = data_obj.many
        self.histos_1D = {'W':None, 'Q2':None}
        self.histos_2D = {'W_vs_Q2':None}
        
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
                if self.many is list:
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
                else:
                    raise(TypeError)
                    
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
    def make_histos(self, weights = False, nbins = 100):
        histos_to_plot = {'W':[self.W,(0.875,1.0)],'Q2':[self.Q2,(2.7,4.6)]}
        
        histos2D_to_plot = {'W_vs_Q2':[self.Q2,self.W,[(2.7,4.6),60],[(0.875,1.0),60]]}
                            

        # make 1D histograms and store them in self.histos_1D
        try:
            for v in histos_to_plot:
                h = {}
                if type(self.many) is list:
                    for m in self.many:
                        histo = histos_to_plot[v][0][m]
                        ran = histos_to_plot[v][1]
                        
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
    
    def make_PS_histos(self,nbins = 100):
        histos_to_plot = {'W_PS':[self.W,(0.875,1.0)],'Q2_PS':[self.Q2,(2.7,4.6)]}
        
        histos2D_to_plot = {'W_vs_Q2_PS':[self.Q2,self.W,[(2.7,4.6),60],[(0.875,1.0),60]]}

        # make 1D phase space histograms
        for v in histos_to_plot:
            histo = histos_to_plot[v][0]
            ran = histos_to_plot[v][1]
            
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
         'P.react.z']

br_sel_sim = ['Em','W','Q2','e_delta','h_delta','Weight','Normfac',
              'e_xptar','e_yptar','h_xptar','h_yptar','tar_x',
              'h_zv','e_zv','h_ytar','e_ytar','sig']
              
    
TSP_sel = ['evNumber','P.BCM4A.scalerCharge','P.BCM4A.scalerChargeCut',
           'P.BCM4A.scalerCurrent','P.1MHz.scalerTime']

br_sel = {'T':T_sel,'TSP':TSP_sel}
t_sel = ['T','TSP']

#%% load runs with selected branches
DATA_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

heep_dp_0 = D.DATA_INIT(data_type='deut23_data',run=20851,
                         select_branches=br_sel,
                         select_trees=t_sel,
                         ROOTfiles_path= DATA_DIR)

#%% load SIMC files
# directory paths for different SIMC files
SIMC_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/worksim" +\
                "/heep_alloffsets/"

SIMC_rad = D.DATA_INIT(data_type='SIMC',setting = 'delta_scan_0',
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='-')
SIMC_norad = D.DATA_INIT(data_type='SIMC',setting = 'delta_scan_0',
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='norad')

#%% apply event selection cuts to data
apply_evt_sel_cuts(heep_dp_0)
apply_evt_sel_cuts(SIMC_norad,is_SIMC=True)
apply_evt_sel_cuts(SIMC_rad,is_SIMC=True)

#%% create necessary plots using helper class
data_plots = yield_plots(heep_dp_0,use_cuts=True,data_type='deut23_data')
data_plots.make_histos(nbins=50)

simc_norad_plots = yield_plots(SIMC_norad,use_cuts=True,data_type='SIMC')
simc_norad_plots.make_histos(nbins=50,weights=True)

simc_rad_plots = yield_plots(SIMC_rad,use_cuts=True,data_type='SIMC')
simc_rad_plots.make_histos(nbins=50,weights=True)

#%%

radcorr_fac = simc_norad_plots.histos_1D['W']/simc_rad_plots.histos_1D['W']

W_radcorr = data_plots.histos_1D['W']*radcorr_fac 

W_vs_Q2 = data_plots.histos_2D['W_vs_Q2']

brange_cut = W_vs_Q2.x_bin_center<4.25 & W_vs_Q2.x_bin_center>3.5
brange = W_vs_Q2.x_bin_center[brange_cut]

for b in brange[:]:
    h = W_vs_Q2.project_x()
    
