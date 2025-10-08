#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:35:44 2025

@author: gvill

Extract Deuteron Cross Sections (PRELIMINARY)
"""

import data_init as D
import LT.box as B
import cut_handler as C

import numpy as np
from matplotlib import colormaps as cmp
import matplotlib.lines as mlines

#%% helper functions

def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    #B.pl.legend(fontsize = fsize)
    
def fit_legend(fit_params):
    ax = B.pl.gca()
    handles, labels = ax.get_legend_handles_labels()
    for par in fit_params:
        fit_info = (
               f"A = {par['A'].value:.2f} $\pm$ {par['A'].err:.2f}\n"
               f"$\mu$ = {par['mean'].value:.4e} $\pm$ {par['mean'].err:.4e}\n"
               f"$\sigma$ = {par['sigma'].value:.4e} $\pm$ {par['sigma'].err:.4e}\n"
           )
        handles.append(mlines.Line2D([], [], color='None', label=fit_info))
        labels.append(fit_info)
    
    #ax.legend(handles=handles, labels=labels)
    return (handles, labels)    

def set_inf(inpt,set_inf_to=0.):
    is_inf = np.where(inpt == np.inf)
    for i in is_inf[0]:
        inpt[i] = set_inf_to
    print(f'Set {is_inf[0].size} inf to {set_inf_to}\n')
    return

#%% Analysis Functions

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
            
            # print('h.bin_content[0] = ',h.bin_content[0])
            h.bin_content = enorm*h.bin_content
            # print('enorm = ',enorm)
            # print('h.bin_content[0] = ',h.bin_content[0])
            h_enorm[m] = h
            # print('h_enorm[m].bin_content[0] = ',h_enorm[m].bin_content[0])
            # print('h_enorm[m].y_bin_center = ',h_enorm[m].y_bin_center)
            
        charge_tot = np.sum(np.array(charge))
        return charge_tot,h_enorm
    else:
        print('Need to know the run # to get the normalization factors')

###
# function to combine histograms from different runs, will be normalized by
# efficiencies per run, then divided by total charge of combined runs.
###        
        
def combine_histos(histos,many=[]):
    tot_q, enorm_h = normalize_histos(histos,many)
    
    loop1 = True 
    for m in many:
        if loop1:
            hsum = enorm_h[m]
            # print('loop1: hsum.y_bin_center = ',hsum.y_bin_center)
            loop1 = False
        else:
            hsum += enorm_h[m]
            # print('hsum.y_bin_center = ',hsum.y_bin_center)
    
    # print('hsum.bin_content[0] = ',hsum.bin_content[0])
    #hsum.bin_content = (1/tot_q)*hsum.bin_content
    # print('hsum.bin_content[0] = ',hsum.bin_content[0])       
    qnorm_hsum = (1/tot_q)*hsum
    # print('qnorm_hsum.y_bin_center =', qnorm_hsum.y_bin_center)
    
    return qnorm_hsum

###
# made a function that applies the event selection cuts and stores the variables
# with cuts back into the orginal data object
###

def apply_evt_sel_cuts(data_obj,is_SIMC=False,show_stats = False):
    if is_SIMC:
        Pm = data_obj.Branches['Pm']
        Th_rq = data_obj.Branches['theta_rq']*radtodeg
        Em = data_obj.Branches['Em']
        WEIGHTS = D.calc_weights(data_obj.Branches)
        WEIGHTS_PS = D.calc_weights_PS(data_obj.Branches)

        cuts_list = C.event_selection_SIMC
        
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

        Pm_cut = Pm[all_cuts_sim]
        Th_rq_cut = Th_rq[all_cuts_sim]
        Em_cut = Em[all_cuts_sim] 
        WEIGHTS_cut = WEIGHTS[all_cuts_sim]
        WEIGHTS_PS_cut = WEIGHTS_PS[all_cuts_sim]

        data_obj.Branches.update({'Pm_cut':Pm_cut,
                                       'theta_rq_cut':Th_rq_cut,
                                       'Em_cut':Em_cut,
                                       'WEIGHTS_cut':WEIGHTS_cut,
                                       'WEIGHTS_PS_cut':WEIGHTS_PS_cut})
        print('SIMC data object updated with the branches:',
              'Pm_cut, theta_rq_cut, Em_cut\n',
              'WEIGHTS_cut, WEIGHTS_PS_cut\n')
        print('Applied the following cuts:\n',
              [c for c in cuts_list + [hcoll_cut,zt_cut]])
        
        if show_stats:
            print('Cuts STATS:\n',
                  [c.stats() for c in cuts_list])

    else:
        try:
            Pm = {}
            Th_rq = {}
            Em = {}

            for m in data_obj.many:
                Pm[m] = data_obj.Branches[m]['H.kin.secondary.pmiss']
                Th_rq[m] = data_obj.Branches[m]['H.kin.secondary.th_bq']*radtodeg
                Em[m] = data_obj.Branches[m]['H.kin.secondary.emiss_nuc']    

            ## making the cut lists
            # first select the cuts to apply, defined already in cut_handler
            
            cuts_list = C.event_selection_cuts # contains hms_delta, shms_delta,
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
                Pm_cut = Pm[m][all_cuts[m]]
                Th_rq_cut = Th_rq[m][all_cuts[m]]
                Em_cut = Em[m][all_cuts[m]]
                
                data_obj.Branches[m].update({'H.kin.secondary.pmiss_cut':Pm_cut,
                                               'H.kin.secondary.th_bq_cut':Th_rq_cut,
                                               'H.kin.secondary.emiss_nuc_cut':Em_cut})
            
            print('Data object updated with the branches:\n',
                  'H.kin.secondary.pmiss_cut, H.kin.secondary.th_bq_cut,',
                  'H.kin.secondary.emiss_nuc_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])
                
        except TypeError:
            
            Pm = data_obj.Branches['H.kin.secondary.pmiss']
            Th_rq = data_obj.Branches['H.kin.secondary.th_bq']*radtodeg
            Em = data_obj.Branches['H.kin.secondary.emiss_nuc']    

            ## making the cut lists
            # first select the cuts to apply, defined already in cut_handler
            cuts_list = C.event_selection_cuts
            
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
  
            Pm_cut = Pm[all_cuts]
            Th_rq_cut = Th_rq[all_cuts]
            Em_cut = Em[all_cuts]
            
            data_obj.Branches[m].update({'H.kin.secondary.pmiss_cut':Pm_cut,
                                           'H.kin.secondary.th_bq_cut':Th_rq_cut,
                                           'H.kin.secondary.emiss_nuc_cut':Em_cut})
        
            print('Data object updated with the branches:\n',
                  'H.kin.secondary.pmiss_cut, H.kin.secondary.th_bq_cut,',
                  'H.kin.secondary.emiss_nuc_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_list + [hcoll_cut,curr_cut,zt_cut,ct_cut]])

###
# function to get the yield based on a y-axis projection,
#   in a specified bin range in x, the projection is done bin by bin and 
#   the integral of the projected histogram is taken from y_min to y_max
# the opposite is true for projection along x-axis, set the flag project_x 
# to True to invert yield projection 
###

def get_yield(h2,bin_range,proj_min,proj_max,norm=1,plot_proj=False,project_x=False):
    
    if project_x:
        counts = []
        for ny in bin_range[:]:
            h = norm*h2.project_x(bins = [ny])
            if plot_proj:
                B.pl.figure()
                h.plot()
                B.pl.vlines([proj_min,proj_max],0,30,
                            linestyles='--',color='black')
                
            counts.append(h.sum(proj_min,proj_max))
        
        counts = np.array(counts)
           
        proj_values = h.bin_center
        yield_values = counts[:,0]
        yield_errors = np.sqrt(yield_values)/yield_values
    
    else:
        counts = []
        y_proj = []
        for nx in bin_range[:]:
            h = norm*h2.project_y(bins = [nx])
            if plot_proj:
                #B.pl.figure()
                set_inf(h.bin_content)
                h.plot_exp(ignore_zeros=True)
                B.pl.vlines([proj_min,proj_max],0,30,
                            linestyles='--',color='black')
                B.pl.yscale('log')
            
            c = h.sum(proj_min,proj_max)
            
            counts.append(c)
            y_proj.append(h)
        
        counts = np.array(counts)
        
        first = True
        for h in y_proj:
            if first:
                h_y_proj = h 
                first = False
            else:
                h_y_proj += h
                
        proj_values = h.bin_center
        yield_values = counts[:,0]
        yield_errors = np.sqrt(yield_values)
    
    return (proj_values,yield_values,yield_errors,y_proj,h_y_proj)

#%% Analysis Classes

## class to make 2d histograms for yield extraction, 
# class members are variables to plot, cuts, histograms with and without cuts,
### just to plot necessary histograms for yield extraction and save them,
# will have 1D histos of the variables as well as the needed 2D histos per run,
# not normalized and no cuts
# But SIMC histos will have weights applied, so -yes normalized

"""
CLASS yield_plots:
    
This class will create an object with the necessary histograms used in yield 
extraction, these are saved in the class member 'histos_2D' as a dictionary.
The histograms can be made with and without cuts but SIMC histograms will
always be weighted.

    data_obj: data_init class object | with loaded runs/simc data, will need to 
                have the necessary hcana/simc variables loaded for cuts and to 
                do the yield extraction
    
    use_cuts: bool | if True will use the saved '*_cut' variables assumed to 
                        be in the data_init object
    
    data_type: string | 'deut23_data' or 'SIMC'
"""

class yield_plots:
    def __init__(self, data_obj, use_cuts = False, data_type = ''):
        d = data_obj.Branches
        self.many = data_obj.many
        self.histos_1D = {'Pm':None, 'Em':None, 'th_rq':None}
        self.histos_2D = {'Pm_vs_th_rq':None, 'Em_vs_Pm':None}
        
        if data_type == 'deut23_data':
            try:
                self.Pm = {}
                self.Em = {}
                self.thrq = {}
                for m in self.many:
                    if use_cuts:
                        self.Pm[m] = d[m]['H.kin.secondary.pmiss_cut']
                        self.Em[m] = d[m]['H.kin.secondary.emiss_nuc_cut']
                        self.thrq[m] = d[m]['H.kin.secondary.th_bq_cut']
                    else:    
                        self.Pm[m] = d[m]['H.kin.secondary.pmiss']
                        self.Em[m] = d[m]['H.kin.secondary.emiss_nuc']
                        self.thrq[m] = d[m]['H.kin.secondary.th_bq']*radtodeg

            except TypeError:
                if use_cuts:
                    self.Pm = d['H.kin.secondary.pmiss_cut']
                    self.Em = d['H.kin.secondary.emiss_nuc_cut']
                    self.thrq = d['H.kin.secondary.th_bq_cut']
                else:
                    self.Pm = d['H.kin.secondary.pmiss']
                    self.Em = d['H.kin.secondary.emiss_nuc']
                    self.thrq = d['H.kin.secondary.th_bq']*radtodeg
        
        elif data_type == 'SIMC':
            try:
                self.Pm = {}
                self.Em = {}
                self.thrq = {}
                self.WEIGHTS = {}
                self.WEIGHTS_PS = {}
                if self.many is list:
                    for m in self.many:
                        if use_cuts:
                            self.Pm[m] = d[m]['Pm_cut']
                            self.Em[m] = d[m]['Em_cut']
                            self.thrq[m] = d[m]['theta_rq_cut']
                            self.WEIGHTS[m] = d[m]['WEIGHTS_cut']
                            self.WEIGHTS_PS[m] = d[m]['WEIGHTS_PS_cut']
                            
                        else:    
                            self.Pm[m] = d[m]['Pm']
                            self.Em[m] = d[m]['Em']
                            self.thrq[m] = d[m]['theta_rq']
                            self.WEIGHTS[m] = D.calc_weights(d[m])
                            self.WEIGHTS_PS[m] = D.calc_weights_PS(d[m])
                else:
                    raise(TypeError)
                    
            except TypeError:               
                if use_cuts:
                    self.Pm = d['Pm_cut']
                    self.Em = d['Em_cut']
                    self.thrq = d['theta_rq_cut']
                    self.WEIGHTS = d['WEIGHTS_cut']
                    self.WEIGHTS_PS = d['WEIGHTS_PS_cut']
                    
                else:    
                    self.Pm = d['Pm']
                    self.Em = d['Em']
                    self.thrq = d['theta_rq']
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
        histos_to_plot = {'Pm':[self.Pm,(-0.1,1.5)], 'Em':[self.Em,(-0.1,0.1)],
            'th_rq':[self.thrq,(0,180)]}
        
        histos2D_to_plot = {'Pm_vs_th_rq':[self.thrq,self.Pm,[(0,180),60],[(-0.1,1.5),38]],
                            'Em_vs_Pm':[self.Pm,self.Em,[(-0.1,1.5),38],[(-0.1,0.1),38]]}

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
        histos_to_plot = {'Pm_PS':[self.Pm,(-0.1,1.5)], 'Em_PS':[self.Em,(-0.1,0.1)],
            'th_rq_PS':[self.thrq,(0,180)]}
        
        histos2D_to_plot = {'Pm_vs_th_rq_PS':[self.thrq,self.Pm,[(0,180),60],[(-0.1,1.5),38]],
                            'Em_vs_Pm_PS':[self.Pm,self.Em,[(-0.1,1.5),38],[(-0.1,0.1),38]]}

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


#%% conversion factors
radtodeg = 180./np.pi

#%% select appropriate HCANA variables 
kin_var = ['H.kin.secondary.th_bq','H.kin.secondary.pmiss',
           'H.kin.secondary.emiss_nuc','P.kin.primary.Q2']

sieve_var = ['H.extcor.xsieve', 'H.extcor.ysieve',
                  'P.extcor.xsieve', 'P.extcor.ysieve']

coinTime_var = ['CTime.epCoinTime_ROC2']

acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph']

calPID_var = ['P.cal.etottracknorm']

z_react_var = ['H.react.z','P.react.z']

T_sel = kin_var + coinTime_var + acc_var + calPID_var + sieve_var + z_react_var

br_sel_sim = ['theta_rq','Pm','Em','Q2','e_delta','h_delta','Weight','Normfac',
              'e_xptar','e_yptar','h_xptar','h_yptar','tar_x','h_zv','e_zv',
              'h_ytar','e_ytar']
    
TSP_sel = ['evNumber','P.BCM4A.scalerCharge','P.BCM4A.scalerChargeCut',
           'P.BCM4A.scalerCurrent',
           'P.BCM4B.scalerCharge','P.BCM4B.scalerCurrent',
           'P.BCM4C.scalerCharge','P.BCM4C.scalerCurrent',
           'P.1MHz.scalerTime']

br_sel = {'T':T_sel,'TSP':TSP_sel}
t_sel = ['T','TSP']


#%% load runs with selected branches
run_list = [20871,20872,20871,20872,21065,21066,21067,21068,21069,21070,
            21071,21072,21073,21074,21075]

DATA_baseDIR = "/media/gvill/Gema's T7/ROOTfiles"
pm_120 = '/pass_2/'

deep_pm120 = D.DATA_INIT(data_type='deut23_data',run=run_list,
                         select_branches=br_sel,
                         select_trees=t_sel,
                         ROOTfiles_path= DATA_baseDIR+pm_120)

#%% load SIMC files
# directory paths for different SIMC files
SIMC_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/worksim" +\
                "/deep_pm120/"

SIMC_FSIrad = D.DATA_INIT(data_type='SIMC',setting = 'pm_120',
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='jmlfsi_rad')
SIMC_FSInorad = D.DATA_INIT(data_type='SIMC',setting = 'pm_120',
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='jmlfsi_norad')
SIMC_PWIAnorad = D.DATA_INIT(data_type='SIMC',setting = 'pm_120',
                              select_branches={'SNT':br_sel_sim},
                              SIMC_ROOTfiles_path=SIMC_baseDIR,
                              simc_type='jmlpwia_norad')

               
#%% Prepare the data and make histograms
# in this part of the code I apply cuts to the needed DATA ans SIMC variables and
# save them back in the DATA_INIT object as '[variable_name]_cut'

apply_evt_sel_cuts(deep_pm120)
apply_evt_sel_cuts(SIMC_FSIrad, is_SIMC=True)
apply_evt_sel_cuts(SIMC_FSInorad, is_SIMC=True)
apply_evt_sel_cuts(SIMC_PWIAnorad, is_SIMC=True)
    
#%% use the class yield_plots to create the data histograms
dataPlots_cut = yield_plots(deep_pm120,use_cuts=True,data_type='deut23_data')
dataPlots_cut.make_histos()

# use combine_histos() function to normalize and combine histograms from 
# different runs, the efficiency normalization is done per run by another 
# function normalize_histos() and calculates the total charge, then the 
# charge normalization is done at the end and that histo is returned.
 
pm120_h = combine_histos(dataPlots_cut.histos_2D['Pm_vs_th_rq'],deep_pm120.many)

#%% create necessary histograms from SIMC

simcPlots_cut_FSIrad = yield_plots(SIMC_FSIrad,use_cuts=True,data_type='SIMC')
simcPlots_cut_FSIrad.make_histos(weights=True)

simcPlots_cut_FSInorad = yield_plots(SIMC_FSInorad,use_cuts=True,data_type='SIMC')
simcPlots_cut_FSInorad.make_histos(weights=True)
simcPlots_cut_FSInorad.make_PS_histos()

simcPlots_cut_PWIAnorad = yield_plots(SIMC_PWIAnorad,use_cuts=True,data_type='SIMC')
simcPlots_cut_PWIAnorad.make_histos(weights=True)
simcPlots_cut_PWIAnorad.make_PS_histos()

FSI_ratio_h = simcPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                simcPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq']
# set nans from division to zero
FSI_ratio_h.set_nans()                                

#%% histogram manipulation to get xsec

pm120_h_radcorr = FSI_ratio_h*pm120_h

pm120_h_Xsec = pm120_h_radcorr/\
                simcPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq_PS']
# set nans from division to zero
pm120_h_Xsec.set_nans()
set_inf(pm120_h_Xsec.bin_content)

SIMC_pm120_h_FSIXsec = simcPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                         simcPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq_PS']
SIMC_pm120_h_FSIXsec.set_nans()
set_inf(SIMC_pm120_h_FSIXsec.bin_content) 
                        
SIMC_pm120_h_PWIAXsec = simcPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                         simcPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq_PS']
SIMC_pm120_h_PWIAXsec.set_nans() 
set_inf(SIMC_pm120_h_PWIAXsec.bin_content)                        

#%% set up histogram projections

angles = np.arange(0,10,9)
# angles = [45]
rows = 4
cols = 4
i=1
# B.pl.figure(figsize=(20,10),layout='constrained')
B.pl.figure()
for a in angles:
    
    bins_x_min = np.where(pm120_h_Xsec.x_bin_center >= a-9)[0][0]
    bins_x_max = np.where(pm120_h_Xsec.x_bin_center <= a+9)[0][-1]
    
    if bins_x_min == bins_x_max:
        brange = [int(bins_x_max)]
    else:    
        brange = np.arange(int(bins_x_min),int(bins_x_max))
    # brange = np.arange(0,pm120_h_Xsec.nbins_y)
    
    pm_values, pmth_xsec, pmth_error, proj_list, proj =\
        get_yield(h2=pm120_h_Xsec, 
                  bin_range=brange, proj_min=0., proj_max=0.5)
        
    pm, pmth_yield, pmth_error, prl, pr0 =\
        get_yield(h2=pm120_h_radcorr, bin_range=brange, proj_min=0., proj_max=0.5)
    
    # ax = B.pl.subplot(rows,cols,i)
    B.pl.figure()
    # if i <= 16:
    #     i+=1
    

    B.plot_exp(plot_proj_x,plot_proj_y,plot_proj_dy,logy=True)
    # proj.plot_exp(ignore_zeros=True,logy=True)
    
    pm_values_sim, pmth_xsec_sim, pmth_error_sim, l2, pr2 =\
        get_yield(h2=SIMC_pm120_h_PWIAXsec, 
                  bin_range=brange, proj_min=0., proj_max=0.5)    
    
    # pr2.plot_exp(ignore_zeros=True, color='r')
    
    pm_values_sim, pmth_xsec_sim, pmth_error_sim, l3, pr3 =\
        get_yield(h2=SIMC_pm120_h_FSIXsec, 
                  bin_range=brange, proj_min=0., proj_max=0.5)    
    # pr3.plot_exp(ignore_zeros=True, color='b')
    B.pl.xlim((-0.01,1.5))
    B.pl.ylim((0,1))
    # B.pl.yscale('log')
    B.pl.title(f'{a}$\pm$9')
    B.pl.xlabel('')
    B.pl.ylabel('')
    
# B.plot_exp(pm_values_sim, pmth_xsec_sim)
# B.pl.yscale('log')

#%% plotting xsec

fig_to_plot = {'Data':[pm_values, pmth_xsec, pmth_error,'#1f77b4'],
'SIMC':[pm_values_sim, pmth_xsec_sim, pmth_error_sim,'#de425b']}
B.pl.figure()
for fig in fig_to_plot:   
    B.plot_exp(x = fig_to_plot[fig][0],
               y = fig_to_plot[fig][1],
               
               c=fig_to_plot[fig][3], marker='+', markersize=8,
                     capsize=0,mew=1,elinewidth=1,label=fig)