#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:35:44 2025

@author: gvill

Extract Deuteron Cross Sections (PRELIMINARY)

This code utilizes modules to load data, apply cuts, then produce histograms
and manipulates them to get cross sections, then saves all relevant plots to
a pickle file

"""

import data_init as D
import LT.box as B
import cut_handler as C
import load_data_init as L
import database_operations as db

import pickle as P
import numpy as np

rtd = 180/np.pi

#%% helper functions
def set_inf(inpt,set_inf_to=0.):
    is_inf = np.where(inpt == np.inf)
    for i in is_inf[0]:
        inpt[i] = set_inf_to
    print(f'Set {is_inf[0].size} inf to {set_inf_to}\n')
    return

def save_histos(d,h):
    hist = h.histos_2D
    filename = f'./yield_n_xsec/average_kin/{d.setting}_{d.simc_type}_2Dhistos_avgKin.pcl'
    with open(filename, 'wb') as f:
        P.dump(hist,f)
    print(f'Saved {filename}\n')
#%% Analysis Functions

###
# function to combine histograms from different runs, will be normalized by
# efficiencies per run, then divided by total charge of combined runs.
# updated to deal with different data sets
###        
        
def comb_n_norm_histos(histos,many=[]):   
    combined_histos = {}
    for i in many:
        char = []
        enorm_h = 0
        for m in many[i]:
            enorm = D.get_eff_norm(m)
            char.append(D.get_charge_norm(m,curr='BCM1'))
            
            enorm_h += histos[m]/enorm
        
        tot_char = np.sum(np.array(char))
        combined_histos[i] = enorm_h/tot_char
    
    return combined_histos

###
# function to combine histograms but does not normalize them.
# updated to handle sets of data
###

def combine_histos(histos,many=[]): 
    combined_histos = {}
    for i in many:
        hsum = 0
        for m in many[i]:
            hsum += histos[m]
        combined_histos[i] = hsum
    return combined_histos



###
# function to project 2D histograms along x or y axis for a given bin range,
# there are options to also get the integral of each projection by setting 
# 'get_counts = True', and the range of integration can also be set in 
# 'int_range'
#       h2 -> 2D histogram to project
#       project_along -> string | options: 'x', 'y' | axis of projection
#       proj_bin_range -> tuple | set bin range on the projected axis, otherwise
#                                    project all bins
#       get_counts ------> bool | sum all bin contents of each projection, return 
#                                    list of counts as well
#       int_range ------> tuple | set range in axis of projection to get counts over
#       proj_title ----> string | set title name of each projection histogram
###

def project_2D(h2,project_along='y',proj_bin_range=(),get_counts=False,int_range=[],
               proj_title=''):
    if project_along == 'y':
        if proj_bin_range:
            x_nbins_min = proj_bin_range[0]
            x_nbins_max = proj_bin_range[1]
        else:
            x_nbins_min = 0
            x_nbins_max = h2.nbins_x
        xbin_width = h2.x_bin_width
        
        y_proj_h = []
        int_counts = []
        # proj_stats_err = []
        for i in range(x_nbins_min,x_nbins_max):
            xbin_center = h2.x_bin_center[i]
            hproj = h2.project_y(bins=[i])            
            hproj.title = f'{proj_title} = {xbin_center}$\pm${xbin_width/2}'
            hproj.xlabel = ''
            hproj.ylabel = ''
            
            #calculate stats error and store
            # if stats_err:
            #     c = hproj.bin_content
            #     err = np.sqrt(c)
            #     proj_stats_err.append(err)                
            
            
            if get_counts:
                if int_range:
                    bin_counts = hproj.sum(int_range)
                else:
                    bin_counts = hproj.sum()
                
                int_counts.append(bin_counts)
                hproj.title = hproj.title + f' (N = {bin_counts[0]:.0f})'
            
            y_proj_h.append(hproj)
        if get_counts:
            return (y_proj_h,np.array(int_counts))
        else:
            return y_proj_h
        
    elif project_along == 'x':
        if proj_bin_range:
            y_nbins_min = proj_bin_range[0]
            y_nbins_max = proj_bin_range[1]
        else:
            y_nbins_min = 0
            y_nbins_max = h2.nbins_y
        ybin_width = h2.y_bin_width
        
        x_proj_h = []
        int_counts = []
        for i in range(y_nbins_min,y_nbins_max):
            ybin_center = h2.y_bin_center[i]
            hproj = h2.project_x(bins=[i])            
            hproj.title = f'{proj_title} = {ybin_center}$\pm${ybin_width/2}'
            hproj.xlabel = ''
            hproj.ylabel = ''
            
            if get_counts:
                if int_range:
                    bin_counts = hproj.sum(int_range)
                else:
                    bin_counts = hproj.sum()
                
                int_counts.append(bin_counts)
                hproj.title = hproj.title + f' (N = {bin_counts[0]:.0f})'
            
            x_proj_h.append(hproj)
        if get_counts:
            return (x_proj_h,np.array(int_counts))
        else:
            return x_proj_h

###    
# this function cuts out data points outside the phase space
# must give phase space histograms in same bins
# function fits phase space histos and cuts out data points
# outside 1.5 sigma of fit
##

def inside_PS(PS_histo,histo,sigma_fac=1.0):
    corr_histos = []
    for i in range(len(histo)):
        print(f'histo {i}')
        if any(PS_histo[i].bin_content > 100.):
            PS_histo[i].fit(ignore_zeros=True,plot_fit=False)
            j = 2.
            while PS_histo[i].chi2_red > 10. and j > 0.5:
                PS_histo[i].fit(PS_histo[i].mean.value - j*PS_histo[i].sigma.value,
                                PS_histo[i].mean.value + j*PS_histo[i].sigma.value,
                                ignore_zeros=True,plot_fit=False)
                j = j - 0.5
            xmin = PS_histo[i].mean.value - sigma_fac*PS_histo[i].sigma.value
            xmax = PS_histo[i].mean.value + sigma_fac*PS_histo[i].sigma.value
        else:
            xmin = -1e6
            xmax = 1e6
        
        is_inside_PS = (histo[i].bin_center > xmin) & (histo[i].bin_center < xmax)
        
        newhisto = B.histo(bin_center=histo[i].bin_center[is_inside_PS],
                           bin_content=histo[i].bin_content[is_inside_PS],
                           bin_error=histo[i].bin_error[is_inside_PS],
                           title=histo[i].title,
                           xlabel=histo[i].xlabel,
                           ylabel=histo[i].ylabel)
        corr_histos.append(newhisto)
    return corr_histos
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
    def __init__(self, data_obj, use_cuts = False, data_type = '',
                 with_weights=True,calc_avg_kin=False):
        d = data_obj.Branches
        self.many = data_obj.many
        self.histos_1D = {'Pm':None, 'Em':None, 'th_rq':None}
        self.histos_2D = {'Pm_vs_th_rq':None, 'Em_vs_Pm':None}
        self.has_weights = with_weights
        
        
        #set histogram limits to default values,
        # set_hist_lim can be called outside after initializing
        # yield_plots to set custom ranges or bins, otherwise
        # default values will be used
        self.set_hist_lim()
        
        if data_type == 'deut23_data':
            self.has_weights = False
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
                        self.thrq[m] = d[m]['H.kin.secondary.th_bq']*rtd

            except TypeError:
                if use_cuts:
                    self.Pm = d['H.kin.secondary.pmiss_cut']
                    self.Em = d['H.kin.secondary.emiss_nuc_cut']
                    self.thrq = d['H.kin.secondary.th_bq_cut']
                else:
                    self.Pm = d['H.kin.secondary.pmiss']
                    self.Em = d['H.kin.secondary.emiss_nuc']
                    self.thrq = d['H.kin.secondary.th_bq']*rtd
        
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
                    
                    if calc_avg_kin:
                        self.Ein_v = d['Ein_v_cut']
                        self.kf_v = d['kf_v_cut']
                        self.the_v = d['the_v_cut']
                        self.pf_v = d['pf_v_cut']
                        self.thp_v = d['thp_v_cut']
                        self.q_v = d['q_v_cut']
                        self.q_v_calc = d['q_v_calc_cut']
                        self.thq_v = d['thq_v_cut']
                        self.Q2_v = d['Q2_v_cut']
                        self.omega_v = d['omega_v_cut']
                        self.xbj_v = d['xbj_v_cut']
                        self.pm_v = d['pm_v_cut']
                        self.pm_v_calc = d['pm_v_calc_cut']
                        self.th_pq_v = d['th_pq_v_cut']
                        self.th_nq_v = d['th_nq_v_cut']
                        self.ph_pq_v = d['ph_pq_v_cut']
                        self.ph_nq_v = d['ph_nq_v_cut']
                        self.sig = d['sig_cut']
                    
                else:    
                    self.Pm = d['Pm']
                    self.Em = d['Em']
                    self.thrq = d['theta_rq']
                    self.WEIGHTS = D.calc_weights(d)
                    self.WEIGHTS_PS = D.calc_weights_PS(d)
                    
                    if calc_avg_kin:
                        self.Ein_v = d['Ein_v_GeV']
                        self.kf_v = d['kf_v']
                        self.the_v = d['the_v']
                        self.pf_v = d['pf_v_GeV']
                        self.thp_v = d['thp_v']
                        self.q_v = d['q_lab_v_GeV']
                        self.q_v_calc = d['q_v_calc']
                        self.thq_v = d['thq_v']
                        self.Q2_v = d['Q2_v_GeV2']
                        self.omega_v = d['nu_v_GeV']
                        self.xbj_v = d['xbj_v']
                        self.pm_v = d['pm_v_GeV']
                        self.pm_v_calc = d['pm_v_calc']
                        self.th_pq_v = d['th_pq_v']
                        self.th_nq_v = d['th_nq_v']
                        self.ph_pq_v = d['ph_pq_v']
                        self.ph_nq_v = d['ph_nq_v']
        else:
            print('No data type chosen.\n',
                  'Available types: "deut23_data", "SIMC"')
    ###
    # method to create desired histograms, note for SIMC histograms will be 
    # weighted if the flag weights is set to True.
    ###    
    def make_histos(self):
        histos_to_plot = {'Pm':[self.Pm,self.pm_range,self.pm_bins], 
                          'Em':[self.Em,self.Em_range,self.Em_bins],
                          'th_rq':[self.thrq,self.thrq_range,self.thrq_bins]}
        
        histos2D_to_plot = {'Pm_vs_th_rq':[self.thrq,self.Pm,
                                           [self.thrq_range,self.thrq_bins],
                                           [self.pm_range,self.pm_bins]],
                            'Em_vs_Pm':[self.Pm,self.Em,
                                        [self.pm_range,self.pm_bins],
                                        [self.Em_range,self.Em_bins]]}

        # make 1D histograms and store them in self.histos_1D
        try:
            for v in histos_to_plot:
                h = {}
                if type(self.many) is list:
                    for m in self.many:
                        histo = histos_to_plot[v][0][m]
                        ran = histos_to_plot[v][1]
                        nbins = histos_to_plot[v][2]
                        
                        if self.has_weights:
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
                
                if self.has_weights:
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
                        
                        if self.has_weights:
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
                
                if self.has_weights:
                    h = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins],
                                     weights = self.WEIGHTS,
                                     calc_w2=True)
                else:
                    h = B.histo2d(hx,hy,range=[ranx,rany],bins=[nxbins,nybins])
                
                self.histos_2D[v] = h
    
    def make_PS_histos(self):
        W = self.WEIGHTS_PS
        histos_to_plot = {'Pm_PS':[self.Pm,W,self.pm_range,self.pm_bins], 
                          'Em_PS':[self.Em,W,self.Em_range,self.Em_bins],
                          'th_rq_PS':[self.thrq,W,self.thrq_range,self.thrq_bins]}
        
        histos2D_to_plot = {'Pm_vs_th_rq_PS':[self.thrq,self.Pm,W,
                                           [self.thrq_range,self.pm_range],
                                           [self.thrq_bins,self.pm_bins]],
                            'Em_vs_Pm_PS':[self.Pm,self.Em,W,
                                        [self.pm_range,self.Em_range],
                                        [self.pm_bins,self.Em_bins]]}
        # make 1D phase space histograms
        self.__make_histos1D(histos_to_plot)

        # make 2D phase space histograms
        self.__make_histos2D(histos2D_to_plot)
    
    ###
    # method to create avg. kinematics histograms
    ###    
    def make_avgKin_histos(self):
        pm_v = self.pm_v
        thrq_v = self.th_nq_v
        
        histos_to_plot = {'pm_v':[pm_v,self.WEIGHTS,(-0.1,1.5),38],
                          'pm_v_calc':[self.pm_v_calc,self.WEIGHTS,(-0.1,1.5),38],
                          'th_nq_v':[thrq_v,self.WEIGHTS,(0,190),19],
                          'Ein_v':[self.Ein_v,self.WEIGHTS,(10.5,10.6),100],
                          'kf_v':[self.kf_v,self.WEIGHTS,(7.75,9),100],
                          'the_v':[self.the_v,self.WEIGHTS,(10,14),40],
                          'pf_v':[self.pf_v,self.WEIGHTS,(2.5,3.5),100],
                          'thp_v':[self.thp_v,self.WEIGHTS,(0,180),100],
                          'q_v':[self.q_v,self.WEIGHTS,(2,3.5),100],
                          'q_v_calc':[self.q_v_calc,self.WEIGHTS,(2,3.5),100],
                          'thq_v':[self.thq_v,self.WEIGHTS,(30,55),25],
                          'Q2_v':[self.Q2_v,self.WEIGHTS,(2.5,5.5),60],
                          'omega_v':[self.omega_v,self.WEIGHTS,(1.5,3),100],
                          'xbj_v':[self.xbj_v,self.WEIGHTS,(0.8,2.0),100],
                          'th_pq_v':[self.th_pq_v,self.WEIGHTS,(0,180),100],
                          'ph_pq_v':[self.ph_pq_v,self.WEIGHTS,(-180,180),100],
                          'ph_nq_v':[self.ph_nq_v,self.WEIGHTS,(-180,180),100],
                          'cph_pq_v':[np.cos(self.ph_pq_v),self.WEIGHTS,(-1,1),100],
                          'sph_pq_v':[np.sin(self.ph_pq_v),self.WEIGHTS,(-1,1),100]}
        
        W = self.WEIGHTS
        Ein_W = self.Ein_v*self.WEIGHTS
        kf_W = self.kf_v*self.WEIGHTS
        the_W = self.the_v*self.WEIGHTS
        pf_W = self.pf_v*self.WEIGHTS
        thp_W = self.thp_v*self.WEIGHTS
        q_W = self.q_v*self.WEIGHTS
        thq_W = self.thq_v*self.WEIGHTS
        Q2_W = self.Q2_v*self.WEIGHTS
        om_W = self.omega_v*self.WEIGHTS
        xbj_W = self.xbj_v*self.WEIGHTS
        pm_W = pm_v*self.WEIGHTS
        thpq_W = self.th_pq_v*self.WEIGHTS
        thnq_W = self.th_nq_v*self.WEIGHTS
        phpq_W = self.ph_pq_v*self.WEIGHTS
        phnq_W = self.ph_nq_v*self.WEIGHTS
        cphpq_W = np.cos(self.ph_pq_v/rtd)*self.WEIGHTS
        sphpq_W = np.sin(self.ph_pq_v/rtd)*self.WEIGHTS
        sig_W = self.sig*W

        histos2D_to_plot = {'Pm_vs_thnq_v':[thrq_v,pm_v,W,
                                            [self.thrq_range,self.pm_range],
                                            [self.thrq_bins,self.pm_bins]],
                            'Ein_2Davg':[thrq_v,pm_v,Ein_W,
                                        [self.thrq_range,self.pm_range],
                                        [self.thrq_bins,self.pm_bins]],
                            'kf_2Davg':[thrq_v,pm_v,kf_W,
                                       [self.thrq_range,self.pm_range],
                                       [self.thrq_bins,self.pm_bins]],
                            'the_2Davg':[thrq_v,pm_v,the_W,
                                        [self.thrq_range,self.pm_range],
                                        [self.thrq_bins,self.pm_bins]],
                            'pf_2Davg':[thrq_v,pm_v,pf_W,
                                       [self.thrq_range,self.pm_range],
                                       [self.thrq_bins,self.pm_bins]],
                            'thp_2Davg':[thrq_v,pm_v,thp_W,
                                        [self.thrq_range,self.pm_range],
                                        [self.thrq_bins,self.pm_bins]],
                            'q_2Davg':[thrq_v,pm_v,q_W,
                                      [self.thrq_range,self.pm_range],
                                      [self.thrq_bins,self.pm_bins]],
                            'thq_2Davg':[thrq_v,pm_v,thq_W,
                                        [self.thrq_range,self.pm_range],
                                        [self.thrq_bins,self.pm_bins]],
                            'Q2_2Davg':[thrq_v,pm_v,Q2_W,
                                       [self.thrq_range,self.pm_range],
                                       [self.thrq_bins,self.pm_bins]],
                            'omega_2Davg':[thrq_v,pm_v,om_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'xbj_2Davg':[thrq_v,pm_v,xbj_W,
                                        [self.thrq_range,self.pm_range],
                                        [self.thrq_bins,self.pm_bins]],
                            'Pm_2Davg':[thrq_v,pm_v,pm_W,
                                       [self.thrq_range,self.pm_range],
                                       [self.thrq_bins,self.pm_bins]],
                            'th_pq_2Davg':[thrq_v,pm_v,thpq_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'th_nq_2Davg':[thrq_v,pm_v,thnq_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'ph_pq_2Davg':[thrq_v,pm_v,phpq_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'ph_nq_2Davg':[thrq_v,pm_v,phnq_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'cph_pq_2Davg':[thrq_v,pm_v,cphpq_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'sph_pq_2Davg':[thrq_v,pm_v,sphpq_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]],
                            'sig_2Davg':[thrq_v,pm_v,sig_W,
                                          [self.thrq_range,self.pm_range],
                                          [self.thrq_bins,self.pm_bins]]}  
        
        # make 1D histograms and store them in self.histos_1D
        self.__make_histos1D(histos_to_plot)
        
        # make 2D histograms and store them in self.histos_2D
        self.__make_histos2D(histos2D_to_plot)
        
        # finish avg histos by dividing by Pm_vs_thnq_v
        self.__divide_avg_histos(self.histos_2D)
    
    def set_hist_lim(self,pm_range=(-0.1,1.5),pm_bins=38,thrq_range=(0,125),
                     thrq_bins=25,Em_range=(-0.1,0.1),Em_bins=38):
        self.pm_range = pm_range
        self.pm_bins = pm_bins
        self.thrq_range = thrq_range
        self.thrq_bins = thrq_bins
        self.Em_range = Em_range
        self.Em_bins = Em_bins
        
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

    def __make_histos1D(self,histo_dict):
        for v in histo_dict:
            histo = histo_dict[v][0]
            w = histo_dict[v][1]
            ran = histo_dict[v][2]
            nbins = histo_dict[v][3]
            
            if self.has_weights:
                h = B.histo(histo,range = ran, bins = nbins,
                            weights = w, calc_w2=True)
            else:    
                h = B.histo(histo,range = ran, bins = nbins)
                
            self.histos_1D[v] = h
    
    def __make_histos2D(self,histo_dict):
        for v in histo_dict:
            hx = histo_dict[v][0]
            hy = histo_dict[v][1]
            w = histo_dict[v][2]
            ran = histo_dict[v][3]
            nbins = histo_dict[v][4]
            
            if self.has_weights:
                h = B.histo2d(hx, hy, range = ran, bins = nbins,
                                 weights = w, calc_w2=True)
            else:
                h = B.histo2d(hx,hy,range=ran,bins=nbins)
            
            self.histos_2D[v] = h
            
    def __divide_avg_histos(self,histo_dict):
        for histo in histo_dict:
            if histo == 'Pm_vs_thnq_v' or histo[-3:] != 'avg':
                continue
            else:
                h = histo_dict[histo]/histo_dict['Pm_vs_thnq_v']
                
                h.title = histo
                h.xlabel = '$\\theta_{nq}$ (vertex)'
                h.ylabel = '$P_m$ (vertex)'
                
                histo_dict[histo] = h
        
#%% load runs with selected branches
DATA_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/pass_3/"
SIMC_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/deep_alloffsets/"

load_simc = ['JML']
# load_simc = ['Paris','CD-Bonn']
# load_simc = ['V18']

load_data = input('Select a setting to load:\npm_120, pm_580, pm_800, pm_900\n')

F = L.load_data_init(setting=load_data,data_dir=DATA_baseDIR,
                   simc_dir=SIMC_baseDIR,simc_models=load_simc,
                   load_comm_data=False,load_data=False)

calc_vertex_kinematics = False

#%% calculate avg kinematics variables
print('** Going to calculate avg kin...')
calc_vertex_kinematics = True

if calc_vertex_kinematics:
    if F.load_JML:
        L.calc_vertex_kin(F.JML_FSInorad)
        L.calc_vertex_kin(F.JML_PWIAnorad)
    if F.load_PAR:
        L.calc_vertex_kin(F.PAR_FSInorad)
        L.calc_vertex_kin(F.PAR_PWIAnorad)
            
#%% Apply cuts to data and SIMC
# in this part of the code I apply cuts to the needed DATA and SIMC variables and
# save them back in the DATA_INIT object as '[variable_name]_cut'

print('** Applying cuts to data...')
if F.load_data:
    C.apply_evt_sel_cuts(F.data,PCOLL=False)
if F.load_JML:
    C.apply_evt_sel_cuts(F.JML_FSIrad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.JML_PWIArad, is_SIMC=True,PCOLL=False)
    if calc_vertex_kinematics:
        C.apply_evt_sel_cuts(F.JML_FSInorad, is_SIMC=True,PCOLL=False,calc_avg_kin=True)
        C.apply_evt_sel_cuts(F.JML_PWIAnorad, is_SIMC=True,PCOLL=False,calc_avg_kin=True)
    else:
        C.apply_evt_sel_cuts(F.JML_FSInorad, is_SIMC=True,PCOLL=False)
        C.apply_evt_sel_cuts(F.JML_PWIAnorad, is_SIMC=True,PCOLL=False)       
if F.load_PAR:
    C.apply_evt_sel_cuts(F.PAR_FSIrad, is_SIMC=True,PCOLL=False)
    if calc_vertex_kinematics:
        C.apply_evt_sel_cuts(F.PAR_FSInorad, is_SIMC=True,PCOLL=False,calc_avg_kin=True)
        C.apply_evt_sel_cuts(F.PAR_PWIAnorad, is_SIMC=True,PCOLL=False,calc_avg_kin=True)
    else:
        C.apply_evt_sel_cuts(F.PAR_FSInorad, is_SIMC=True,PCOLL=False)
        C.apply_evt_sel_cuts(F.PAR_PWIAnorad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.PAR_PWIArad, is_SIMC=True,PCOLL=False)
if F.load_V18:
    C.apply_evt_sel_cuts(F.V18_FSIrad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.V18_FSInorad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.V18_PWIAnorad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.V18_PWIArad, is_SIMC=True,PCOLL=False)
if F.load_CDB:
    C.apply_evt_sel_cuts(F.CDB_FSIrad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.CDB_FSInorad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.CDB_PWIAnorad, is_SIMC=True,PCOLL=False)
    C.apply_evt_sel_cuts(F.CDB_PWIArad, is_SIMC=True,PCOLL=False)
    
#%% use the class yield_plots to create the data histograms
print('** Creating data plots...')
if F.load_data:
    dataPlots_cut = yield_plots(F.data,use_cuts=True,data_type='deut23_data')
    dataPlots_cut.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                               pm_range=(3.469446951953614e-18,1.2),
                               pm_bins=30)
    dataPlots_cut.make_histos()
    
    dataPlots_nocut = yield_plots(F.data,use_cuts=False,
                                  data_type='deut23_data')
    dataPlots_nocut.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                               pm_range=(3.469446951953614e-18,1.2),
                               pm_bins=30)
    dataPlots_nocut.make_histos()

# create necessary histograms from SIMC
print('** Creating SIMC plots...')
# JML
if F.load_JML:
    jmlPlots_cut_FSIrad = yield_plots(F.JML_FSIrad,use_cuts=True,
                                       data_type='SIMC')
    jmlPlots_cut_FSIrad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                     pm_range=(3.469446951953614e-18,1.2),
                                     pm_bins=30)
    jmlPlots_cut_FSIrad.make_histos()
    
    if calc_vertex_kinematics:
        jmlPlots_cut_FSInorad = yield_plots(F.JML_FSInorad,use_cuts=True,
                                             data_type='SIMC',calc_avg_kin=True)
    else:
        jmlPlots_cut_FSInorad = yield_plots(F.JML_FSInorad,use_cuts=True,
                                             data_type='SIMC')
        
    # jmlPlots_cut_FSInorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
    #                                    pm_range=(3.469446951953614e-18,1.2),
    #                                    pm_bins=30)
    jmlPlots_cut_FSInorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                       pm_range=(3.469446951953614e-18,1.2),
                                       pm_bins=30)
    jmlPlots_cut_FSInorad.make_histos()
    jmlPlots_cut_FSInorad.make_PS_histos()
    if calc_vertex_kinematics:
        jmlPlots_cut_FSInorad.make_avgKin_histos()
    
    if calc_vertex_kinematics:
        jmlPlots_cut_PWIAnorad = yield_plots(F.JML_PWIAnorad,use_cuts=True,
                                              data_type='SIMC',calc_avg_kin=True)
    else:
        jmlPlots_cut_PWIAnorad = yield_plots(F.JML_PWIAnorad,use_cuts=True,
                                              data_type='SIMC')
    jmlPlots_cut_PWIAnorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                        pm_range=(3.469446951953614e-18,1.2),
                                        pm_bins=30)
    jmlPlots_cut_PWIAnorad.make_histos()
    jmlPlots_cut_PWIAnorad.make_PS_histos()
    if calc_vertex_kinematics:
        jmlPlots_cut_PWIAnorad.make_avgKin_histos()
    
    jmlPlots_cut_PWIArad = yield_plots(F.JML_PWIArad,use_cuts=True,
                                          data_type='SIMC')
    jmlPlots_cut_PWIArad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                      pm_range=(3.469446951953614e-18,1.2),
                                      pm_bins=30)
    jmlPlots_cut_PWIArad.make_histos()
# Paris
if F.load_PAR:
    parPlots_cut_FSIrad = yield_plots(F.PAR_FSIrad,use_cuts=True,
                                       data_type='SIMC')
    parPlots_cut_FSIrad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                     pm_range=(3.469446951953614e-18,1.2),
                                     pm_bins=30)
    parPlots_cut_FSIrad.make_histos()
    
    if calc_vertex_kinematics:
        parPlots_cut_FSInorad = yield_plots(F.PAR_FSInorad,use_cuts=True,
                                             data_type='SIMC',calc_avg_kin=True)
    else:
        parPlots_cut_FSInorad = yield_plots(F.PAR_FSInorad,use_cuts=True,
                                             data_type='SIMC')
    parPlots_cut_FSInorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                       pm_range=(3.469446951953614e-18,1.2),
                                       pm_bins=30)
    parPlots_cut_FSInorad.make_histos()
    parPlots_cut_FSInorad.make_PS_histos()
    if calc_vertex_kinematics:
        parPlots_cut_FSInorad.make_avgKin_histos()
    
    if calc_vertex_kinematics:
        parPlots_cut_PWIAnorad = yield_plots(F.PAR_PWIAnorad,use_cuts=True,
                                              data_type='SIMC',calc_avg_kin=True)
    else:
        parPlots_cut_PWIAnorad = yield_plots(F.PAR_PWIAnorad,use_cuts=True,
                                              data_type='SIMC')
    parPlots_cut_PWIAnorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                        pm_range=(3.469446951953614e-18,1.2),
                                        pm_bins=30)
    parPlots_cut_PWIAnorad.make_histos() 
    parPlots_cut_PWIAnorad.make_PS_histos()
    if calc_vertex_kinematics:
        parPlots_cut_PWIAnorad.make_avgKin_histos()
    
    parPlots_cut_PWIArad = yield_plots(F.PAR_PWIArad,use_cuts=True,
                                          data_type='SIMC')
    parPlots_cut_PWIArad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                        pm_range=(3.469446951953614e-18,1.2),
                                        pm_bins=30)
    parPlots_cut_PWIArad.make_histos()
    
# V18
if F.load_V18:
    v18Plots_cut_FSIrad = yield_plots(F.V18_FSIrad,use_cuts=True,
                                       data_type='SIMC')
    v18Plots_cut_FSIrad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                     pm_range=(3.469446951953614e-18,1.2),
                                     pm_bins=30)
    v18Plots_cut_FSIrad.make_histos()
    
    v18Plots_cut_FSInorad = yield_plots(F.V18_FSInorad,use_cuts=True,
                                         data_type='SIMC')
    v18Plots_cut_FSInorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                       pm_range=(3.469446951953614e-18,1.2),
                                       pm_bins=30)
    v18Plots_cut_FSInorad.make_histos()
    v18Plots_cut_FSInorad.make_PS_histos()
    
    v18Plots_cut_PWIAnorad = yield_plots(F.V18_PWIAnorad,use_cuts=True,
                                          data_type='SIMC')
    v18Plots_cut_PWIAnorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                        pm_range=(3.469446951953614e-18,1.2),
                                        pm_bins=30)
    v18Plots_cut_PWIAnorad.make_histos() 
    v18Plots_cut_PWIAnorad.make_PS_histos()
# CD-Bonn
if F.load_CDB:
    cdbPlots_cut_FSIrad = yield_plots(F.CDB_FSIrad,use_cuts=True,
                                       data_type='SIMC')
    cdbPlots_cut_FSIrad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                     pm_range=(3.469446951953614e-18,1.2),
                                     pm_bins=30)
    cdbPlots_cut_FSIrad.make_histos()
    
    cdbPlots_cut_FSInorad = yield_plots(F.CDB_FSInorad,use_cuts=True,
                                         data_type='SIMC')
    cdbPlots_cut_FSInorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                       pm_range=(3.469446951953614e-18,1.2),
                                       pm_bins=30)
    cdbPlots_cut_FSInorad.make_histos()
    cdbPlots_cut_FSInorad.make_PS_histos()
    
    cdbPlots_cut_PWIAnorad = yield_plots(F.CDB_PWIAnorad,use_cuts=True,
                                          data_type='SIMC')
    cdbPlots_cut_PWIAnorad.set_hist_lim(thrq_range=(0,190),thrq_bins=19,
                                        pm_range=(3.469446951953614e-18,1.2),
                                        pm_bins=30)
    cdbPlots_cut_PWIAnorad.make_histos()  
    cdbPlots_cut_PWIAnorad.make_PS_histos()                       

print('Finished creating plots...') 
#%%
# manipulate histograms to get yields and xsec
# use combine_histos() function to normalize and combine histograms from 
# different runs, the efficiency normalization is done per run by another 
# function normalize_histos() and calculates the total charge, then the 
# charge normalization is done at the end and that histo is returned.

# set filename here
save = False
filename = f'yield_n_xsec/{F.setting}_data_PAR_histos.pcl'

print('** Going to get yields and xsec...')

if F.load_data:
    many = F.data.many
    
    sets = db.retrieve('deuteron_db.db', 'run,data_set', 'RUN_LIST_UPDATED', 
        where= f"setting=\'{F.setting}\'")
    if F.setting == 'pm_120' or F.setting == 'pm_580':
        many = {'set_1':[x[0]for x in sets if x[1] == 'set_1'],
                'set_2':[x[0]for x in sets if x[1] == 'set_2']}
    else:
        many = {'set_1':[x[0]for x in sets if x[1] == 'set_1']}
    
    # data raw counts
    data_raw_h = combine_histos(dataPlots_nocut.histos_2D['Pm_vs_th_rq'],
                                many)
    
    # not normalized yield 
    yield_h = combine_histos(dataPlots_cut.histos_2D['Pm_vs_th_rq'],
                                   many)
    
    # normalized yield
    yield_norm_h = comb_n_norm_histos(dataPlots_cut.histos_2D['Pm_vs_th_rq'],
                                            many)
    
    norm_yield_proj = {}
    raw_yield_proj = {} ; raw_yield_counts = {}
    relative_error = {}
    for s in yield_norm_h:
        # Not Normalized and Norm Yield Projections
        norm_yield_proj[s] = project_2D(yield_norm_h[s],
                                           proj_title='$\\theta_{rq}$') 
        
        raw_yield_proj[s], raw_yield_counts[s] = project_2D(yield_h[s], 
                                                            get_counts=True,
                                                  proj_title='$\\theta_{rq}$')
        # relative error projections
        temp = {}
        for h in raw_yield_proj[s]:
            stats_err = np.sqrt(h.bin_content)
            rel_err = (stats_err/h.bin_content)*100 # in [%]
            
            temp[h.title] = rel_err
            relative_error[s] = temp

# JML SIMC plots
if F.load_JML:
    # Phase Space and PS projection
    JML_FSI_PS = jmlPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq_PS']
    JML_FSI_PS_proj = project_2D(JML_FSI_PS, proj_title='$\\theta_{rq}$')
    
    JML_PWIA_PS = jmlPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq_PS']
    JML_PWIA_PS_proj = project_2D(JML_PWIA_PS, proj_title='$\\theta_{rq}$')
    
    # JML SIMC FSI Cross Section
    JML_fsi_xsec_h = jmlPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                                                                JML_FSI_PS   
    # JML FSI Cross Section Projection
    JML_fsi_xsec_proj = project_2D(JML_fsi_xsec_h,
                                                proj_title='$\\theta_{rq}$')
    # JML SIMC PWIA Cross Section
    JML_pwia_xsec_h = jmlPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                                                                JML_PWIA_PS
    # JML PWIA Cross Section Projection
    JML_pwia_xsec_proj = project_2D(JML_pwia_xsec_h,
                                                proj_title='$\\theta_{rq}$')
     
    # JML Yield Projections                   
    JML_FSI_norad_proj = project_2D(jmlPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                                    proj_title='$\\theta_{rq}$') 
    
    JML_FSI_rad_proj = project_2D(jmlPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                                  proj_title='$\\theta_{rq}$')
    
    JML_PWIA_norad_proj = project_2D(jmlPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                                     proj_title='$\\theta_{rq}$') 

    JML_PWIA_rad_proj = project_2D(jmlPlots_cut_PWIArad.histos_2D['Pm_vs_th_rq'],
                                   proj_title='$\\theta_{rq}$')   
    
### Cross Section Calculation using JML #######################################
    
    # JML radiative correction histogram Y_norad/Y_rad
    JML_FSI_ratio_h = jmlPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                    jmlPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq']
    JML_FSI_radcorr_proj = project_2D(JML_FSI_ratio_h,
                                      proj_title='$\\theta_{rq}$')
    
    JML_PWIA_ratio_h = jmlPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                    jmlPlots_cut_PWIArad.histos_2D['Pm_vs_th_rq'] 
    JML_PWIA_radcorr_proj = project_2D(JML_PWIA_ratio_h,
                                      proj_title='$\\theta_{rq}$')

    if F.load_data:
        # Data Cross Section for each data set
        JML_yield_FSI_radcorr_h = {}
        JML_norm_yield_FSI_radcorr_proj = {}
        JML_yield_PWIA_radcorr_h = {}
        JML_norm_yield_PWIA_radcorr_proj = {}
        JML_FSI_dataXsec_h = {}
        JML_FSI_dataXsec_proj = {}
        JML_PWIA_dataXsec_h = {} 
        JML_PWIA_dataXsec_proj = {}
        for s in yield_norm_h:
            # Radiative corrected yield
            JML_yield_FSI_radcorr_h[s] = JML_FSI_ratio_h*yield_norm_h[s]
            JML_norm_yield_FSI_radcorr_proj[s] = project_2D(JML_yield_FSI_radcorr_h[s],
                                                   proj_title='$\\theta_{rq}$')
            
            JML_yield_PWIA_radcorr_h[s] = JML_PWIA_ratio_h*yield_norm_h[s]
            JML_norm_yield_PWIA_radcorr_proj[s] = project_2D(JML_yield_PWIA_radcorr_h[s],
                                               proj_title='$\\theta_{rq}$')
            
            # Data Differential Cross Section
            JML_FSI_dataXsec_h[s] = JML_yield_FSI_radcorr_h[s]/JML_FSI_PS 
            JML_FSI_dataXsec_proj[s] = project_2D(JML_FSI_dataXsec_h[s],proj_title='$\\theta_{rq}$')
            
            JML_PWIA_dataXsec_h[s] = JML_yield_PWIA_radcorr_h[s]/JML_PWIA_PS
            JML_PWIA_dataXsec_proj[s] = project_2D(JML_PWIA_dataXsec_h[s],proj_title='$\\theta_{rq}$')
                 
# PAR SIMC plots
if F.load_PAR:
    # Phase Space and PS projection
    PAR_FSI_PS = parPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq_PS']
    PAR_FSI_PS_proj = project_2D(PAR_FSI_PS, proj_title='$\\theta_{rq}$')
    
    PAR_PWIA_PS = parPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq_PS']
    PAR_PWIA_PS_proj = project_2D(PAR_PWIA_PS, proj_title='$\\theta_{rq}$')
    
    # PAR SIMC FSI Cross Section
    PAR_fsi_xsec_h = parPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                                                                PAR_FSI_PS   
    # PAR FSI Cross Section Projection
    PAR_fsi_xsec_proj = project_2D(PAR_fsi_xsec_h,
                                                proj_title='$\\theta_{rq}$')
    # PAR SIMC PWIA Cross Section
    PAR_pwia_xsec_h = parPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                                                                PAR_PWIA_PS
    
    # PAR PWIA Cross Section Projection
    PAR_pwia_xsec_proj = project_2D(PAR_pwia_xsec_h,
                                                proj_title='$\\theta_{rq}$')
     
    # PAR Yield Projections                   
    PAR_FSI_norad_proj = project_2D(parPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                                    proj_title='$\\theta_{rq}$') 
    
    PAR_FSI_rad_proj = project_2D(parPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                                  proj_title='$\\theta_{rq}$')
    
    PAR_PWIA_norad_proj = project_2D(parPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                                     proj_title='$\\theta_{rq}$')
    
    PAR_PWIA_rad_proj = project_2D(parPlots_cut_PWIArad.histos_2D['Pm_vs_th_rq'],
                                     proj_title='$\\theta_{rq}$')
    
### Cross Section Calculation using PAR #######################################
    
    # PAR radiative correction histogram Y_norad/Y_rad
    PAR_FSI_ratio_h = parPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                    parPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq']
    PAR_FSI_radcorr_proj = project_2D(PAR_FSI_ratio_h,
                                      proj_title='$\\theta_{rq}$')
    
    PAR_PWIA_ratio_h = parPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                    parPlots_cut_PWIArad.histos_2D['Pm_vs_th_rq'] 
    PAR_PWIA_radcorr_proj = project_2D(PAR_PWIA_ratio_h,
                                      proj_title='$\\theta_{rq}$')

    if F.load_data:
        # Data Cross Section for each data set
        PAR_yield_FSI_radcorr_h = {}
        PAR_norm_yield_FSI_radcorr_proj = {}
        PAR_yield_PWIA_radcorr_h = {}
        PAR_norm_yield_PWIA_radcorr_proj = {}
        PAR_FSI_dataXsec_h = {}
        PAR_FSI_dataXsec_proj = {}
        PAR_PWIA_dataXsec_h = {} 
        PAR_PWIA_dataXsec_proj = {}
        for s in yield_norm_h:
            # Radiative corrected yield
            PAR_yield_FSI_radcorr_h[s] = PAR_FSI_ratio_h*yield_norm_h[s]
            PAR_norm_yield_FSI_radcorr_proj[s] = project_2D(PAR_yield_FSI_radcorr_h[s],
                                                   proj_title='$\\theta_{rq}$')
            
            PAR_yield_PWIA_radcorr_h[s] = PAR_PWIA_ratio_h*yield_norm_h[s]
            PAR_norm_yield_PWIA_radcorr_proj[s] = project_2D(PAR_yield_PWIA_radcorr_h[s],
                                               proj_title='$\\theta_{rq}$')
            
            # Data Differential Cross Section
            PAR_FSI_dataXsec_h[s] = PAR_yield_FSI_radcorr_h[s]/PAR_FSI_PS 
            PAR_FSI_dataXsec_proj[s] = project_2D(PAR_FSI_dataXsec_h[s],proj_title='$\\theta_{rq}$')
            
            PAR_PWIA_dataXsec_h[s] = PAR_yield_PWIA_radcorr_h[s]/PAR_PWIA_PS
            PAR_PWIA_dataXsec_proj[s] = project_2D(PAR_PWIA_dataXsec_h[s],proj_title='$\\theta_{rq}$')

# V18 SIMC plots
if F.load_V18:
    # Phase Space and PS projection
    V18_FSI_PS = v18Plots_cut_FSInorad.histos_2D['Pm_vs_th_rq_PS']
    V18_FSI_PS_proj = project_2D(V18_FSI_PS, proj_title='$\\theta_{rq}$')
    
    V18_PWIA_PS = v18Plots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq_PS']
    V18_PWIA_PS_proj = project_2D(V18_PWIA_PS, proj_title='$\\theta_{rq}$')
    
    # V18 SIMC FSI Cross Section
    V18_fsi_xsec_h = v18Plots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                                                                V18_FSI_PS   
    # V18 FSI Cross Section Projection
    V18_fsi_xsec_proj = project_2D(V18_fsi_xsec_h,
                                                proj_title='$\\theta_{rq}$')
    # V18 SIMC PWIA Cross Section
    V18_pwia_xsec_h = v18Plots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                                                                V18_PWIA_PS
    
    # V18 PWIA Cross Section Projection
    V18_pwia_xsec_proj = project_2D(V18_pwia_xsec_h,
                                                proj_title='$\\theta_{rq}$')
     
    # V18 Yield Projections                   
    V18_FSI_norad_proj = project_2D(v18Plots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                                    proj_title='$\\theta_{rq}$') 
    
    V18_FSI_rad_proj = project_2D(v18Plots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                                  proj_title='$\\theta_{rq}$')
    
    V18_PWIA_norad_proj = project_2D(v18Plots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                                     proj_title='$\\theta_{rq}$') 
# CDB SIMC plots
if F.load_CDB:
    # Phase Space and PS projection
    CDB_FSI_PS = cdbPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq_PS']
    CDB_FSI_PS_proj = project_2D(CDB_FSI_PS, proj_title='$\\theta_{rq}$')
    
    CDB_PWIA_PS = cdbPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq_PS']
    CDB_PWIA_PS_proj = project_2D(CDB_FSI_PS, proj_title='$\\theta_{rq}$')
    
    # CDB SIMC FSI Cross Section
    CDB_fsi_xsec_h = cdbPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq']/\
                                                                CDB_FSI_PS   
    # CDB FSI Cross Section Projection
    CDB_fsi_xsec_proj = project_2D(CDB_fsi_xsec_h,
                                                proj_title='$\\theta_{rq}$')
    # CDB SIMC PWIA Cross Section
    CDB_pwia_xsec_h = cdbPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq']/\
                                                                CDB_PWIA_PS
    
    # CDB PWIA Cross Section Projection
    CDB_pwia_xsec_proj = project_2D(CDB_pwia_xsec_h,
                                                proj_title='$\\theta_{rq}$')
     
    # CDB Yield Projections                   
    CDB_FSI_norad_proj = project_2D(cdbPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                                    proj_title='$\\theta_{rq}$') 
    
    CDB_FSI_rad_proj = project_2D(cdbPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                                  proj_title='$\\theta_{rq}$')
    
    CDB_PWIA_norad_proj = project_2D(cdbPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                                     proj_title='$\\theta_{rq}$')    
# get commissioning data Cross Section and proj

if F.load_cd:
    if F.setting == 'pm_120':
        pm80_Xsec_h = F.pm80_data.histos['H_data2DXsec']
        pm80_xsec_proj = project_2D(pm80_Xsec_h,proj_title='$\\theta_{rq}$')
    elif F.setting == 'pm_580':
        h1 = F.pm580_comm_data1.histos['H_data2DXsec']
        h2 = F.pm580_comm_data2.histos['H_data2DXsec']
        # h = {'data1':h1,'data2':h2}
        # pm580_comm_Xsec_h = combine_histos(h,['data1','data2'])
        pm580_comm_Xsec_h = h1
        pm580_comm_xsec_proj = project_2D(pm580_comm_Xsec_h,proj_title='$\\theta_{rq}$')
    elif F.setting == 'pm_800' or F.setting == 'pm_900':
        h1 = F.pm580_comm_data1.histos['H_data2DXsec']
        h2 = F.pm580_comm_data2.histos['H_data2DXsec']
        # h = {'data1':h1,'data2':h2}
        # pm580_comm_Xsec_h = combine_histos(h,['data1','data2'])
        pm580_comm_Xsec_h = h1
        pm580_comm_xsec_proj = project_2D(pm580_comm_Xsec_h,proj_title='$\\theta_{rq}$')
        h1 = F.pm750_data1.histos['H_data2DXsec']
        h2 = F.pm750_data2.histos['H_data2DXsec']
        h3 = F.pm750_data3.histos['H_data2DXsec']
        # h = {'data1':h1,'data2':h2,'data3':h3}
        # pm750_Xsec_h = combine_histos(h,['data1','data2','data3'])
        pm750_Xsec_h = h1
        pm750_xsec_proj = project_2D(pm750_Xsec_h,proj_title='$\\theta_{rq}$')

# save histograms to a pickle file 
# define dict of histograms to save
hist = {}
if F.load_data:
    hist.update({'data_raw': data_raw_h,
            f'{F.setting}_yield': yield_h,
            f'{F.setting}_yield_norm': yield_norm_h,   
            f'{F.setting}_yield_norm_proj': norm_yield_proj,
            f'{F.setting}_yield_proj': raw_yield_proj,
            f'{F.setting}_yield_counts':raw_yield_counts,
            f'{F.setting}_Relative_Stats_Error': relative_error
            })

if F.load_JML:
    hist.update({f'{F.setting}_JML_FSI_PS': JML_FSI_PS,
                 f'{F.setting}_JML_FSI_PS_proj': JML_FSI_PS_proj,
                 f'{F.setting}_JML_PWIA_PS': JML_PWIA_PS,
                 f'{F.setting}_JML_PWIA_PS_proj': JML_PWIA_PS_proj,
                 f'{F.setting}_JML_fsi_xsec': JML_fsi_xsec_h,
                 f'{F.setting}_JML_fsi_xsec_proj': JML_fsi_xsec_proj,
                 f'{F.setting}_JML_pwia_xsec': JML_pwia_xsec_h,
                 f'{F.setting}_JML_pwia_xsec_proj': JML_pwia_xsec_proj,
                 f'{F.setting}_JML_Yield_FSI_norad':jmlPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_JML_Yield_FSI_norad_proj':JML_FSI_norad_proj,
                 f'{F.setting}_JML_Yield_FSI_rad':jmlPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_JML_Yield_FSI_rad_proj':JML_FSI_rad_proj,
                 f'{F.setting}_JML_Yield_PWIA_norad':jmlPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_JML_Yield_PWIA_norad_proj':JML_PWIA_norad_proj,
                 f'{F.setting}_JML_Yield_PWIA_rad':jmlPlots_cut_PWIArad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_JML_Yield_PWIA_rad_proj':JML_PWIA_rad_proj,
                 f'{F.setting}_JML_FSI_radcorr_factor':JML_FSI_ratio_h,
                 f'{F.setting}_JML_FSI_radcorr_factor_proj':JML_FSI_radcorr_proj,
                 f'{F.setting}_JML_PWIA_radcorr_factor':JML_PWIA_ratio_h,
                 f'{F.setting}_JML_PWIA_radcorr_factor_proj':JML_PWIA_radcorr_proj
                 })
    if F.load_data:
        hist.update({f'{F.setting}_yield_JML_FSI_radcorr': JML_yield_FSI_radcorr_h,
            f'{F.setting}_yield_JML_FSI_radcorr_proj': JML_norm_yield_FSI_radcorr_proj,
            f'{F.setting}_yield_JML_PWIA_radcorr': JML_yield_PWIA_radcorr_h,
            f'{F.setting}_yield_JML_PWIA_radcorr_proj': JML_norm_yield_PWIA_radcorr_proj,
            f'{F.setting}_JML_FSI_data_xsec': JML_FSI_dataXsec_h,
            f'{F.setting}_JML_FSI_data_xsec_proj': JML_FSI_dataXsec_proj,
            f'{F.setting}_JML_PWIA_data_xsec': JML_PWIA_dataXsec_h,
            f'{F.setting}_JML_PWIA_data_xsec_proj': JML_PWIA_dataXsec_proj,
            })
        
if F.load_PAR:
    hist.update({f'{F.setting}_PAR_FSI_PS': PAR_FSI_PS,
                 f'{F.setting}_PAR_FSI_PS_proj': PAR_FSI_PS_proj,
                 f'{F.setting}_PAR_PWIA_PS': PAR_PWIA_PS,
                 f'{F.setting}_PAR_PWIA_PS_proj': PAR_PWIA_PS_proj,
                 f'{F.setting}_PAR_fsi_xsec': PAR_fsi_xsec_h,
                 f'{F.setting}_PAR_fsi_xsec_proj': PAR_fsi_xsec_proj,
                 f'{F.setting}_PAR_pwia_xsec': PAR_pwia_xsec_h,
                 f'{F.setting}_PAR_pwia_xsec_proj': PAR_pwia_xsec_proj,
                 f'{F.setting}_PAR_Yield_FSI_norad':parPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_PAR_Yield_FSI_norad_proj':PAR_FSI_norad_proj,
                 f'{F.setting}_PAR_Yield_FSI_rad':parPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_PAR_Yield_FSI_rad_proj':PAR_FSI_rad_proj,
                 f'{F.setting}_PAR_Yield_PWIA_norad':parPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_PAR_Yield_PWIA_norad_proj':PAR_PWIA_norad_proj,
                 f'{F.setting}_PAR_Yield_PWIA_rad':parPlots_cut_PWIArad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_PAR_Yield_PWIA_rad_proj':PAR_PWIA_rad_proj,
                 f'{F.setting}_PAR_FSI_radcorr_factor':PAR_FSI_ratio_h,
                 f'{F.setting}_PAR_FSI_radcorr_factor_proj':PAR_FSI_radcorr_proj,
                 f'{F.setting}_PAR_PWIA_radcorr_factor':PAR_PWIA_ratio_h,
                 f'{F.setting}_PAR_PWIA_radcorr_factor_proj':PAR_PWIA_radcorr_proj
                 })
    if F.load_data:
        hist.update({f'{F.setting}_yield_PAR_FSI_radcorr': PAR_yield_FSI_radcorr_h,
            f'{F.setting}_yield_PAR_FSI_radcorr_proj': PAR_norm_yield_FSI_radcorr_proj,
            f'{F.setting}_yield_PAR_PWIA_radcorr': PAR_yield_PWIA_radcorr_h,
            f'{F.setting}_yield_PAR_PWIA_radcorr_proj': PAR_norm_yield_PWIA_radcorr_proj,
            f'{F.setting}_PAR_FSI_data_xsec': PAR_FSI_dataXsec_h,
            f'{F.setting}_PAR_FSI_data_xsec_proj': PAR_FSI_dataXsec_proj,
            f'{F.setting}_PAR_PWIA_data_xsec': PAR_PWIA_dataXsec_h,
            f'{F.setting}_PAR_PWIA_data_xsec_proj': PAR_PWIA_dataXsec_proj,
            })
if F.load_V18:
    hist.update({f'{F.setting}_V18_Yield_FSI_norad':v18Plots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_V18_Yield_FSI_rad':v18Plots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_V18_Yield_PWIA_norad':v18Plots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_V18_Yield_FSI_norad_proj':V18_FSI_norad_proj,
                 f'{F.setting}_V18_Yield_FSI_rad_proj':V18_FSI_rad_proj,
                 f'{F.setting}_V18_Yield_PWIA_norad_proj':V18_PWIA_norad_proj,
                 f'{F.setting}_V18_fsi_xsec': V18_fsi_xsec_h,
                 f'{F.setting}_V18_pwia_xsec': V18_pwia_xsec_h,
                 f'{F.setting}_V18_fsi_xsec_proj': V18_fsi_xsec_proj,
                 f'{F.setting}_V18_pwia_xsec_proj': V18_pwia_xsec_proj,
                 f'{F.setting}_V18_FSI_PS': V18_FSI_PS,
                 f'{F.setting}_V18_FSI_PS_proj': V18_FSI_PS_proj})
if F.load_CDB:
    hist.update({f'{F.setting}_CDB_Yield_FSI_norad':cdbPlots_cut_FSInorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_CDB_Yield_FSI_rad':cdbPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_CDB_Yield_PWIA_norad':cdbPlots_cut_PWIAnorad.histos_2D['Pm_vs_th_rq'],
                 f'{F.setting}_CDB_Yield_FSI_norad_proj':CDB_FSI_norad_proj,
                 f'{F.setting}_CDB_Yield_FSI_rad_proj':CDB_FSI_rad_proj,
                 f'{F.setting}_CDB_Yield_PWIA_norad_proj':CDB_PWIA_norad_proj,
                 f'{F.setting}_CDB_fsi_xsec': CDB_fsi_xsec_h,
                 f'{F.setting}_CDB_pwia_xsec': CDB_pwia_xsec_h,
                 f'{F.setting}_CDB_fsi_xsec_proj': CDB_fsi_xsec_proj,
                 f'{F.setting}_CDB_pwia_xsec_proj': CDB_pwia_xsec_proj,
                 f'{F.setting}_CDB_FSI_PS': CDB_FSI_PS,
                 f'{F.setting}_CDB_FSI_PS_proj': CDB_FSI_PS_proj})

if F.load_cd:
    if F.setting == 'pm_120':
        hist.update({'pm80_xsec':pm80_Xsec_h,
                     'pm80_xsec_proj':pm80_xsec_proj})
    elif F.setting == 'pm_580':
        hist.update({'pm580_comm_xsec':pm580_comm_Xsec_h,
                     'pm580_comm_xsec_proj':pm580_comm_xsec_proj})
    elif F.setting == 'pm_800' or F.setting == 'pm_900':
        hist.update({'pm580_comm_xsec':pm580_comm_Xsec_h,
                     'pm580_comm_xsec_proj':pm580_comm_xsec_proj,
                     'pm750_xsec':pm750_Xsec_h,
                     'pm750_xsec_proj':pm750_xsec_proj})
if save:
    with open(filename, 'wb') as f:
        P.dump(hist,f)
    print(f'Saved {filename}\n')

#%% save 2D avg histograms into pickle file
if F.load_JML:
    save_histos(F.JML_FSInorad,jmlPlots_cut_FSInorad)
    save_histos(F.JML_PWIAnorad,jmlPlots_cut_PWIAnorad)
if F.load_PAR:
    save_histos(F.PAR_FSInorad,parPlots_cut_FSInorad)
    save_histos(F.PAR_PWIAnorad,parPlots_cut_PWIAnorad)
















#%% histogram manipulation to get xsec
#####################
# OLD CODE GRAVEYARD
#####################
###
# function to get the yield based on a y-axis projection,
#   in a specified bin range in x, the projection is done bin by bin and 
#   the integral of the projected histogram is taken from y_min to y_max
# the opposite is true for projection along x-axis, set the flag project_x 
# to True to invert yield projection 
# NOT USED
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


pm120_h_radcorr = FSI_ratio_h*pm120_h_norm

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

#%%
xbin_range = np.arange(0,120,3)
ybin_range = np.arange(0.04,.2,0.04)

data_yield = []
data_yield_err = []
simc_yield = []
simc_yield_err = []
for xbin in xbin_range:
    pm120_yield_proj = pm120_h.project_y(range = (xbin-3,xbin+3))
    yi = pm120_yield_proj.sum(0,.2)
    
    simc_yield_proj = simcPlots_cut_FSIrad.histos_2D['Pm_vs_th_rq'].project_y(range = (xbin-3,xbin+3))
    syi = simc_yield_proj.sum(0,.2)
    
    data_yield.append(yi[0])
    data_yield_err.append(np.sqrt(yi[0])/yi[0])
    simc_yield.append(syi[0])
    simc_yield_err.append(np.sqrt(syi[0])/syi[0])
   
data_yield = []
data_yield_err = []
data_xsec = []
simc_xsec = []
for ybin in ybin_range:
    pm120_yield_proj = pm120_h.project_x(range = (ybin-0.04,ybin+0.04))
    yi = pm120_yield_proj.sum(0,180)
    
    pm120_xsec_proj = pm120_h_Xsec.project_x(range = (ybin-0.04,ybin+0.04))
    xs = pm120_xsec_proj.sum(42,48)
    
    simc_xsec_proj = SIMC_pm120_h_FSIXsec.project_x(range = (ybin-0.04,ybin+0.04))
    sxs = simc_yield_proj.sum(42,48)
    
        
    data_yield.append(yi[0])
    data_yield_err.append(np.sqrt(yi[0])/yi[0])
    
    data_xsec.append(xs[0])
    simc_xsec.append(sxs[0])

B.plot_exp(ybin_range,data_xsec,data_yield_err)
B.plot_line(ybin_range,simc_xsec)


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