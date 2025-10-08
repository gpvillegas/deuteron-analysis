#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import style
#style.use('seaborn-v0_8-poster')
#import uproot as upr
#import re
import os

import root_util as R
import LT.box as B

#import sqlite3 as lite
import database_operations as db

from matplotlib import colormaps as cmp
"""
Analyze root files using uproot based tool: root_util
"""

#%% CLASS DEFINITIONS

# class for heep run information
class load_run:
    
    def __init__(self, run, dbfile = 'dbfile', load = True):
        self.run = run
        self.dbfile = dbfile
        if load:
            self.load_settings()
            
    def load_settings(self):
        self.kin_study, self.setting, self.file =  \
            db.retrieve(self.dbfile, 'kin_study, setting, filename', 
                        'ROOTfiles', where = f"run=\'{self.run}\'")[0]
 

# make_plots class definition
"""
Using make_plots:
    
Description: Plot variable *v* from runs in *runs_list* with 
            or without default cuts (see below) and calculate normalization
            factors for each run. 
    make_plots(v,branches,runs_list,apply_cuts=True,custom_cuts=None)
            
            v: str or list(str1,str2,...);
                variable to plot using names from hcana.
            
            branches: dict; 
                        dictionary with loaded branches from root
                        trees to plot i. e. where v is located.
            
            runs_list: int or list(int1,int2,...);
                        list of runs for which you wish to see v plotted.
             
            apply_cuts: bool = True(default);
                        do you wish to see minimal cuts on v?,
                        standard cuts applied by default are:
                            -15 < shms delta < 22
                            -10 < hms delta < 10
                            -3 < z_react_diff < 3
                            -0.03 < shms xptar(th) < 0.03 
                            -0.03 < shms yptar(ph) < 0.03 
                            -0.08 < hms xptar < 0.08 
                            -0.03 < hms yptar < 0.03 

            custom_cuts: dict =\
                        {(str)'var-to-cut':(tuple)(min,max) or value}; 
                            apply additional cuts to v by setting custom_cuts
                            to a dictionary entry or entries with the hcana var
                            name as key and a tuple with the cut range or number 
                            as value.
                        as value
 Methods:
    make_1Dhistos(histo_params=[(-100,100),100],fig_title=None)
        --> Make 1D histos of v
        histo_params: list(range=(tuple),bins=int), parameters of the histo
    
    make_2Dhistos(histo_params=[(-100,100),100,(-100,100),100],
                  fig_title='')
        --> Make 2D histos of loaded variables in v in the order [x,y]
        histo_params: list(xrange=(tuple),xbins=int,yrange=(tuple),ybins=int), 
        parameters of the histogram    
                
"""

class make_plots:
    def __init__(self,v,branches,runs_list,apply_cuts=True,custom_cuts={}):
        self.v = dict.fromkeys(cast_to_list(v))     # variables to plot
        self.branches = branches                    # dict where the loaded 
                                                    #  variables are located
        self.runs_list = cast_to_list(runs_list)    # list of runs to plot vars from
        self.custom_cuts = custom_cuts              # add more cuts
        self.cuts = dict.fromkeys(self.runs_list)
        self.norm = {}
        #self.SIMC = SIMC
       
        # define cuts & apply
        
        # define standard cuts
        if apply_cuts:
            std_cuts = {'H.gtr.dp':(-10.,10.),'P.gtr.dp':(-10.,22.),
                        'H.gtr.th':(-0.08,0.08),'H.gtr.ph':(-0.03,0.03),
                        'P.gtr.th':(-0.03,0.03),'P.gtr.ph':(-0.03,0.03)}
            
            s_cuts = make_cuts(self.runs_list,std_cuts,self.branches)
            
            # standard cuts with operations
            op_cuts = {}
            for run in self.runs_list:
                #load cut variables
                H_react_z = self.branches[run]['H.react.z']
                P_react_z = self.branches[run]['P.react.z']
                
                #do operations
                z_diff = abs(H_react_z-P_react_z)
                
                #define cuts
                z_diff_cut_max = z_diff < 3.
                z_diff_cut_min = z_diff > -3.
                z_diff_cut = z_diff_cut_min & z_diff_cut_max
                
                op_cuts[run] = z_diff_cut
        
        # custom cuts
        if self.custom_cuts:
            c_cuts = make_cuts(self.runs_list,self.custom_cuts,self.branches)
        
        # add up all cuts to apply to v then apply them
        if apply_cuts and self.custom_cuts:
            for run in self.runs_list:
                self.cuts[run] = s_cuts[run] & op_cuts[run] & c_cuts[run]
                
            var = load_branches(self.runs_list,self.v,self.branches,self.cuts)               
        elif apply_cuts:
            for run in self.runs_list:
                self.cuts[run] = s_cuts[run] & op_cuts[run]
                
            var = load_branches(self.runs_list,self.v,self.branches,self.cuts)
        elif self.custom_cuts:
            for run in self.runs_list:
                self.cuts[run] = c_cuts[run]
                
            var = load_branches(self.runs_list,self.v,self.branches,self.cuts)  
        else:
            var = load_branches(self.runs_list,self.v,self.branches)               
        
        self.v = var 
                   
        # calculate normalization factor for each run 
        for run in self.runs_list:
            (q, h_teff, p_teff, lt) =\
                db.retrieve(deutDB_name, 
                            'BCM4A_charge, HMS_TrkEff, SHMS_TrkEff, T6_tLT', 
                            'RUN_LIST', where = f"run=\'{run}\'")[0]
                
            self.norm[run] = 1/(q*h_teff*p_teff*lt)
# =============================================================================

    def make_1Dhistos(self,histo_params=[(-100,100),100],fig_title=None,
                      apply_norm=False,SIMC=None,plot=True,save_as={}):
        self.fig_title = fig_title
        self.histo_params = histo_params
        self.SIMC = SIMC
        
        self.fname = save_as
        
        # plot variable(s)
        h_list={}
        for var in self.v:
            fig = plt.figure(layout='constrained')
            if self.fig_title != None:
                fig.suptitle(self.fig_title, fontsize=16)
            else:
                fig.suptitle(var, fontsize=16)
            
            rows= 2
            tot= len(self.runs_list)
            cols = int(tot/2)
            if tot == 1:
                rows = 1
                cols = 1                
            i=1
            j=1
            h2_list={}
            for runs in self.runs_list:                
                v = self.v[var][runs]
                ax = plt.subplot(rows,cols,i)
                if apply_norm:
                    h = self.norm[runs]*B.histo(v,histo_params[0],histo_params[1],
                                title=f'Run {runs}', xlabel='', ylabel='',
                                calc_w2=True)
                else:
                    h = B.histo(v,histo_params[0],histo_params[1],
                                title=f'Run {runs}', xlabel='', ylabel='')
                
                h2_list[runs] = h
                    
                
                if plot:
                    h.plot(label='Data',edgecolor='white')
                    
                    if self.SIMC != None:
                        hS = self.SIMC.make_1Dhistos(
                            histo_params=self.histo_params,
                            fig_flag=False,apply_weights=True)
                    
                # h.plot(color='#203080', edgecolor='white',
                #         label='Data')    
                                    
                ax.autoscale(enable=True, axis='both', tight=None)
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title(f'Run {runs}')
                ax.legend()
          
                i+=1
                j+=1
                if tot == 1:
                    continue
                if i > (rows*cols):
                    i=1
                    if j == tot:
                        plt.figure(layout='constrained')
                        
            h_list[var]=h2_list
        
            if save_as:                
                B.pl.savefig(self.fname[var])
        
        return h_list
        
    def make_2Dhistos(self,histo_params=[(-100,100),100,(-100,100),100],
                      fig_title='',apply_norm=False,plot=True,save_as=""):
        self.fig_title = fig_title
        
        self.fname = save_as
        # plot variable(s)
        if plot:
            fig = plt.figure(layout='constrained')
            fig.suptitle(self.fig_title, fontsize=16)
        
        rows= 2
        tot= len(self.runs_list)
        cols = int(tot/2)
        if tot == 1:
            rows = 1
            cols = 1             
        i=1
        j=1
        plot2= list(self.v.keys())
        hx= plot2[0]
        hy= plot2[1]
        h_list={}
        for runs in self.runs_list:
            vx = self.v[hx][runs]
            vy = self.v[hy][runs]
            ax = plt.subplot(rows,cols,i)
            if apply_norm:
                h = self.norm*B.histo2d(vx,vy,[histo_params[0],histo_params[2]],
                          [histo_params[1],histo_params[3]],logz=True,
                          calc_w2=True)
            else:
                h = B.histo2d(vx,vy,[histo_params[0],histo_params[2]],
                              [histo_params[1],histo_params[3]],logz=True)
            h_list[runs] = h
            
            # plt.subplot(rows,cols,i)
            if plot:
                h.plot(colormap=cmp['viridis'])
                        
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title(f'Run {runs}')
            
            i+=1
            j+=1
            if tot == 1:
                continue
            if i > (rows*cols):
                i=1
                if j == tot:
                    plt.figure(layout='constrained')
                    
        if save_as:
            B.pl.savefig(self.fname)
            
        return h_list

"""
make_plots_SIMC : same as make_plots for SIMC variables
"""                

class make_plots_SIMC:
    def __init__(self,v,branches,apply_cuts=False,custom_cuts=None):
        self.v = dict.fromkeys(cast_to_list(v))     # variables to plot
        self.branches = branches                    # dict where the loaded 
                                                    #  variables are located
        self.custom_cuts = custom_cuts              # add more cuts
        
        # calculate weights
        n = self.calc_weights()  
        
        # define cuts & apply
        if apply_cuts:
            std_cuts = {'h_delta':(-10.,10.),'e_delta':(-10.,22.),
                        'h_xptar':(-0.08,0.08),'h_yptar':(-0.03,0.03),
                        'e_xptar':(-0.03,0.03),'e_yptar':(-0.03,0.03)}
            
            s_cuts = make_cuts([],std_cuts,self.branches)

        # custom cuts
        if self.custom_cuts:
            c_cuts = make_cuts([],self.custom_cuts,self.branches)
        
        # add up all cuts to apply to v then apply them
        if apply_cuts and self.custom_cuts:
            self.cuts = s_cuts & c_cuts
                
            var = load_branches([],self.v,self.branches,self.cuts)
            wei = n[self.cuts]               
        elif apply_cuts:
            self.cuts = s_cuts
                
            var = load_branches([],self.v,self.branches,self.cuts)
            wei = n[self.cuts]               
        elif self.custom_cuts:
            self.cuts = c_cuts
                
            var = load_branches([],self.v,self.branches,self.cuts)  
            wei = n[self.cuts]               
        else:
            var = load_branches([],self.v,self.branches)
            wei = n              
        
        self.v = var 
        self.weights = wei
      
    def make_1Dhistos(self,histo_params=[(-100,100),100],fig_title='',
                      fig_flag=True,plot=True,apply_weights=False):
        self.fig_title = fig_title
        
        # plot variable(s)
        h_list={}
        for var in self.v.keys():
            if plot and fig_flag:
                fig = plt.figure(layout='constrained')
                if self.fig_title != None:
                    fig.suptitle(self.fig_title, fontsize=16)
                else:
                    fig.suptitle(var, fontsize=16)
            
            v = self.v[var]
            if apply_weights:
                hS = B.histo(v,histo_params[0],histo_params[1],
                            weights=self.weights,density=False,
                            calc_w2=True)
            else:
                hS = B.histo(v,histo_params[0],histo_params[1])
            
            h_list[var] = hS
            
            if plot:
                # hS.plot(label='SIMC',color='#de425b',edgecolor='white')
                hS.plot(label='SIMC',filled=False,color='#de425b')
            
        return h_list  
        
    def make_2Dhistos(self,histo_params=[(-100,100),100,(-100,100),100],
                      fig_title='',plot=True,apply_weights=False):
        # plot variable(s)
        fig = plt.figure(layout='constrained')
        fig.suptitle(fig_title, fontsize=16)

        plot2= list(self.v.keys())
        hx= plot2[0]
        hy= plot2[1]
        vx = self.v[hx]
        vy = self.v[hy]
        if apply_weights:
            h = B.histo2d(vx,vy,[histo_params[0],histo_params[2]],
                      [histo_params[1],histo_params[3]],logz=True,
                      title='', xlabel='', ylabel='', weights=self.weights,
                      calc_w2=True)
        else:
            h = B.histo2d(vx,vy,[histo_params[0],histo_params[2]],
                          [histo_params[1],histo_params[3]],logz=True,
                          title='', xlabel='', ylabel='')
        if plot:
            h.plot(colormap=cmp['viridis'])

        return h
        
    def calc_weights(self):
        w = self.branches['Weight']
        nf = self.branches['Normfac']
        nevt = w.size
        
        return w*nf/nevt

'''
Either 
    1. DEUT_DATA(kin_study= '', (opt)setting='', branches_sel=[])
    2. DEUT_DATA(branches_sel=[],Runs=[])
'''

class DEUT_DATA:
    def __init__(self,kin_study=None,setting=None,branches_sel=[],
                 trees_sel=['T','TSP','TSH','E'],Runs=None,load=True):
        self.kin_study = kin_study
        self.setting = setting
        self.branches_sel = branches_sel
        self.trees_sel = trees_sel
        self.Runs = cast_to_list(Runs)
        self.ROOTfiles_path = ROOTfiles_DIR
                 
        if kin_study != None and setting != None:
            k = self.kin_study
            s = self.setting
            
            self.Runs = get_list(
                db.retrieve(deutDB_name, 'run', 'ROOTfiles', 
                            where= f"kin_study=\'{k}\' AND setting=\'{s}\'",
                            distinct=True))
            
            self.ROOTfiles = {}
            for run_num in self.Runs:   
               self.ROOTfiles[run_num] = load_run(run_num, deutDB_name)
            
            T = self.TREES(self.Runs,self.ROOTfiles,self.ROOTfiles_path,
                           self.trees_sel)               
            self.BRANCHES(self.Runs,T,self.branches_sel,load)
            
        elif kin_study != None and setting == None:
            k = self.kin_study
            
            self.Runs = get_list(
                db.retrieve(deutDB_name, 'run', 'ROOTfiles', 
                            where= f"kin_study=\'{k}\'",distinct=True))
           
            s = {}
            self.ROOTfiles = {}
            for run_num in self.Runs:   
               self.ROOTfiles[run_num] = load_run(run_num, deutDB_name)
               s[run_num] = self.ROOTfiles[run_num].setting
            
            self.setting = s
            
            T = self.TREES(self.Runs,self.ROOTfiles,self.ROOTfiles_path,
                           self.trees_sel)               
            self.BRANCHES(self.Runs,T,self.branches_sel,load)
        
        elif kin_study == None and setting != None:
            s = self.setting
            
            self.Runs = get_list(
                db.retrieve(deutDB_name, 'run', 'ROOTfiles', 
                            where= f"setting=\'{s}\'",distinct=True))
           
            k = {}
            self.ROOTfiles = {}
            for run_num in self.Runs:   
               self.ROOTfiles[run_num] = load_run(run_num, deutDB_name)
               k[run_num] = self.ROOTfiles[run_num].kin_study
            
            self.kin_study = k
            
            T = self.TREES(self.Runs,self.ROOTfiles,self.ROOTfiles_path,
                           self.trees_sel)               
            self.BRANCHES(self.Runs,T,self.branches_sel,load)    
        
        elif kin_study == None and setting == None and Runs != None:
            k = {}
            s = {}
            self.ROOTfiles = {}
            for run_num in self.Runs:   
               self.ROOTfiles[run_num] = load_run(run_num, deutDB_name)
               (k[run_num], s[run_num]) = (self.ROOTfiles[run_num].kin_study,
                                           self.ROOTfiles[run_num].setting)
            self.kin_study = k
            self.setting = s
            
            T = self.TREES(self.Runs,self.ROOTfiles,self.ROOTfiles_path,
                           self.trees_sel)               
            self.BRANCHES(self.Runs,T,self.branches_sel,load)
        
        else:
            print('USAGE INFO:')
            print('OPTION 1 --> DEUT_DATA(kin_study= '', (opt)setting='',\
                   branches_sel=[])')
            print('Select a *kin_study* (heep_coin, heep_singles, deep, lumi) \
                  or *setting* and group of branches to load using leaf names \
                      from HCANA')    
            print('OPTION 2 --> DEUT_DATA(branches_sel=[],Runs=[])')
            print('Specify *Runs* and group of branches to load')
            
            if kin_study == None and setting == None and Runs == None:
                raise TypeError('kin_study, setting or Run not specified')
            else:
                raise TypeError('Something went wrong.\n Check input types:\n\
                                kin_study = string,\n setting = string,\n\
                                    branches_sel = list of strings*,\n\
                                        Runs = int or list of ints')
            
    def TREES(self,r,root,p,t):
        self.Trees = {}
        
        for run_num in r:
            rootfile = p + root[run_num].file
            if os.path.isfile(rootfile) :
                self.Trees[run_num]= R.root_tree(rootfile, trees=t)
                self.Trees[run_num].load_trees()
            else:
                continue
        return self.Trees         
    
    def BRANCHES(self,r,t,bs,l=False):
        sel = bs
        
        self.Branches = {}
        if l:
            for run_num in r:
               self.Branches[run_num] = t[run_num].get_branches(sel)    
            
        return self.Branches  

'''
    DEUT_SIMC(kin_study= '', (opt)setting='', branches_sel=[])
'''

class DEUT_SIMC:
    def __init__(self,kin_study=None,setting=None,branches_sel=[],load=True):
        self.kin_study = kin_study
        self.setting = cast_to_list(setting)
        self.branches_sel = branches_sel
        self.ROOTfiles_path = SIMCworksim_DIR 
        
        if kin_study != None and setting != None:
            k = self.kin_study
            s = self.setting
             
            self.ROOTfile = {}
            for s in self.setting:   
                self.ROOTfile[s] = get_list(
                    db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                    where = f"kin_study=\'{k}\' AND setting =\'{s}\'"))[0]
           
            if load:
                T = self.TREES(self.setting,self.ROOTfile,self.ROOTfiles_path)               
                self.BRANCHES(self.setting,T,self.branches_sel)
            
        elif kin_study != None and setting == None:
            k = self.kin_study
            
            self.setting = get_list(
                db.retrieve(deutDB_name, 'setting', 'SIMC_ROOTfiles', 
                            where = f"kin_study=\'{k}\'",distinct=True))
            
            self.ROOTfile = {}
            for s in self.setting:   
                self.ROOTfile[s] = get_list(
                    db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                    where = f"kin_study=\'{k}\' AND setting =\'{s}\'"))[0]
           
            if load:
                T = self.TREES(self.setting,self.ROOTfile,self.ROOTfiles_path)               
                self.BRANCHES(self.setting,T,self.branches_sel)
        
        elif kin_study == None and setting != None:
            s = self.setting
            
            self.kin_study = get_list(
                db.retrieve(deutDB_name, 'kin_study', 'SIMC_ROOTfiles', 
                            where = f"setting=\'{s}\'",distinct=True))
            k = self.kin_study
            
            self.ROOTfile = {}
            for s in self.setting:   
                self.ROOTfile[s] = get_list(
                    db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                    where = f"kin_study=\'{k}\' AND setting =\'{s}\'"))[0]
           
            if load:
                T = self.TREES(self.setting,self.ROOTfile,self.ROOTfiles_path)               
                self.BRANCHES(self.setting,T,self.branches_sel)  
        
        else:
            print('USAGE INFO:')
            print('DEUT_SIMC(kin_study= '', (opt)setting='',\
                   branches_sel=[])')
            print('Select a *kin_study* (heep_coin, heep_singles, deep, lumi) \
                  or *setting* and group of branches to load using leaf names \
                      from SIMC')
            print('Depending on how many settings are chosen there will be 1\
                      or more root files loaded')          
                      
            if kin_study == None and setting == None:
                raise TypeError('kin_study or setting not specified')
            else:
                raise TypeError('Something went wrong.\n Check input types:\n\
                                kin_study = string,\n setting = string,\n\
                                    branches_sel = list of strings*')

    def TREES(self,s,root,p,l=False):
        self.Trees = {}
        
        for setting in s:
            rootfile = p + root[setting]
            if os.path.isfile(rootfile) :
                self.Trees[setting]= R.root_tree(rootfile, trees=['SNT'])
                self.Trees[setting].load_trees()
            else :
                continue
        return self.Trees         
    
    def BRANCHES(self,s,t,bs):
        sel = bs
        
        self.Branches = {}
        for setting in s:
           self.Branches[setting] = t[setting].get_branches(sel)    
            
        return self.Branches     
                                    
#%% FUNCTION DEFINITIONS
# To get rid of unhelpful groupings (lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0] 

def cut_eff(c):
    return c.sum()/len(c)

def cast_to_list(obj):
    if type(obj) == list:
        return obj
    else:
        return [obj]   
# Branch selection
# sel: selection of branches to load

def branches_types(b=None):
        # assess calibration status:
        # Reference Time branches
        hREF_branches = ['T.coin.hDCREF1_tdcMultiplicity','T.coin.hDCREF1_tdcTime',
                'T.coin.hDCREF1_tdcTimeRaw','T.coin.hDCREF2_tdcMultiplicity',
                'T.coin.hDCREF2_tdcTime','T.coin.hDCREF2_tdcTimeRaw',
                'T.coin.hDCREF3_tdcMultiplicity','T.coin.hDCREF3_tdcTime',
                'T.coin.hDCREF3_tdcTimeRaw','T.coin.hDCREF4_tdcMultiplicity',
                'T.coin.hDCREF4_tdcTime','T.coin.hDCREF4_tdcTimeRaw',
                'T.coin.hDCREF5_tdcMultiplicity','T.coin.hDCREF5_tdcTime',
                'T.coin.hDCREF5_tdcTimeRaw', 'T.coin.hT1_tdcMultiplicity',
                'T.coin.hT1_tdcTime','T.coin.hT1_tdcTimeRaw',
                'T.coin.hT2_tdcMultiplicity','T.coin.hT2_tdcTime',
                'T.coin.hT2_tdcTimeRaw','T.coin.hFADC_TREF_ROC1_adcMultiplicity',
                'T.coin.hFADC_TREF_ROC1_adcPulseTime',
                'T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw']
        
        pREF_branches = ['T.coin.pDCREF10_tdcMultiplicity','T.coin.pDCREF10_tdcTime',
                  'T.coin.pDCREF10_tdcTimeRaw','T.coin.pDCREF1_tdcMultiplicity',
                  'T.coin.pDCREF1_tdcTime','T.coin.pDCREF1_tdcTimeRaw',
                  'T.coin.pDCREF2_tdcMultiplicity','T.coin.pDCREF2_tdcTime',
                  'T.coin.pDCREF2_tdcTimeRaw','T.coin.pDCREF3_tdcMultiplicity',
                  'T.coin.pDCREF3_tdcTime','T.coin.pDCREF3_tdcTimeRaw',
                  'T.coin.pDCREF4_tdcMultiplicity','T.coin.pDCREF4_tdcTime',
                  'T.coin.pDCREF4_tdcTimeRaw','T.coin.pDCREF5_tdcMultiplicity',
                  'T.coin.pDCREF5_tdcTime','T.coin.pDCREF5_tdcTimeRaw',
                  'T.coin.pDCREF6_tdcMultiplicity','T.coin.pDCREF6_tdcTime',
                  'T.coin.pDCREF6_tdcTimeRaw','T.coin.pDCREF7_tdcMultiplicity',
                  'T.coin.pDCREF7_tdcTime','T.coin.pDCREF7_tdcTimeRaw',
                  'T.coin.pDCREF8_tdcMultiplicity','T.coin.pDCREF8_tdcTime',
                  'T.coin.pDCREF8_tdcTimeRaw','T.coin.pDCREF9_tdcMultiplicity',
                  'T.coin.pDCREF9_tdcTime','T.coin.pDCREF9_tdcTimeRaw',
                  'T.coin.pT1_tdcMultiplicity','T.coin.pT1_tdcTime',
                  'T.coin.pT1_tdcTimeRaw','T.coin.pT2_tdcMultiplicity',
                  'T.coin.pT2_tdcTime','T.coin.pT2_tdcTimeRaw',
                  'T.coin.pFADC_TREF_ROC2_adcMultiplicity',
                  'T.coin.pFADC_TREF_ROC2_adcPulseTime',
                  'T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw']
        
        
            
        # HODOSCOPE branches
        Phodo_branches = ['P.hod.beta','P.hod.betanotrack','P.hod.betachisqnotrack',
                         'P.hod.adctdc_offset','P.hod.2y.totNumTdcHits',
                         'P.hod.2y.totNumAdcHits','P.hod.2y.TrackXPos',
                         'P.hod.2y.TrackYPos','P.hod.2x.totNumTdcHits',
                         'P.hod.2x.totNumAdcHits','P.hod.2x.TrackXPos',
                         'P.hod.2x.TrackYPos','P.hod.1y.totNumTdcHits',
                         'P.hod.1y.totNumAdcHits','P.hod.1y.TrackXPos',
                         'P.hod.1y.TrackYPos','P.hod.1x.totNumTdcHits',
                         'P.hod.1x.totNumAdcHits','P.hod.1x.TrackXPos',
                         'P.hod.1x.TrackYPos']
        
        Hhodo_branches = ['H.hod.beta','H.hod.betanotrack','H.hod.betachisqnotrack',
                         'H.hod.adctdc_offset','H.hod.2y.totNumTdcHits',
                         'H.hod.2y.totNumAdcHits','H.hod.2y.TrackXPos',
                         'H.hod.2y.TrackYPos','H.hod.2x.totNumTdcHits',
                         'H.hod.2x.totNumAdcHits','H.hod.2x.TrackXPos',
                         'H.hod.2x.TrackYPos','H.hod.1y.totNumTdcHits',
                         'H.hod.1y.totNumAdcHits','H.hod.1y.TrackXPos',
                         'H.hod.1y.TrackYPos','H.hod.1x.totNumTdcHits',
                         'H.hod.1x.totNumAdcHits','H.hod.1x.TrackXPos',
                         'H.hod.1x.TrackYPos']
        
        # DRIFT CHAMBER branches
        Pdc_branches = ['P.dc.x_fp','P.dc.y_fp','P.dc.residualExclPlane',
                        'P.dc.2x2.dist', 'P.dc.2x1.dist', 'P.dc.2v2.dist', 
                        'P.dc.2v1.dist', 'P.dc.2u2.dist', 'P.dc.2u1.dist',
                        'P.dc.1x2.dist', 'P.dc.1x1.dist', 'P.dc.1v2.dist',
                        'P.dc.1v1.dist', 'P.dc.1u2.dist', 'P.dc.1u1.dist'] 
        
        Hdc_branches = ['H.dc.x_fp','H.dc.y_fp','H.dc.residualExclPlane',
                        'H.dc.2x2.dist', 'H.dc.2x1.dist', 'H.dc.2v2.dist', 
                        'H.dc.2v1.dist', 'H.dc.2u2.dist', 'H.dc.2u1.dist',
                        'H.dc.1x2.dist', 'H.dc.1x1.dist', 'H.dc.1v2.dist',
                        'H.dc.1v1.dist', 'H.dc.1u2.dist', 'H.dc.1u1.dist']
        
        # CALORIMETER branches
        Pcal_branches = ['P.cal.etottracknorm', 'P.cal.etotnorm', 'P.cal.ytrack',
                         'P.cal.xtrack']
        
        Hcal_branches = ['H.cal.etottracknorm', 'H.cal.etotnorm', 'H.cal.ytrack',
                         'H.cal.xtrack']
        
        # CERENKOV branches
        Phgcer_branches = ['P.hgcer.npeSum', 'P.hgcer.npe']
        
        Pngcer_branches = ['P.ngcer.npeSum', 'P.ngcer.npe']
        
        Hcer_branches = ['H.cer.npeSum', 'H.cer.npe']
        
        # KINEMATIC branches
        
        kin_branches= ['H.kin.secondary.Erecoil','H.kin.secondary.MMK',
                       'H.kin.secondary.MMp','H.kin.secondary.MMpi',
                       'H.kin.secondary.MandelS','H.kin.secondary.MandelT',
                       'H.kin.secondary.MandelU','H.kin.secondary.Mrecoil',
                       'H.kin.secondary.Prec_x','H.kin.secondary.Prec_y',
                       'H.kin.secondary.Prec_z','H.kin.secondary.emiss',
                       'H.kin.secondary.emiss_nuc','H.kin.secondary.ph_bq',
                       'H.kin.secondary.ph_xq','H.kin.secondary.phb_cm',
                       'H.kin.secondary.phx_cm','H.kin.secondary.pmiss',
                       'H.kin.secondary.pmiss_x','H.kin.secondary.pmiss_y',
                       'H.kin.secondary.pmiss_z','H.kin.secondary.px_cm',
                       'H.kin.secondary.t_tot_cm','H.kin.secondary.tb',
                       'H.kin.secondary.tb_cm','H.kin.secondary.th_bq',
                       'H.kin.secondary.th_xq','H.kin.secondary.thb_cm',
                       'H.kin.secondary.thx_cm','H.kin.secondary.tx',
                       'H.kin.secondary.tx_cm','H.kin.secondary.xangle',
                       'P.kin.primary.Q2','P.kin.primary.W','P.kin.primary.W2',
                       'P.kin.primary.epsilon','P.kin.primary.nu',
                       'P.kin.primary.omega','P.kin.primary.ph_q','P.kin.primary.q3m',
                       'P.kin.primary.q_x','P.kin.primary.q_y','P.kin.primary.q_z',
                       'P.kin.primary.scat_ang_deg','P.kin.primary.scat_ang_rad',
                       'P.kin.primary.th_q','P.kin.primary.x_bj']
        
        # GOLDEN TRACK branches
        
        gtr_branches= ['H.dc.gtrack_nsp','H.gtr.beta','H.gtr.dp','H.gtr.index',
                       'H.gtr.ok','H.gtr.p','H.gtr.ph','H.gtr.px','H.gtr.py',
                       'H.gtr.pz','H.gtr.th','H.gtr.x','H.gtr.y','P.dc.gtrack_nsp',
                       'P.gtr.beta','P.gtr.dp','P.gtr.index','P.gtr.ok','P.gtr.p',
                       'P.gtr.ph','P.gtr.px','P.gtr.py','P.gtr.pz','P.gtr.th',
                       'P.gtr.x','P.gtr.y']
        
        # REACTION VERTEX branches
        
        react_branches= ['H.react.ok','H.react.x','H.react.y','H.react.z',
                         'P.react.ok','P.react.x','P.react.y','P.react.z']
        
        # CUT VARIABLES (not included in the groups above)
        cut_branches= ['P.gtr.dp','H.gtr.dp','P.react.z','H.react.z',
                       'P.gtr.th','P.gtr.ph','H.gtr.th','H.gtr.ph']
        
        # optics variables: x/yfp, x/yptar, dp, x/ysieve
        optics_branches= ['P.dc.x_fp','P.dc.y_fp','P.dc.xp_fp','P.dc.yp_fp',
                          'H.dc.x_fp','H.dc.y_fp','H.dc.xp_fp','H.dc.yp_fp',
                          'P.gtr.th','P.gtr.ph','H.gtr.th','H.gtr.ph',
                          'P.extcor.xsieve','P.extcor.ysieve',
                          'H.extcor.xsieve','H.extcor.ysieve',
                          'P.gtr.dp','H.gtr.dp']
        
        branch_types= {'phodo':Phodo_branches,
                         'hhodo':Hhodo_branches,
                         'pdc':Pdc_branches,
                         'hdc':Hdc_branches,
                         'pcal':Pcal_branches,
                         'hcal':Hcal_branches,
                         'phgcer':Phgcer_branches,
                         'pngcer':Pngcer_branches,
                         'hcer':Hcer_branches,
                         'kin':kin_branches,
                         'gtr': gtr_branches,
                         'react': react_branches,
                         'cut':cut_branches,
                         'hREF':hREF_branches,
                         'pREF':pREF_branches,
                         'optics':optics_branches}

        branches_groups= {'calib':branch_types['phodo']+
                                   branch_types['hhodo']+
                                   branch_types['pdc']+
                                   branch_types['hdc']+
                                   branch_types['pcal']+
                                   branch_types['hcal']+
                                   branch_types['phgcer']+
                                   branch_types['pngcer']+
                                   branch_types['hcer'],
                          'kins':branch_types['kin']+
                                 branch_types['cut'],
                          'reftime':branch_types['hREF']+
                                    branch_types['pREF']}
        branches_groups.update(branch_types)
        
        list_b = cast_to_list(b)
        if list_b:
            branches_sel = []
            for name in list_b:
                branches_sel = branches_sel + branches_groups[name]        
            return branches_sel
        else:
            return []    

def load_branches(run=[],var={},bran={},cuts={}):
    
    flag = False
    if cuts != {} or cuts != []:
        flag = True
    
    if run:    
        for v in var:
            v_load = {}
            for r in run:
                if flag:
                    v_load[r] = bran[r][v][cuts[r]]
                else:
                    v_load[r] = bran[r][v] 
            var[v] = v_load
    else:       
        for v in var:
            v_load = []
            if flag:
                v_load = bran[v][cuts]
            else:
                v_load = bran[v] 
            var[v] = np.array(v_load)      
    return var

def make_cuts(run=[],cus_cut={},bran={}):
    c_var = dict.fromkeys(cus_cut.keys())
    c_var = load_branches(run,c_var,bran,{})
    c_cuts = {}
    
    if run:
        for r in run:
            c_max = []
            c_min = []
            c_val = []
            this_cut = []
            flag = True
    
            for cv in c_var:
                if type(cus_cut[cv]) != tuple:
                    c_val = c_var[cv][r] == cus_cut[cv]
                    this_cut = np.array(c_val)
                else:    
                    c_max =\
                        c_var[cv][r] < cus_cut[cv][1] 
                    c_min =\
                        c_var[cv][r] > cus_cut[cv][0]
                    this_cut = c_min & c_max
                    this_cut = np.array(this_cut)
                    
                if flag:
                    c_cuts[r] = this_cut
                    flag = False 
                else:
                    c_cuts[r] = c_cuts[r] & this_cut
    else:
        c_max = []
        c_min = []
        c_val = []
        this_cut = []
        flag = True     # used to initialize c_cuts 
                        # then apply consecutive cuts on top

        for cv in c_var:
            if type(cus_cut[cv]) != tuple:
                c_val = c_var[cv] == cus_cut[cv]
                this_cut = np.array(c_val)
            else:    
                c_max =\
                    c_var[cv] < cus_cut[cv][1] 
                c_min =\
                    c_var[cv] > cus_cut[cv][0]
                this_cut = c_min & c_max
                this_cut = np.array(this_cut)
                
            if flag:
                c_cuts = this_cut
                flag = False 
            else:
                c_cuts = c_cuts & this_cut    
    return c_cuts
#%% setting up directory locations

ROOTfiles_DIR= '/home/gvill/deuteron/ROOTfiles/prod/' # location of hcana root files
# ROOTfiles_DIR= '/home/gvill/deuteron/ROOTfiles/prod/heep_coin_optics_nooffsets_2018optsmatrx/' 
#ROOTfiles_DIR= '/home/gvill/deuteron/ROOTfiles/scalers/' # location of hcana root files
#ROOTfiles_DIR="/media/gvill/Gema's T7/ROOTfiles/"
# SIMCworksim_DIR= '/home/gvill/deuteron/hallc_simulations/worksim/oldOffsets/' 
SIMCworksim_DIR= '/home/gvill/deuteron/hallc_simulations/worksim/' # location of SIMC root files

# database name
deutDB_name= 'deuteron_db.db'

# gather table names = ['RUN_LIST', 'ROOTfiles', 'DataFiles', 'SIMC_ROOTfiles']
dbtable_names= db.get_list_of_tables(deutDB_name)

# gather row information
RUN_LIST_rows= db.get_table_information(deutDB_name, dbtable_names[0]) 
ROOTfiles_rows= db.get_table_information(deutDB_name, dbtable_names[1]) 
SIMC_ROOTfiles_rows= db.get_table_information(deutDB_name, dbtable_names[3]) 

#%% LOADING DEUTERON 2023 ROOTFILES
# BRANCHES ARE STORED IN THE pm120, pm580, pm800, pm900, DEUT_DATA objects.

# DATA ========================================================================
#pm120 = DEUT_DATA(kin_study='deep',setting='pm_120',branches_sel=branches_types('kins'))
#pm580 = DEUT_DATA(kin_study='deep',setting='pm_580',branches_sel=['kins'])
#pm800 = DEUT_DATA(kin_study='deep',setting='pm_800',branches_sel=['kins'])
#pm900 = DEUT_DATA(kin_study='deep',setting='pm_900',branches_sel=['kins'])
   
# SIMC ========================================================================
#deep = DEUT_SIMC(kin_study='deep')
#pm120_SIMC = DEUT_SIMC(kin_study='deep',setting='pm_120')
#pm120_SIMC = pm120_SIMC.Branches['pm_120']
#pm580_SIMC = DEUT_SIMC(kin_study='deep',setting='pm_580')
#pm800_SIMC = DEUT_SIMC(kin_study='deep',setting='pm_800')
#pm900_SIMC = DEUT_SIMC(kin_study='deep',setting='pm_900')