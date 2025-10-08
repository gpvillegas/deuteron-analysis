#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:49:00 2024

@author: gvill

deut_analysis_tools v2: reworked data classes
"""
import os
import root_util as R
import database_operations as db
import LT.box as B

import numpy as np

#%% setting up directory locations

# location of hcana root files
# ROOTfiles_DIR= '/home/gvill/deuteron/ROOTfiles/prod/deut18_sieve_slit/'
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step0_floorvalues_nooffsets/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step1_kf_offset/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step2_pf_offset/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step3_excludehighdelta/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step4_finetuning/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step5_kfoff_pfoffperset/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_offsetstudies_step6_pthetaoffset/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/deep_pm120/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_defaultpoffsets/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim/" 
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim2/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim2_pfoff/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim2_pfoff_hangleoffs/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim2_1_09966/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim2_1_pfofflinefit/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/deep_testing0/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim3_allfp/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim3_thpoff/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim3_nodeltam8/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim3_theoff/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_hyptaroff/"
# ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_hyptaroff_hthetacentraloff/"
ROOTfiles_DIR= "/media/gvill/Gema's T7/ROOTfiles/heepcoin_deltaoptim3_alloffsets/"


# location of SIMC root files
# SIMCworksim_DIR= '/home/gvill/deuteron/hallc_simulations/worksim/'
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_offsetstudies_step0_floorvalues_nooffsets/"  
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_offsetstudies_step1_kf_offset/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_offsetstudies_step2_pf_offset/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_offsetstudies_step3_excludehighdelta/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_offsetstudies_step4_finetuning/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_offsetstudies_step5_kfoff_pfoffperset/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_defaultpoffsets/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim2/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim2_1_pfofflinefit/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim2_pfoff_hangleoffs/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim2_1/"
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim2_1_kfoffpersett/"
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/deep_testing0/" 
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim3/"
SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim3_thpoff/"
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_deltaoptim3_theoff/"
# SIMCworksim_DIR= "/media/gvill/Gema's T7/ROOTfiles/worksim/heepcoin_hyptaroff_hthetacentraloff/"


# database name
deutDB_name= 'deuteron_db.db'

#%% Helpful Functions
def cast_to_list(obj):
    if type(obj) == list:
        return obj
    else:
        return [obj]
    
# To get rid of unhelpful groupings (lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0]
    
def get_eff_norm(run,run_type='coin'):
    if run_type == 'coin':
        h_teff, p_teff, lt =\
            db.retrieve('deuteron_db.db', 
                        'HMS_TrkEff, SHMS_TrkEff, T6_tLT', 
                        'RUN_LIST', where = f"run=\'{run}\'")[0]
        norm = h_teff*p_teff*lt
    elif run_type == 'singles':
        p_teff, lt =\
            db.retrieve('deuteron_db.db', 
                        'SHMS_TrkEff, T1_tLT', 
                        'RUN_LIST', where = f"run=\'{run}\'")[0]
        norm = p_teff*lt
    else:
        print('Need run_type to know what trigger efficiency is needed.')
        norm = 1        
    # print(h_teff,p_teff,lt)    
    return norm

def get_charge_norm(run):
    q =\
        get_list(db.retrieve('deuteron_db.db', 
                    'BCM4A_charge', 
                    'RUN_LIST', where = f"run=\'{run}\'"))
    
    return q

def calc_weights(simc_branches):
    w = simc_branches['Weight']
    nf = simc_branches['Normfac']
    nevt = w.size
    
    return w*nf/nevt

def calc_weights_PS(simc_branches):
    nf = simc_branches['Normfac']
    nevt = nf.size
    
    return nf/nevt  

def get_coinTimeOffset(coinTime):

    CTime = coinTime  
    CTime_histo_0 = B.histo(CTime,range=(-100,100),bins=200,
                              title = 'CTime no cuts')
    CTime_histo_0.fit(plot_fit=False)
    
    # CTime_histo_0.plot()
    # CTime_histo_0.plot_fit()
    
    mu_0 = CTime_histo_0.mean.value    
    sigma_0 = abs(CTime_histo_0.sigma.value)
    
    # x1_min = mu_0 - 2.*sigma_0
    # x1_max = mu_0 + 2.*sigma_0
    # CTime_histo_1 = B.histo(CTime,range=(x1_min,x1_max),bins=200,
    #                           title = 'CTime no cuts')    
    # CTime_histo_1.fit(plot_fit=False)   
    # mu_1 = CTime_histo_1.mean.value
    # sigma_1 = CTime_histo_1.sigma.value
    
    x_min = mu_0 - 2.*sigma_0
    x_max = mu_0 + 2.*sigma_0
    
    CTime_histo = B.histo(CTime,range=(x_min,x_max),bins=200,
                              title = 'CTime no cuts')    
    CTime_histo.fit(plot_fit=False)
    
    # B.pl.figure()
    # CTime_histo.plot()
    # CTime_histo.plot_fit()
    
    mu = CTime_histo.mean.value    
    CTime_corr = CTime - mu    
    return CTime_corr
#%% DATA_INIT CLASS DEFINITION
class DATA_INIT:
    def __init__(self,data_type=None,kin_study='',setting='',
                 run = int(),select_branches={},select_trees=[],simc_type = '',
                 load=True, 
                 ROOTfiles_path=ROOTfiles_DIR, 
                 SIMC_ROOTfiles_path=SIMCworksim_DIR):
        
        self.dtype = data_type
        
#=============================================================================
#       'deut23_data' DEFINITION 
#=============================================================================
        
        if self.dtype == 'deut23_data':
            self.kin_study = kin_study
            self.setting = setting
            
            if not select_branches and load:
                print('No specific branches selected.',
                      'Loading all available branches!',
                      'This might take a while...',sep=' ')
            self.branches_sel = select_branches
             
            if not select_trees:
                print('Only data tree: T will be loaded.')
                self.trees_sel = ['T']
            else:    
                self.trees_sel = select_trees
            
            self.many = cast_to_list(run)       
            self.ROOTfiles_path = ROOTfiles_path
            
            # Based on the kinematic study, setting, or both, the code will
            # look in the db for the appropriate run numbers and their root 
            # file name.
            # Alternatively, specific run numbers could be given for the 
            # same end.
            # Then it will load the selected trees (the first one in the 
            # file if none chosen) and selected branches (all if an empty 
            # list is passed)
            
            # Case 1: kin_study and setting are given 
            # e.g. kin_study = 'heep_coin', setting = 'delta_scan_0'
            
            if kin_study and setting:
                k = self.kin_study
                s = self.setting
                
                self.many = get_list(
                    db.retrieve(deutDB_name, 'run', 'ROOTfiles', 
                        where= f"kin_study=\'{k}\' AND setting=\'{s}\'",
                        distinct=True))
                                
                try:
                    self.ROOTfiles = {}
                    for run_num in self.many:   
                       self.ROOTfiles[run_num] = get_list(
                           db.retrieve(deutDB_name, 'filename', 'ROOTfiles', 
                                       where = f"run=\'{run_num}\'"))
                except TypeError:
                    self.ROOTfiles = get_list(
                        db.retrieve(deutDB_name, 'filename', 'ROOTfiles', 
                                    where = f"run=\'{self.many}\'"))
                
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                               self.trees_sel) 
                ##### trying something new ####
                # 06/18 changed BRANCHES method to allow for a dict of
                # branches_sel where you specify branches from each tree.
                self.BRANCHES(self.many,T,self.branches_sel,load)
                
            # Case 2: only kin_study is given
            # e.g. kin_study = 'heep_coin'
              
            elif kin_study:
                k = self.kin_study
                
                self.many = get_list(
                    db.retrieve(deutDB_name, 'run', 'ROOTfiles', 
                        where= f"kin_study=\'{k}\'",
                        distinct=True))
                                
                try:
                    s = {}
                    self.ROOTfiles = {}
                    for run_num in self.many:   
                       self.ROOTfiles[run_num], s[run_num] =\
                           db.retrieve(deutDB_name, 
                                       'filename, setting', 'ROOTfiles', 
                                       where = f"run=\'{run_num}\'")[0]
                    self.setting = s       
                           
                except TypeError:
                    self.ROOTfiles, s =\
                           db.retrieve(deutDB_name, 
                                       'filename, setting', 'ROOTfiles', 
                                       where = f"run=\'{self.many}\'")[0]
                    self.setting = s       
                
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                               self.trees_sel)
                self.BRANCHES(self.many,T,self.branches_sel,load)
            
            # Case 3: only setting is given
            # e.g. setting = 'delta_scan_0'
            
            elif setting:
                s = self.setting
                
                self.many = get_list(
                    db.retrieve(deutDB_name, 'run', 'ROOTfiles', 
                        where= f"setting=\'{s}\'",
                        distinct=True))
                                
                try:
                    k = {}
                    self.ROOTfiles = {}
                    for run_num in self.many:   
                       self.ROOTfiles[run_num], k[run_num] =\
                           db.retrieve(deutDB_name, 
                                       'filename, kin_study', 'ROOTfiles', 
                                       where = f"run=\'{run_num}\'")[0]
                    self.kin_study = k 
                           
                except TypeError:
                    self.ROOTfiles[run_num], k[run_num] =\
                           db.retrieve(deutDB_name, 
                                       'filename, kin_study', 'ROOTfiles', 
                                       where = f"run=\'{self.many}\'")[0]
                    self.kin_study = k   
                
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                               self.trees_sel)               
                self.BRANCHES(self.many,T,self.branches_sel,load)

            # Case 4: specific runs are given
            # e.g. run = [20851,20870,20871] or run = 20851
            
            elif not kin_study and not setting and run != 0:
                self.many = run
                
                try:
                    k = {}
                    s = {}
                    self.ROOTfiles = {}
                    for run_num in self.many:   
                       self.ROOTfiles[run_num], k[run_num], s[run_num] =\
                           db.retrieve(
                               deutDB_name, 'filename, kin_study, setting',
                               'ROOTfiles', where = f"run=\'{run_num}\'")[0]
                    self.kin_study = k
                    self.setting = s
                           
                except TypeError:
                    self.ROOTfiles, k, s =\
                           db.retrieve(
                               deutDB_name, 'filename, kin_study, setting',
                               'ROOTfiles', where = f"run=\'{self.many}\'")[0]
                    self.kin_study = k
                    self.setting = s
                
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                                self.trees_sel)  
            
                self.BRANCHES(self.many,T,self.branches_sel,load)
            
            # Case 5: nothing is given, do you want ALL runs? maybe not...
            
            elif not kin_study and not setting and not run:
                print('======================================================')
                print(' NO kin_study, setting, or run SELECTED')
                print('======================================================')
                print('Run USAGE() for more information')
                 
#=============================================================================
#=============================================================================
#       'SIMC' DEFINITION 
#=============================================================================
            
        elif self.dtype == 'SIMC':
            self.kin_study = kin_study
            self.setting = setting
            self.simc_type = simc_type
            
            if not select_branches and load:
                print('No specific branches selected.',
                      'Loading all available branches!',
                      'This might take a while...',sep='\n')
            self.branches_sel = select_branches            
            
            if not select_trees:
                print('First tree: SNT will be loaded.')
                self.trees_sel = ['SNT']
            else:    
                self.trees_sel = select_trees
            
            self.many = setting    
            self.ROOTfiles_path = SIMC_ROOTfiles_path
            
            # Case 1: kin_study and setting are given 
            # e.g. kin_study = 'heep_coin', setting = 'delta_scan_0'
            
            if kin_study and setting and simc_type:
                
                if type(self.many) is list: 
                    self.ROOTfiles = {}
                    for sett in self.many:   
                       self.ROOTfiles[sett] = get_list(
                           db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                                       where = f"setting=\'{sett}\' AND"+\
                                        f" type=\'{self.simc_type}\'"))[0]
                else:    
                    self.ROOTfiles =get_list(
                        db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                                    where = f"setting=\'{self.many}\' AND"+\
                                        f" type=\'{self.simc_type}\'"))
                
                # print(self.many,'\n',self.ROOTfiles,'\n',
                #       self.ROOTfiles_path,'\n',self.trees_sel)
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                               self.trees_sel)
                # print(T)
                self.BRANCHES(self.many,T,self.branches_sel,load)
                
            # Case 2: only kin_study is given
            # e.g. kin_study = 'heep_coin'
              
            elif kin_study and simc_type:
                k = self.kin_study
                
                self.many = get_list(
                    db.retrieve(deutDB_name, 'setting', 'SIMC_ROOTfiles', 
                        where= f"kin_study=\'{k}\'",
                        distinct=True))
                self.setting = self.many
                
                if type(self.many) is list:
                    self.ROOTfiles = {}
                    for setting in self.many:   
                        self.ROOTfiles[setting] = get_list(
                            db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                                        where = f"setting=\'{setting}\' AND"+\
                                         f" type=\'{self.simc_type}\' AND"+\
                                         f" kin_study=\'{self.kin_study}\'"))    
                else:    
                    self.ROOTfiles = get_list(
                        db.retrieve(deutDB_name, 'filename', 'SIMC_ROOTfiles', 
                                    where = f"setting=\'{self.many}\' AND"+\
                                     f" type=\'{self.simc_type}\' AND"+\
                                     f" kin_study=\'{self.kin_study}\'"))
                # print(self.ROOTfiles,'\n')
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                               self.trees_sel)               
                self.BRANCHES(self.many,T,self.branches_sel,load)
                
            # Case 3: only setting is given
            # e.g. setting = 'delta_scan_0'
            
            elif setting and simc_type:
                if type(self.many) is list:
                    k = {}
                    self.ROOTfiles = {}
                    for setting in self.many:   
                        self.ROOTfiles[setting], k[setting] =\
                            db.retrieve(deutDB_name, 'filename, kin_study', 'SIMC_ROOTfiles', 
                                        where = f"setting=\'{setting}\' AND"+\
                                         f" type=\'{self.simc_type}\'")[0]    
                else:    
                    self.ROOTfiles, k =\
                            db.retrieve(deutDB_name, 'filename, kin_study', 'SIMC_ROOTfiles', 
                                        where = f"setting=\'{self.many}\' AND"+\
                                         f" type=\'{self.simc_type}\'")[0]    
                self.kin_study = k
                
                T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                               self.trees_sel)
                self.BRANCHES(self.many,T,self.branches_sel,load)

            # Case 4: nothing is given, do you want ALL SIMC files? maybe not..
            
            elif not kin_study and not setting:
                print('======================================================')
                print(' NO kin_study or setting SELECTED')
                print('======================================================')
                print('Run USAGE() for more information')
                print('Do you want to load ALL available SIMC files? Yes/No')
                
                load_all = input('Type "LOAD ALL" for Yes or "NO"/Return for No\n')
                
                if load_all == 'LOAD ALL':
                    print('This might take a while...')
                    
                    self.many = get_list(db.retrieve(
                                deutDB_name, 'setting', 'SIMC_ROOTfiles',
                                distinct=True))
                    self.setting = self.many
                    
                    if type(self.many) is list:                      
                        k = {}
                        self.ROOTfiles = {}                        
                        for setting in self.many:   
                           self.ROOTfiles[setting], k[setting] =\
                               db.retrieve(deutDB_name, 'filename, kin_study',
                                   'SIMC_ROOTfiles', 
                                   where = f"setting=\'{setting}\'")[0]

                    else:
                      self.ROOTfiles, k =\
                             db.retrieve(deutDB_name, 'filename, kin_study',
                                 'SIMC_ROOTfiles', 
                                 where = f"setting=\'{self.many}\'")[0]
                             
                    self.kin_study = k
                  
                    T = self.TREES(self.many,self.ROOTfiles,self.ROOTfiles_path,
                                    self.trees_sel)               
                    self.BRANCHES(self.many,T,self.branches_sel,load)
                    
            # Case 5: no simc_type chosen, this is a must!
            else:
                print('No simc_type given. Must choose a simc_type,\n',
                      'this is written in the column "type" of the\n',
                      'SIMC_ROOTfiles table of deuteron_db')
                
            
        else:
            print('=========================================================')
            print('NO data_type SELECTED')
            print('=========================================================')
            print('data_type must be defined',
                  'The current defined data types are:',
                  '"deut23_data" -> root files from replayed deuteron data.',
                  '\tAvailable root trees: ["T","TSP","TSH","E"]',
                  '\tNote that if you are loading more than 1 tree, the first',
                  '\tone will be selected for analysis.',
                  '"SIMC" -> root files from SIMC data.',sep='\n')
            
    def TREES(self,many,root,path,t):
        self.Trees = {}
        
        try:
            for m in many:
                rootfile = path + root[m]
                if os.path.isfile(rootfile) :
                    self.Trees[m]= R.root_tree(rootfile, trees=t)
                    self.Trees[m].load_trees()
                else:
                    continue
        except TypeError:
                rootfile = path + root
                if os.path.isfile(rootfile) :
                    self.Trees = R.root_tree(rootfile, trees=t)
                    self.Trees.load_trees()
                else:
                    self.Trees = None
                    
        return self.Trees         
    
    def BRANCHES(self,many,t,b,load=False):
        self.Branches = {}
        # print(t)
        if load:
            try:
                temp = {}
                for m in many:
                    tree = t[m]
                    temp[m] = {}
                    for T in tree.trees:
                       tree.select_tree(T)
                       sel = b[T] 
                       temp[m].update(tree.get_branches(sel)) 
            except TypeError:
                tree = t
                temp = {}
                for T in tree.trees:
                   tree.select_tree(T)
                   sel = b[T] 
                   temp.update(tree.get_branches(sel))   
            
            self.Branches = temp           
        
        else:
            self.Branches = None
            
        return self.Branches
    
    # method to define and create boolean arrays (cut arrays)
    # def CUT(self,iterable=self.many,cuts=[],where_branches=self.Branches):
    #     # Inputs
    #     # ----------------------------------
    #     # iterable = list, are we looping over runs? settings? 
    #     #               if empty it will be assumed only 1 is given
    #     # cuts = cut object, see cut_handler for class definition
    #     # where_branches = dict, is a dict containing the cut variable 
        
    #     many = iterable
    #     self.cuts_applied = cuts
    #     br = where_branches
        
        # if many:
            
        # else:
            
        
    
    # def CUT_BRANCHES(self,)
    
    def USAGE(self):
            print('USAGE INFO FOR "DATA_INIT" CLASS:')
            print('This class takes .root trees and leaf data and packages',
                  ' them into "numpy.ndarrays" to prepare them for easy',
                  ' analysis in python.')
            #WIP---------------------------------------------------------------
            print('data_type must be defined.\n\
                  What kind of data are you loading?')
            print('data_type = "deut23_data" \
                  -> root files from deuteron data\n\
                  \t available root trees: ["T","TSP","TSH","E"]\n\
                      \t note if loading more than 1 tree,\
                          the first one will be selected for analysis.\n')

        

