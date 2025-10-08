#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:39:53 2025

@author: gvill

Module to handle cuts from hcana branches
"""
import numpy as np
        
def CUT(iterator=None,cut_dict={},where_branches={}):
    # declare cut variables
    many = iterator
    br = where_branches
    c_dict = cut_dict # where cuts are defined in the following format:
                        # c_dict = {'a':(-1,1),'b':1}
                        # this means I want a cut for variable 'a' 
                        #                                   where -1 < a < 1
                        #  and a cut for variable 'b' where b == 1
    
    if many:
        # initialize variables for later use
        c_var = []
        c_range = tuple()
        c_value = float() or int()
        all_cuts = {}
        
        
        for i in many:
            first = True
            for cv in c_dict:
                c_var = br[i][cv]

                # which type of cut are we doing?
                #   1. range cut e. g. -1 <= a <= 1  
                #   or 2. value cut e. g. a == 1 
                if type(c_dict[cv]) == tuple:
                    c_range = c_dict[cv]

                    c_max = c_var <= c_range[1]
                    c_min = c_var >= c_range[0]

                    this_cut = c_min & c_max
                    this_cut = np.array(this_cut)
                else:
                    c_value = c_dict[cv]
                    
                    this_cut = c_var == c_value
                    this_cut = np.array(this_cut)
                    
                if first:
                    all_cuts[i] = this_cut
                    first = False 
                else:
                    all_cuts[i] = all_cuts[i] & this_cut

    else:        
        # initialize variables for later use
        c_var = []
        c_range = tuple()
        c_value = float() or int()
        all_cuts = []
        first = True
        
        for cv in c_dict:
            c_var = br[cv]
    
            # which type of cut are we doing?
            #   1. range cut e. g. -1 <= a <= 1  
            #   or 2. value cut e. g. a == 1 
            if type(c_dict[cv]) == tuple:
                c_range = c_dict[cv]
    
                c_max = c_var <= c_range[1]
                c_min = c_var >= c_range[0]
    
                this_cut = c_min & c_max
                this_cut = np.array(this_cut)
            else:
                c_value = c_dict[cv]
                
                this_cut = c_var == c_value
                this_cut = np.array(this_cut)
                
            if first:
                all_cuts = this_cut
                first = False 
            else:
                all_cuts = all_cuts & this_cut
            
    return all_cuts