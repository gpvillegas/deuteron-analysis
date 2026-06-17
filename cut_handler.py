#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:39:35 2025

@author: gvill
"""

import numpy as np
import LT.box as B
import database_operations as db
import data_init as D
# import re

#%%
"""
Tests for later:
    -Check class by comparing doing cuts with if statements and the old way.
        should yield same results

"""


"""
WINDOW CUT DEFINTION
-----------------------------
USAGE

An instance of this class will define a window/range 
cut on a range (xmin,xmax) and store simple stats. 

Calling this class will return a boolean array 
that can be used to place a cut on another variable 
or upon itself.

example:
    myvar = np.linspace(0,100)

    cut1 = WCUT(xmin=10,xmax=50,name='10<=x<=50')

    cut1(myvar)
    Out[4]: 
    array([False, False, False, False, False,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True, False, False,
           False, False, False, False, False, False, False, False, False,
           False, False, False, False, False, False, False, False, False,
           False, False, False, False, False])

    cut1.stats()
    10<=x<=50 evaluated 50 times
    40.00% passed
    50.00% were above the cut
    10.00% were below the cut

"""
class WCUT:
    def __init__(self, xmin = -100., xmax = 100., name='window_cut'):
        self.name = name
        self.xmin = xmin
        self.xmax = xmax
        self.init()
        
    def init(self):
        self.n_pass = 0
        self.n_high = 0
        self.n_low = 0
        self.n_reject = 0
        self.tot = 0 
        self.eval = None
        
    def __call__(self, val):
        # deal with nans at the beginning
        self.has_nan = np.isnan(val)     
        if self.has_nan.any():
            val[self.has_nan] = 1e38
            print(f'{val[self.has_nan].size} nan counted')
            
        cut_max = val <= self.xmax
        cut_min = self.xmin <= val        
        self.eval = cut_min & cut_max
        
        #handle nans
        if np.isscalar(val):
            if val == np.nan:
                self.eval = False
            if self.eval:
                self.n_pass = 1
                if val > self.xmax:
                    self.n_high = 1
                else:
                    self.n_low = 1
            return self.eval
    
        else:
           self.eval[self.has_nan] = False
           self.tot = self.eval.size
           self.n_pass = self.eval.sum()
           self.n_high = (~cut_max).sum()
           self.n_low = (~cut_min).sum()
           self.n_reject = (~self.eval).sum()
           return self.eval
    
    def stats(self):
        tot = self.n_pass + self.n_high + self.n_low
        
        if tot != self.tot:
            print('Warning: totals not equal, you might be missing events!')
        
        frac = 0.
        frac_high = 0.
        frac_low = 0.
        frac_drop = 0.
        if tot > 0.:
            frac = self.n_pass/tot*100
            frac_high = self.n_high/tot*100
            frac_low = self.n_low/tot*100
            frac_drop = self.n_reject/tot*100
            print(f'{self.name} evaluated {tot} times')
            print(f'{frac:.2f}% passed: {self.n_pass} events')
            print(f'{frac_high:.2f}% above the cut: {self.n_high} events')
            print(f'{frac_low:.2f}% below the cut: {self.n_low} events')
            print(f'{frac_drop:.2f}% rejected: {self.n_reject} events\n')
        
        return 

    
    def __repr__(self):
        return \
            f"Window cut: {self.name}; xmin = {self.xmin}, xmax = {self.xmax}\n"
    
    def report(self, fname=''):
        with open(fname,'a') as f:
            f.write(f"Window cut: {self.name}; xmin = {self.xmin}, xmax = {self.xmax}\n")
        return
"""
VALUE CUT DEFINITION
-----------------------------
USAGE
An instance of this class will define a value/number 
cut which will pass numbers equal to the value a
nd store simple stats. 

Calling this class will return a boolean array 
that can be used to place a cut on another variable 
or upon itself.
"""

class VCUT:
    def __init__(self, value=100., name='value_cut'):
        self.name = name
        self.value = value
        self.init()
        
    def init(self):
        self.n_pass = 0
        self.n_reject = 0
        self.eval = None
        
    def __call__(self, val):
        self.eval = val == self.value
        #handle nans
        if np.isscalar(val):
            if val == np.nan:
                self.eval = False
                if self.eval:
                    self.n_pass = 1
                else:
                    self.n_reject = 1
            return self.eval
    
        else:       
            self.has_nan = np.isnan(val) 
            self.eval[self.has_nan] = False
            self.tot = self.eval.size
            self.n_pass = self.eval.sum()
            self.n_reject = (~self.eval).sum()
            return self.eval
    
    def stats(self):
        tot = self.n_pass + self.n_reject
        # is tot same as self.tot?
        if tot != self.tot:
            print('Warning: totals not equal\nyou might be missing events!')
        
        frac = 0.
        if tot > 0.:
            frac = self.n_pass/tot*100
            frac_drop = self.n_reject/tot*100
        print(f'{self.name} evaluated {tot} times')
        print(f'{frac:.2f}% passed')
        print(f'{frac_drop:.2f}% rejected')

    
    def __repr__(self):
        return f"Value cut: {self.name}; value = {self.value}\n"
    
    def report(self, fname=''):
        with open(fname,'a') as f:
            f.write(f"Value cut: {self.name}; value = {self.value}\n")
        return
#%% class for collimator cut
## tested - works well
# on run 20851 74.3% of events pass hms coll_cut and 80.75% shms coll_cut
# plots of y_coll:x_coll recreate the collimator pattern well
# -- added many functionality to work with multiple runs

# database name
deutDB_name= 'deuteron_db.db'

dtr = np.pi/180.

# To get rid of unhelpful groupings (lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0]

class coll_cut:
    def __init__(self,data_obj,spec='HMS',many=False,is_SIMC=False):
        if many:
            self.many = data_obj.many
        else:
            self.many = False
            
        self.spec = spec
        
        #HMS collimator limits
        self.hcoll_lim = np.array([(-2.2875,11.646),(2.2875,11.646),
                                   (4.575,5.823),(4.575,-5.823),
                                   (2.2875,-11.646),(-2.2875,-11.646),
                                   (-4.575,-5.823),(-4.575,5.823),
                                   (-2.2875,11.646)])
        #SHMS collimator limits
        self.pcoll_lim = np.array([(-4.25,12.5),(4.25,12.5),(8.5,6.25),
                                   (8.5,-6.25),(4.25,-12.5),(-4.25,-12.5),
                                   (-8.5,-6.25),(-8.5,6.25),(-4.25,12.5)])
        
        if self.spec == 'HMS':
            if is_SIMC:
                try:
                    if type(data_obj.many) == list:
                        self.xcoll = {}
                        self.ycoll = {}
                        self.n = {}
                        for m in data_obj.many:
                            tar_x = data_obj.Branches[m]['tar_x']
                            h_xptar = data_obj.Branches[m]['h_xptar']
                            h_yptar = data_obj.Branches[m]['h_yptar']
                            htar_z = data_obj.Branches[m]['h_zv']
                            h_ytar = data_obj.Branches[m]['h_ytar']
                            
                            h_th = get_list(
                                db.retrieve(deutDB_name, 'HMS_Angle', 'RUN_LIST_UPDATED', 
                                    where= f"kin_study=\'{data_obj.kin_study}\'"+\
                                            f" AND setting=\'{m}\'",
                                            distinct=True))
                        
                            # calculate corr tarx
                            htarx_corr = tar_x - h_xptar*htar_z*np.cos(h_th*dtr)
    
                            	# Define Collimator (same as in HCANA)
                            hXColl = htarx_corr + h_xptar*168.   #in cm
                            hYColl = h_ytar + h_yptar*168.
                        
                            # hXColl = data_obj.Branches['xcoll']   #in cm
                            # hYColl = data_obj.Branches['ycoll']
                            n = hXColl.size
                            
                            self.xcoll[m] = hXColl
                            self.ycoll[m] = hYColl
                            self.n[m] = n
                            
                        self.nlim = self.hcoll_lim[:,0].size
                        self.xlim = self.hcoll_lim[:,0]
                        self.ylim = self.hcoll_lim[:,1]
                    else:
                        raise(TypeError)
                        
                except TypeError:
                    tar_x = data_obj.Branches['tar_x']
                    h_xptar = data_obj.Branches['h_xptar']
                    h_yptar = data_obj.Branches['h_yptar']
                    htar_z = data_obj.Branches['h_zv']
                    h_ytar = data_obj.Branches['h_ytar']
                    
                    h_th = get_list(
                        db.retrieve(deutDB_name, 'HMS_Angle', 'RUN_LIST_UPDATED', 
                            where= f"kin_study=\'{data_obj.kin_study}\'"+\
                                    f" AND setting=\'{data_obj.setting}\'",
                                    distinct=True))
                    
                    # calculate corr tarx
                    htarx_corr = tar_x - h_xptar*htar_z*np.cos(h_th*dtr)
    
                    	# Define Collimator (same as in HCANA)
                    hXColl = htarx_corr + h_xptar*168.   #in cm
                    hYColl = h_ytar + h_yptar*168.
                    
                    # hXColl = data_obj.Branches['xcoll']   #in cm
                    # hYColl = data_obj.Branches['ycoll']
                    n = hXColl.size
                    
                    self.xcoll = hXColl
                    self.ycoll = hYColl
                    self.n = n
                    self.nlim = self.hcoll_lim[:,0].size
                    self.xlim = self.hcoll_lim[:,0]
                    self.ylim = self.hcoll_lim[:,1]
	  
            else:    
                if many:
                    self.xcoll = {}
                    self.ycoll = {}
                    self.n = {}
                    for m in self.many:
                        xc = data_obj.Branches[m]['H.extcor.xsieve']
                        yc = data_obj.Branches[m]['H.extcor.ysieve']
                        nc = xc.size
                        
                        self.xcoll[m] = xc
                        self.ycoll[m] = yc
                        self.n[m] = nc
                        
                    self.nlim = self.hcoll_lim[:,0].size
                    self.xlim = self.hcoll_lim[:,0]
                    self.ylim = self.hcoll_lim[:,1]
                else:    
                    self.xcoll = data_obj.Branches['H.extcor.xsieve']
                    self.ycoll = data_obj.Branches['H.extcor.ysieve']
                    self.n = self.xcoll.size
                    self.nlim = self.hcoll_lim[:,0].size
                    self.xlim = self.hcoll_lim[:,0]
                    self.ylim = self.hcoll_lim[:,1]
                
            
        if self.spec == 'SHMS':
            if is_SIMC:
                tar_x = data_obj.Branches['tar_x']
                e_xptar = data_obj.Branches['e_xptar']
                e_yptar = data_obj.Branches['e_yptar']
                e_delta = data_obj.Branches['e_delta']
                etar_z = data_obj.Branches['e_zv']
                e_ytar = data_obj.Branches['e_ytar']
                
                e_th = get_list(
                    db.retrieve(deutDB_name, 'SHMS_Angle', 'RUN_LIST_UPDATED', 
                        where= f"kin_study=\'{data_obj.kin_study}\'"+\
                                f" AND setting=\'{data_obj.setting}\'",
                                distinct=True))
                
                # calculate corrected tar x    
                etarx_corr = tar_x - e_xptar*etar_z*np.cos(e_th*dtr)
                
                # calculate XColl YColl as defined in HCANA in cm
                eXColl = etarx_corr + e_xptar*253.
                
                # correct for HB horizontal bend
                eYColl = e_ytar + e_yptar*253.-\
                    (0.019+40.*.01*0.052)*e_delta+\
                        (0.00019+40*.01*.00052)*e_delta*e_delta
                n = eXColl.size
                
                self.xcoll = eXColl
                self.ycoll = eYColl
                self.n = n
                self.nlim = self.pcoll_lim[:,0].size
                self.xlim = self.pcoll_lim[:,0]
                self.ylim = self.pcoll_lim[:,1]        
                
            else:    
                if many:
                    self.xcoll = {}
                    self.ycoll = {}
                    self.n = {}
                    for m in self.many:
                        xc = data_obj.Branches[m]['P.extcor.xsieve']
                        yc = data_obj.Branches[m]['P.extcor.ysieve']
                        nc = xc.size
                        
                        self.xcoll[m] = xc
                        self.ycoll[m] = yc
                        self.n[m] = nc
                        
                    self.nlim = self.hcoll_lim[:,0].size
                    self.xlim = self.hcoll_lim[:,0]
                    self.ylim = self.hcoll_lim[:,1]
                else:
                    self.xcoll = data_obj.Branches['P.extcor.xsieve']
                    self.ycoll = data_obj.Branches['P.extcor.ysieve']
                    self.n = self.xcoll.size
                    self.nlim = self.pcoll_lim[:,0].size
                    self.xlim = self.pcoll_lim[:,0]
                    self.ylim = self.pcoll_lim[:,1]
        
    def is_inside(self,n, x, y, nlim, xlim, ylim):
        is_point_inside = []
        for i in range(n):
            ncross = 0
            for j in range(nlim-1):
                if (ylim[j] < y[i]) and (ylim[j+1] < y[i]): continue
                if (x[i] == xlim[j]): continue        
                t = x[i] - xlim[j]
                s = xlim[j+1] - x[i]
                if (t*s < 0.): continue       
                di = (ylim[j+1] - ylim[j])/(xlim[j+1] - xlim[j])
                f = ylim[j] + di*(x[i] - xlim[j])
                if (f < y[i]): continue       
                ncross += 1   
            is_point_inside.append((ncross % 2) == 1)    
        return np.array(is_point_inside)
    
    # collimator cut
    def __call__(self):
        if self.many:
            self.res = {}
            for m in self.many:
                self.res[m] = self.is_inside(self.n[m],
                                             self.ycoll[m],self.xcoll[m],
                                             self.nlim,self.xlim,self.ylim)
        else:        
            self.res = self.is_inside(self.n,self.ycoll,self.xcoll,
                                      self.nlim,self.xlim,self.ylim)
        return self.res
    
    def stats(self):
        if self.many:
            for m in self.many:
                passed = self.res[m].sum()/self.res[m].size*100
                print(f'Run {m} {self.spec} collimator cut stats:')
                print(f'{passed:.2f}% passed\n')
        else:        
            passed = self.res.sum()/self.res.size*100
            print(f'{self.spec} collimator cut stats:')
            print(f'{passed:.2f}% passed\n')

    def __repr__(self):
        return f"{self.spec} collimator cut ON\n"
    
    def report(self, fname=''):
        with open(fname,'a') as f:
            f.write(f"{self.spec} collimator cut ON\n")
        return
#%% class for current cut

##checked for interval accuracy by manually determining good/bad intervals
# and checking for their inclusion in .intervals after cut.
# -- added many functionality to handle multiple runs
# -- added cut_lim options for tighter cuts, options are: 'std' --> >5muA and
#       'tgt_boil' --> cuts 3muA around max curr

class current_cut:
    def __init__(self, data_obj, current='BCM4A',many=False,cut_lim = 'std'):
        if many:
            self.many = data_obj.many
            
            self.evNumber = {}
            self.current = {}
            self.c_cut_min = {}
            self.c_cut_max = {}
            for m in self.many:
                evnum = data_obj.Branches[m]['evNumber']
                curr = data_obj.Branches[m][f'P.{current}.scalerCurrent']
                
                # cut_lim options: std, tgt_boil
                if cut_lim == 'tgt_boil':
                    # read avg current
                    # avg_current = get_list(db.retrieve('deuteron_db.db', 
                    #                           f'{current}_current', 
                    #                           'RUN_LIST_UPDATED', 
                    #                           where = f"run=\'{m}\'"))
                    hc = B.histo(curr,range=(0,100),bins=100)
                    max_loc = np.where(hc.bin_content == max(hc.bin_content))[0]
                    max_c = hc.bin_center[max_loc]  

                    self.c_cut_min[m] = max_c[0] - 5.
                    self.c_cut_max[m] = max_c[0] + 5.
                    
                elif cut_lim == 'std':
                    self.c_cut_min[m] = 5.
                    self.c_cut_max[m] = np.inf
                else:
                    print('Please specify current cut limits.')
                    print('Options: "std", "tgt_boil"')
                    
                self.evNumber[m] = evnum 
                self.current[m] = curr
                
                
        else:
            self.many = False
            self.evNumber = data_obj.Branches['evNumber'] 
            self.current = data_obj.Branches[f'P.{current}.scalerCurrent']
            
            if cut_lim == 'tgt_boil':
                hc = B.histo(self.current,range=(0,100),bins=100)
                max_loc = np.where(hc.bin_content == max(hc.bin_content))[0]
                max_c = hc.bin_center[max_loc]  

                self.c_cut_min = max_c[0] - 5.
                self.c_cut_max = max_c[0] + 5.
                
            elif cut_lim == 'std':
                self.c_cut_min = 5.
                self.c_cut_max = np.inf
            else:
                print('Please specify current cut limits.')
                print('Options: "std", "tgt_boil"')
        
        self.cut_name = f'{current}_cut'
  
    def __call__(self):
        if self.many:
            self.intervals = {}
            self.res = {}
            self.cut_array = {}
            for m in self.many:
                curr_cut = WCUT(self.c_cut_min[m],self.c_cut_max[m],
                                name=f'{self.cut_name}')
                self.cut_array[m] = curr_cut(self.current[m])
                cut_array = curr_cut(self.current[m])
                cut_array = np.roll(cut_array,-1)[:-1]
                intervals = np.vstack((self.evNumber[m],
                                       np.roll(self.evNumber[m],-1)))
                
                inter = intervals.T[:-1][cut_array]
                
                n = self.evNumber[m][-1] - 1 
                                            
                r = np.full(int(n),0,dtype=bool)
                
                for i in inter:
                    for j in range(int(i[0])-1,int(i[1])):
                        if j < n:
                            r[j] = True
                
                self.intervals[m] = inter
                self.res[m] = r
            
        else:    
            curr_cut = WCUT(self.c_cut_min,self.c_cut_max,name=f'{self.cut_name}')
            
            # we ignore the first element in the cut so 
            # it can be applied to ranges
            # in the ranges it is the next number that
            # determines if it is a good range or not
            self.cut_array = curr_cut(self.current)
            
            cut_array = curr_cut(self.current)
            cut_array = np.roll(cut_array,-1)[:-1] 
    
            # create array of intervals using the event numbers and 
            # the same array rolled over by 1
            
            intervals = np.vstack((self.evNumber,np.roll(self.evNumber,-1)))
            # then take the transpose to have an array where each line
            # is a range, the last range is [last item, first item]
            # so its not needed.
            self.intervals = intervals.T[:-1] 
            
            # apply cut to the intervals to get only good intervals             
            self.intervals = self.intervals[cut_array]
            
            n = self.evNumber[-1] - 1 #evNumber last element is total elements in 
                                        # phys branch + 1
                                        
            self.res = np.full(int(n),0,dtype=bool)
            
            for i in self.intervals:
                for j in range(int(i[0])-1,int(i[1])):
                    if j < n:
                        self.res[j] = True  
                        
        return self.res
    
    def scaler(self,scaler_var,cut_array):
        # this method has to be called after calling the current cut by itself, 
        # so that the boolean cut arrays are defined.
        # size of scaler array
        n = scaler_var.size
        # array positions in scaler array
        pos = np.arange(0,n,1)
        # intervals made from boolean array, intervals [True,False] or 
        # [False,True] or [False, False] show ranges where scaler counts should
        # not be counted 
        
        inter = np.vstack((cut_array,
                           np.roll(cut_array,-1)))
        intervals = inter.T[:-1] # do not need last one [Last, First]
        
        j=0
        pos_intervals = []
        t = [-1,-1]
        for i in intervals:
            if i[0] == True and i[1] == False or i[0] == False and i[1] == True\
                or j==0 and i[0] == False:
                # get the positions in the scaler array that correspond
                # to [False,True] or [True,False] (the start and end of a
                # bad current period)
                if i[0] == True:
                    t[0] = int(pos[j])
                elif j==0 and i[0] == False:
                    t[0] = int(pos[j])
                elif j != 0 and i[0] == False:
                    t[1] = int(pos[j]) +1
            # here we are basically transforming the [True,False], [False,True]
            # intervals into intervals of their positions, saving them in
            # 'pos_inter'
            if t[0]!= -1 and t[1] != -1:
                pos_intervals.append(t)
                t = [-1,-1]
            j+=1
        
        # loop over all 'pos_intervals' and add up all scaler events counted
        # during those intervals, these are the bad scaler events to be 
        # subtracted
        bad_scaler_var = 0
        for i in pos_intervals:
            try:
                bad_scaler_var += scaler_var[i[1]] - scaler_var[i[0]]
            except IndexError:
                bad_scaler_var += scaler_var[i[1]-1] - scaler_var[i[0]] 
        
        newScalerVarTot = scaler_var[-1] - bad_scaler_var 
        
        return newScalerVarTot
    
    def stats(self):
        if self.many:
            for m in self.many:
                passed = self.res[m].sum()/self.res[m].size*100
                print(f'Run {m} {self.cut_name} stats:')
                print(f'{passed:.2f}% passed\n')
        else:        
            passed = self.res.sum()/self.res.size*100
            print(f'{self.cut_name} stats:')
            print(f'{passed:.2f}% passed\n')

    def __repr__(self):
            return f"{self.cut_name}\n"
    
    def report(self, fname=''):
        with open(fname,'a') as f:
            f.write(f"{self.cut_name}\n")
        return
#%% z_tar difference cut
# cuts around max zdiff
# tested 
class ztar_cut:
    def __init__(self,data_obj,many=False,is_SIMC=False):
        if is_SIMC:
            if many:
                self.many = data_obj.many
                self.hms_ztar = {}
                self.shms_ztar = {}
                for m in self.many:
                    self.hms_ztar[m] = data_obj.Branches[m]['h_zv']
                    self.shms_ztar[m] = data_obj.Branches[m]['e_zv']
            else:
                self.many = False
                self.hms_ztar = data_obj.Branches['h_zv']
                self.shms_ztar = data_obj.Branches['e_zv']
        else:    
            if many:
                self.many = data_obj.many
                self.hms_ztar = {}
                self.shms_ztar = {}
                for m in self.many:
                    self.hms_ztar[m] = data_obj.Branches[m]['H.react.z']
                    self.shms_ztar[m] = data_obj.Branches[m]['P.react.z']
            else:
                self.many = False
                self.hms_ztar = data_obj.Branches['H.react.z']
                self.shms_ztar = data_obj.Branches['P.react.z']
        
    def __call__(self):
        if self.many:
            self.cut = {}
            self.res = {}
            self.cut_min = {}
            self.cut_max = {}
            
            for m in self.many:
                ztar_diff = self.hms_ztar[m] - self.shms_ztar[m]
                
                hc = B.histo(ztar_diff,range=(-5.,5.),bins=50)
                hc.fit(plot_fit=False) ## maybe modify LTbox to not print fit results
                # max_loc = np.where(hc.bin_content == max(hc.bin_content))[0]
                # max_c = hc.bin_center[max_loc]  
                
                cmin = hc.mean.value - 2.0
                cmax = hc.mean.value + 2.0
                
                cut = WCUT(cmin,cmax,name='ztar_diff_cut')
                r = cut(ztar_diff)
                
                self.cut[m] = cut
                self.res[m] = r
                self.cut_min[m] = cmin
                self.cut_max[m] = cmax
        else:
            ztar_diff = self.hms_ztar - self.shms_ztar
            
            hc = B.histo(ztar_diff,range=(-5.,5.),bins=50)
            hc.fit(plot_fit=False)
            # max_loc = np.where(hc.bin_content == max(hc.bin_content))[0]
            # max_c = hc.bin_center[max_loc]  
            
            self.cut_min = hc.mean.value - 2.0
            self.cut_max = hc.mean.value + 2.0
            self.cut = WCUT(self.cut_min,self.cut_max,name='ztar_diff_cut')
            self.res = self.cut(ztar_diff)
        
        return self.res
    
    def stats(self):
        if self.many:
            for m in self.many:
                self.cut[m].stats()
        else:    
            self.cut.stats()
        
    def __repr__(self):
            return "ztar_diff_cut\n"
    
    def report(self, fname=''):
        with open(fname,'a') as f:
            f.write("ztar_diff_cut\n")
        return   
#%% coincidence time cut - corrects for offset

class CTime_cut:
    def __init__(self,data_obj,many=False):
        if many:
            self.many = data_obj.many
            
            self.coinTime = {}
            for m in self.many:
                self.coinTime[m] = data_obj.Branches[m]['CTime.epCoinTime_ROC2']
        else:    
            self.many = many
            self.coinTime = data_obj.Branches['CTime.epCoinTime_ROC2']

    def get_coinTimeOffset(self,coinTime):    
        CTime = coinTime  
        CTime_histo_0 = B.histo(CTime,range=(0,1000),bins=1000,
                                  title = 'CTime no cuts')
        # CTime_histo_0.fit(plot_fit=False,print_res=False)
        # mu_0 = CTime_histo_0.mean.value  
        # sigma_0 = abs(CTime_histo_0.sigma.value)
        # print(mu_0,sigma_0)
        max_loc = np.where(CTime_histo_0.bin_content ==\
                           max(CTime_histo_0.bin_content))[0]
        max_c = CTime_histo_0.bin_center[max_loc][0]

        x_min = max_c - 10.
        x_max = max_c + 10.
        
        CTime_histo = B.histo(CTime,range=(x_min,x_max),bins=200,
                                  title = 'CTime no cuts')    
        CTime_histo.fit(max_c-5,max_c+5,plot_fit=False)
        mu = CTime_histo.mean.value    
        CTime_corr = CTime - mu
        
        # B.pl.figure()
        # CTime_histo_0.plot()
        # B.pl.figure()
        # CTime_histo.plot()
        # B.pl.vlines(mu,ymin=0.,ymax=CTime_histo.A.value,color='black',linestyles='--')
        # CTime_histo.plot_fit()
        
        return CTime_corr

    def __call__(self):
        if self.many:
            self.cut = {}
            self.res = {}
            self.CTime_corr = {}
            for m in self.many:
                coinTime_corr = self.get_coinTimeOffset(self.coinTime[m]) 
                
                cut = WCUT(-2.,2.,name='epCoinTime_cut')
                res = cut(coinTime_corr)
                
                self.cut[m] = cut
                self.res[m] = res
                self.CTime_corr[m] = coinTime_corr
                
                # h = B.histo(coinTime_corr,range=(-5,5),bins=200)
                # h1 = B.histo(coinTime_corr[res],range=(-5,5),bins=200)
                # B.pl.figure()
                # h.plot()
                # h1.plot()
                
        else:
            coinTime_corr = self.get_coinTimeOffset(self.coinTime)
            
            self.cut = WCUT(-2.,2.,name='epCoinTime_cut')
            self.res = self.cut(coinTime_corr)
            self.CTime_corr = coinTime_corr
        
        return self.res
    
    def stats(self):
        if self.many:
            for m in self.many:
                self.cut[m].stats()
        else:    
            self.cut.stats()
        
        
    def __repr__(self):
            return 'epCoinTime_cut\n'
    
    def report(self, fname=''):
        with open(fname,'a') as f:
            f.write('epCoinTime_cut\n')
        return
#%% Commonly used cuts are defined here

"""
STANDARD CUTS
std_cuts = {'H.gtr.dp':(-10.,10.),'P.gtr.dp':(-10.,22.),
            'H.gtr.th':(-0.08,0.08),'H.gtr.ph':(-0.03,0.03),
            'P.gtr.th':(-0.03,0.03),'P.gtr.ph':(-0.03,0.03)}
"""
HCANA_names = {'hms_delta':'H.gtr.dp','shms_delta':'P.gtr.dp',
               'hms_xptar':'H.gtr.th','hms_yptar':'H.gtr.ph',
               'shms_xptar':'P.gtr.th','shms_yptar':'P.gtr.ph',
               'shms_calPID':'P.cal.etottracknorm',
               'Em_cut':'H.kin.secondary.emiss_nuc',
               'epCoinTime':'CTime.epCoinTime_ROC2_corr',
               'shms_calPID':'P.cal.etottracknorm',
               'Q2_cut':'P.kin.primary.Q2',
               'xbj_cut':'P.kin.primary.x_bj',
               'Em_cut_hc':'H.kin.secondary.emiss',
               'noEDTM':'T.coin.pEDTM_tdcTimeRaw',
               'W_cut':'P.kin.primary.W',
               'xbj_cut':'P.kin.primary.x_bj'}

SIMC_names = {'hms_delta':'h_delta','shms_delta':'e_delta',
               'hms_xptar':'h_xptar','hms_yptar':'h_yptar',
               'shms_xptar':'e_xptar','shms_yptar':'e_yptar',
               'Em_cut':'Em','Q2_cut':'Q2','Em_cut_hc':'Em',
               'W_cut':'W'}

#no edtm cut
noEDTM = VCUT(0.,name='noEDTM')

#whole acceptance
hms_delta = WCUT(-8.,8.,name='hms_delta')
shms_delta = WCUT(-10.,22.,name='shms_delta')

#tight acceptance
# hms_delta = WCUT(-5.,5.,name='hms_delta')
# shms_delta = WCUT(-5.,5.,name='shms_delta')

## whole acceptance
hms_xptar = WCUT(-0.1,0.1,name='hms_xptar')
hms_yptar = WCUT(-0.04,0.04,name='hms_yptar')
shms_xptar = WCUT(-0.05,0.05,name='shms_xptar')
shms_yptar = WCUT(-0.03,0.03,name='shms_yptar')

## half acceptance
# hms_xptar = WCUT(-0.05,0.05,name='hms_xptar')
# hms_yptar = WCUT(-0.02,0.02,name='hms_yptar')
# shms_xptar = WCUT(-0.025,0.025,name='shms_xptar')
# shms_yptar = WCUT(-0.015,0.015,name='shms_yptar')

## quarter acceptance
# hms_xptar = WCUT(-0.025,0.025,name='hms_xptar')
# hms_yptar = WCUT(-0.01,0.01,name='hms_yptar')
# shms_xptar = WCUT(-0.0125,0.0125,name='shms_xptar')
# shms_yptar = WCUT(-0.0075,0.0075,name='shms_yptar')

shms_calPID = WCUT(0.7,1.3,name='shms_calPID')

acceptance_cuts = [hms_delta,shms_delta,hms_xptar,hms_yptar,shms_xptar,shms_yptar] 
    
## Event Selection Cuts
# Coinidence Time cut
# ct_cut = WCUT(xmin=-2.,xmax=2.,name='epCoinTime')
# Missing Energy
Em_cut = WCUT(-0.05,0.05,name='Em_cut')
Em_cut_hc = WCUT(-0.01,0.01,name='Em_cut_hc')
# z-target difference
# Ztar_diff = WCUT()
# Momentum Transfer Squared
Q2_cut = WCUT(4.0,np.inf,name='Q2_cut')
# xBj cut for heep
xbj_cut = WCUT(1.0,np.inf,name='xbj_cut')
# invariant mass (W) cut for heep
W_cut = WCUT(0.85,1.05,name='W_cut')

event_selection_cuts = [hms_delta,shms_delta,shms_calPID,Em_cut,Q2_cut]

heep_event_selection_cuts = [hms_delta,shms_delta,Em_cut_hc,shms_calPID]

event_selection_SIMC = [hms_delta,shms_delta,Em_cut,Q2_cut]

heep_event_selection_SIMC = [hms_delta,shms_delta,Em_cut_hc]

#%% 
rtd = 180/np.pi

###
# made a function that applies the event selection cuts and stores the variables
# with cuts back into the orginal data object
###

def apply_evt_sel_cuts(data_obj,is_SIMC=False,show_stats = False,save_report=False,
                       STD=True, HCOLL=True, PCOLL=True, CURR=True, ZT=True, CT=True,
                       calc_avg_kin=False):
    if is_SIMC:
       
        cuts_applied_list = []
        
        if STD:
            cuts_list = [hms_delta,shms_delta,Em_cut,Q2_cut]
            cuts_applied_list = cuts_list.copy()
        
        if HCOLL:
            hcoll_cut = coll_cut(data_obj,spec='HMS',is_SIMC=True)
            hcoll_cut_arrays = hcoll_cut()
            cuts_applied_list.append(hcoll_cut)
        if PCOLL:
            pcoll_cut = coll_cut(data_obj,spec='SHMS',is_SIMC=True)
            pcoll_cut_arrays = pcoll_cut()
            cuts_applied_list.append(pcoll_cut)
        if ZT:
            zt_cut = ztar_cut(data_obj,is_SIMC=True)
            zt_cut_arrays = zt_cut()
            cuts_applied_list.append(zt_cut)
        
        if STD:
            cuts_to_apply_list = []
            for cut in cuts_list:
                br = data_obj.Branches[SIMC_names[cut.name]]
                cut_array = cut(br)
                # cut.stats()
                
                cuts_to_apply_list.append(cut_array)
                
        # add special cut arrays to clist
        if HCOLL: cuts_to_apply_list.append(hcoll_cut_arrays)
        if PCOLL: cuts_to_apply_list.append(pcoll_cut_arrays)
        if ZT: cuts_to_apply_list.append(zt_cut_arrays)
        
        all_cuts_sim = cuts_to_apply_list[0]
        for arr in cuts_to_apply_list:    
            all_cuts_sim = all_cuts_sim & arr 
            
        Pm_cut = data_obj.Branches['Pm'][all_cuts_sim]
        Th_rq_cut = data_obj.Branches['theta_rq'][all_cuts_sim]*rtd
        em_cut = data_obj.Branches['Em'][all_cuts_sim]
        WEIGHTS_cut = D.calc_weights(data_obj.Branches)[all_cuts_sim]
        WEIGHTS_PS_cut = D.calc_weights_PS(data_obj.Branches)[all_cuts_sim]        

        data_obj.Branches.update({'Pm_cut':Pm_cut,
                                       'theta_rq_cut':Th_rq_cut,
                                       'Em_cut':em_cut,
                                       'WEIGHTS_cut':WEIGHTS_cut,
                                       'WEIGHTS_PS_cut':WEIGHTS_PS_cut})
        print('SIMC data object updated with the branches:',
              'Pm_cut, theta_rq_cut, Em_cut\n',
              'WEIGHTS_cut, WEIGHTS_PS_cut\n')
        print('Applied the following cuts:\n',
              [c for c in cuts_applied_list])
        
        if calc_avg_kin:
            Ein_v_cut = data_obj.Branches['Ein_v_GeV'][all_cuts_sim]
            kf_v_cut = data_obj.Branches['kf_v'][all_cuts_sim]
            the_v_cut = data_obj.Branches['the_v'][all_cuts_sim]
            Pf_v_cut = data_obj.Branches['pf_v_GeV'][all_cuts_sim]
            thp_v_cut = data_obj.Branches['thp_v'][all_cuts_sim]
            q_v_cut = data_obj.Branches['q_lab_v_GeV'][all_cuts_sim]
            q_v_calc_cut = data_obj.Branches['q_v_calc'][all_cuts_sim]
            thq_v_cut = data_obj.Branches['thq_v'][all_cuts_sim]
            Q2_v_cut = data_obj.Branches['Q2_v_GeV2'][all_cuts_sim]
            omega_v_cut = data_obj.Branches['nu_v_GeV'][all_cuts_sim]
            xbj_v_cut = data_obj.Branches['xbj_v'][all_cuts_sim]
            pm_v_cut = data_obj.Branches['pm_v_GeV'][all_cuts_sim]
            pm_v_calc_cut = data_obj.Branches['pm_v_calc'][all_cuts_sim]
            thpq_v_cut = data_obj.Branches['th_pq_v'][all_cuts_sim]
            thnq_v_cut = data_obj.Branches['th_nq_v'][all_cuts_sim]
            phpq_v_cut = data_obj.Branches['ph_pq_v'][all_cuts_sim]
            phnq_v_cut = data_obj.Branches['ph_nq_v'][all_cuts_sim]
            sig_cut = data_obj.Branches['sig'][all_cuts_sim]
            
            kin = {'Ein_v_cut':Ein_v_cut,
                   'kf_v_cut':kf_v_cut,
                   'the_v_cut':the_v_cut,
                   'pf_v_cut':Pf_v_cut,
                   'thp_v_cut':thp_v_cut,
                   'q_v_cut':q_v_cut,
                   'q_v_calc_cut':q_v_calc_cut,
                   'thq_v_cut':thq_v_cut,
                   'Q2_v_cut':Q2_v_cut,
                   'omega_v_cut':omega_v_cut,
                   'xbj_v_cut':xbj_v_cut,
                   'pm_v_cut':pm_v_cut,
                   'pm_v_calc_cut':pm_v_calc_cut,
                   'th_pq_v_cut':thpq_v_cut,
                   'th_nq_v_cut':thnq_v_cut,
                   'ph_pq_v_cut':phpq_v_cut,
                   'ph_nq_v_cut':phnq_v_cut,
                   'sig_cut':sig_cut
                    }
            
            data_obj.Branches.update(kin)
            print('SIMC data object updated with the branches:',
                      kin.keys(),'\n')
        
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
            Pm = {}
            Th_rq = {}
            Em = {}

            for m in data_obj.many:
                Pm[m] = data_obj.Branches[m]['H.kin.secondary.pmiss']
                Th_rq[m] = data_obj.Branches[m]['H.kin.secondary.th_bq']*rtd
                Em[m] = data_obj.Branches[m]['H.kin.secondary.emiss_nuc']    

            ## making the cut lists
            cuts_applied_list = []
            # first select the cuts to apply, defined already in cut_handler
            
            if STD:
                cuts_list = [hms_delta, shms_delta, shms_calPID, 
                             Em_cut, Q2_cut, noEDTM]
                
                cuts_applied_list = cuts_list.copy()
            
            # initialize special cut classes: coll_cut, current_cut, ztar_cut,
            # CTime_cut 
            if HCOLL:
                hcoll_cut = coll_cut(data_obj,spec='HMS',many=True)
                hcoll_cut_arrays = hcoll_cut()
                cuts_applied_list.append(hcoll_cut)
            
            if PCOLL:
                pcoll_cut = coll_cut(data_obj,spec='SHMS',many=True)
                pcoll_cut_arrays = pcoll_cut()
                cuts_applied_list.append(pcoll_cut)
            
            if CURR:
                curr_cut = current_cut(data_obj,current='BCM1',many=True)
                curr_cut_arrays = curr_cut()
                cuts_applied_list.append(curr_cut)
                       
            if ZT:
                zt_cut = ztar_cut(data_obj,many=True)
                zt_cut_arrays = zt_cut()
                cuts_applied_list.append(zt_cut)
            
            if CT:
                ct_cut = CTime_cut(data_obj,many=True)
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
                
                if STD:
                    clist = []
                    for cut in cuts_list:
                        br = data_obj.Branches[m][HCANA_names[cut.name]]
                        cut_array = cut(br)
                        # cut.stats()
                        
                        clist.append(cut_array)
                
                # add special cut arrays to clist
                if HCOLL: {clist.append(hcoll_cut_arrays[m])}
                if PCOLL: {clist.append(pcoll_cut_arrays[m])}
                if CURR: {clist.append(curr_cut_arrays[m])}
                if ZT: {clist.append(zt_cut_arrays[m])}
                if CT: {clist.append(ct_cut_arrays[m])}
                
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
                em_cut = Em[m][all_cuts[m]]
                
                data_obj.Branches[m].update({'H.kin.secondary.pmiss_cut':Pm_cut,
                                               'H.kin.secondary.th_bq_cut':Th_rq_cut,
                                               'H.kin.secondary.emiss_nuc_cut':em_cut})
            
            print('Data object updated with the branches:\n',
                  'H.kin.secondary.pmiss_cut, H.kin.secondary.th_bq_cut,',
                  'H.kin.secondary.emiss_nuc_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_applied_list])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_applied_list])
            if save_report:
                filename = 'yield_n_xsec/current_cuts_DATA.txt'
                with open(filename,'w') as f:
                    f.write('CUTS REPORT\n')
                a = [c.report(filename) for c in cuts_applied_list] 
                        
        except TypeError:
            
            Pm = data_obj.Branches['H.kin.secondary.pmiss']
            Th_rq = data_obj.Branches['H.kin.secondary.th_bq']*rtd
            Em = data_obj.Branches['H.kin.secondary.emiss_nuc']    

            ## making the cut lists
            cuts_applied_list = []
            # first select the cuts to apply, defined already in cut_handler
            if STD:
                cuts_list = [hms_delta, shms_delta, shms_calPID, 
                             Em_cut, Q2_cut, noEDTM]
                
                cuts_applied_list = cuts_list.copy()
            
            # initialize special cut classes: coll_cut, current_cut, ztar_cut,
            # CTime_cut 
            if HCOLL:
                hcoll_cut = coll_cut(data_obj,spec='HMS')
                hcoll_cut_arrays = hcoll_cut()
                cuts_applied_list.append(hcoll_cut)
            if PCOLL:
                pcoll_cut = coll_cut(data_obj,spec='SHMS')
                pcoll_cut_arrays = pcoll_cut()
                cuts_applied_list.append(pcoll_cut)
            if CURR:
                curr_cut = current_cut(data_obj,current='BCM4A')
                curr_cut_arrays = curr_cut()
                cuts_applied_list.append(curr_cut)
            if ZT:
                zt_cut = ztar_cut(data_obj)
                zt_cut_arrays = zt_cut()
                cuts_applied_list.append(zt_cut)
            if CT:
                ct_cut = CTime_cut(data_obj)
                ct_cut_arrays = ct_cut() 
                cuts_applied_list.append(ct_cut)
             
            #then make the cut arrays by getting the desired array from the DATA_INIT 
            # object and print stats for said cut, the result is a list of boolean arrays
            # for each cut, e.g. a Q2 cut will make a cut on the 'P.kin.primary.Q2' 
            # variable and the result will be a boolean array
            # note you will have a list of boolean arrays PER RUN, this is necessary to 
            # preserve the length of the arrays, since each run has a different # of array
            # elements.
            if STD:
                cuts_to_apply_list = []
                for cut in cuts_list:
                    br = data_obj.Branches[HCANA_names[cut.name]]
                    cut_array = cut(br)
                    # cut.stats()
                    
                    cuts_to_apply_list.append(cut_array)
                
            # add special cut arrays to cuts_to_apply_list
            if HCOLL: cuts_to_apply_list.append(hcoll_cut_arrays)
            if PCOLL: cuts_to_apply_list.append(pcoll_cut_arrays)
            if CURR: cuts_to_apply_list.append(curr_cut_arrays)
            if ZT: cuts_to_apply_list.append(zt_cut_arrays)
            if CT: cuts_to_apply_list.append(ct_cut_arrays)
                
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
            em_cut = Em[all_cuts]
            
            data_obj.Branches[m].update({'H.kin.secondary.pmiss_cut':Pm_cut,
                                           'H.kin.secondary.th_bq_cut':Th_rq_cut,
                                           'H.kin.secondary.emiss_nuc_cut':em_cut})
        
            print('Data object updated with the branches:\n',
                  'H.kin.secondary.pmiss_cut, H.kin.secondary.th_bq_cut,',
                  'H.kin.secondary.emiss_nuc_cut\n')
            print('Applied the following cuts:\n',
                  [c for c in cuts_applied_list])
            
            if show_stats:
                print('Cuts STATS:\n',
                      [c.stats() for c in cuts_applied_list])
            if save_report:
                filename = 'yield_n_xsec/current_cuts_DATA.txt'
                with open(filename,'w') as f:
                    f.write('CUTS REPORT\n')
                a = [c.report(filename) for c in cuts_applied_list] 