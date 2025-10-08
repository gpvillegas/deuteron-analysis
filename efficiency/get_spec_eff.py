#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 10:55:56 2025

@author: gvill
"""

import data_init as D
import LT.box as B
import cut_handler as C
import database_operations as db

import numpy as np


#%% CLASS to get detector efficiencies
"""
class get_spec_eff:
    data_obj: (DATA_INIT) this is a DATA_INIT object with HCANA variables loaded
                            from wanted runs or SIMC file
    eff_type: (string) efficiency type: 'HMS', 'SHMS', or 'LT', 
                The HMS and SHMS efficiencies are calculated as follows:
                    * define the events that should have made a track as those
                        that pass the cuts:
                        goodScinHit == 1 : requires a good scintillator hit 
                        0.5 > Betanotrk > 1.5 : beta calcualted without tracking
                        Caletotnorm > 0.6 for electrons or
                        0. < Caletotnorm < 0.6 for protons
                    
                    * require an additional cut to define those events that did
                    leave a track in the drift chamber.
                        DCntrack >= 1 : at least 1 track in DC
                        
                    * then the HMS/SHMS tracking efficiency is simply
                            trk_eff = did/should
                
                The Livetime calculated here is the total Livetime:
                    * number of events that had a EDTM TDC signal > 0.
                    divided by the total registered EDTM scaler events.
    curr: (string) which current reader are we using? 
                     'BCM1', 'BCM2', 'BCM4A', 'BCM4B', 'BCM4C'
                     this is to apply a current cut on events
    many: (bool) are we calcualting the efficiency of multiple runs? 
                (the data_obj will contain the list of runs if this is so)

"""
class get_spec_eff:
    def __init__(self,data_obj = None, eff_type = '', curr = 'BCM4A',many = False):
        """
            HMS, SHMS, LT efficiencies
        """
        self.eff_type = eff_type
        self.curr = curr
        if many:
            self.many = data_obj.many
        else:
            self.many = many
        
        #initialize variables needed for efficiency calculation
        self.init(data_obj)
        #initialize cuts
        self.init_cuts(data_obj)
        
        #apply current cut
        cut = self.curr_cut()
        scaler_cut = self.curr_cut.cut_array
        
        if eff_type == 'LT':
            if many:
                for m in data_obj.many:
                    self.var[m][0] = self.var[m][0][cut[m]]
                    self.var[m][1] = self.var[m][1][scaler_cut[m]]
                    self.var[m][2] = self.var[m][2][scaler_cut[m]]
                    self.var[m][3] = self.var[m][3][scaler_cut[m]]
            else:
                self.var[0] = self.var[0][cut]
                self.var[1] = self.var[1][scaler_cut]
                self.var[2] = self.var[2][scaler_cut]
                self.var[3] = self.var[3][scaler_cut]
        else:
            if many:
                for m in data_obj.many:
                    for i,v in enumerate(self.var[m]):
                        self.var[m][i] = v[cut[m]]
            else:
                for i,v in enumerate(self.var):
                    self.var[i] = v[cut]
    
    def calc_eff(self):
        if self.eff_type == 'SHMS' or self.eff_type == 'HMS':
            if self.many:
                self.trk_eff = {}
                self.trk_eff_err = {}
                self.res = {}
                for m in self.many:
                    n_goodScinHit = self.cuts[0](self.var[m][0])
                    n_Betanotrk = self.cuts[1](self.var[m][1])
                    n_Caletotnorm = self.cuts[2](self.var[m][2])
                    n_DCntrack = self.cuts[3](self.var[m][3])
                        
                    good_should_list = [n_goodScinHit,
                                             n_Betanotrk,
                                             n_Caletotnorm]
                                    
                    first = True
                    for l in good_should_list:
                        if first:
                            good_should = l
                            first = False
                        else:
                            good_should = good_should & l
                    
                    good_did = good_should & n_DCntrack
                    
                    te = good_did.sum()/good_should.sum() 
                    teerr =  np.sqrt(good_should.sum() -\
                                                good_did.sum())/\
                                                good_should.sum()
                    self.trk_eff[m] = te 
                    self.trk_eff_err[m] =  teerr
                    self.res[m] = (te, teerr)
                return self.trk_eff
            else:
                n_goodScinHit = self.cuts[0](self.var[0])
                n_Betanotrk = self.cuts[1](self.var[1])
                n_Caletotnorm = self.cuts[2](self.var[2])
                n_DCntrack = self.cuts[3](self.var[3])
                    
                good_should_list = [n_goodScinHit,
                                         n_Betanotrk,
                                         n_Caletotnorm]
                                
                first = True
                for l in good_should_list:
                    if first:
                        good_should = l
                        first = False
                    else:
                        good_should = good_should & l
                
                good_did = good_should & n_DCntrack
                
                self.trk_eff = good_did.sum()/good_should.sum() 
                self.trk_eff_err =  np.sqrt(good_should.sum() -\
                                            good_did.sum())/\
                                            good_should.sum()
                return self.trk_eff                            
        
        elif self.eff_type == 'LT':
            if self.many:
                self.tLT = {}
                tLT_corr = {}
                for m in self.many:
                    pEDTM_scalerRate = self.var[m][1]/self.var[m][2]
                    pTRIG6_scalerRate = (self.var[m][3]-self.var[m][1])/\
                                            self.var[m][2]

                    noEDTMcut_array = self.cuts[0](self.var[m][0])
                    EDTMcut_array = ~self.cuts[0](self.var[m][0])

                    tLT_corr_factor = 1 - (pTRIG6_scalerRate[-1] +\
                                           pEDTM_scalerRate[-1])*250e-06 +\
                                            pTRIG6_scalerRate[-1]*250e-06*\
                                            (1 + pEDTM_scalerRate[-1]/\
                                             (pTRIG6_scalerRate[-1]+\
                                              pEDTM_scalerRate[-1]))

                    self.tLT[m] = EDTMcut_array.sum()/self.var[m][1][-1]

                    tLT_corr[m] = self.tLT[m]*tLT_corr_factor
                      
                return self.tLT
            
            else:
                pEDTM_scalerRate = self.var[1]/self.var[2]
                pTRIG6_scalerRate = (self.var[3]-self.var[1])/self.var[2]
    
                noEDTMcut_array = self.cuts[0](self.var[0])
                EDTMcut_array = ~noEDTMcut_array
    
                tLT_corr_factor = 1 - (pTRIG6_scalerRate[-1] +\
                                       pEDTM_scalerRate[-1])*250e-06 +\
                                        pTRIG6_scalerRate[-1]*250e-06*\
                                        (1 + pEDTM_scalerRate[-1]/\
                                         (pTRIG6_scalerRate[-1]+pEDTM_scalerRate[-1]))
    
                self.tLT = EDTMcut_array.sum()/self.var[1][-1]
    
                tLT_corr = tLT*tLT_corr_factor
                  
                return self.tLT
            
    def init(self,data_obj):
        if self.eff_type == 'SHMS':
            if self.many:
                self.var = {}
                for m in self.many:
                    v = []
                    pgoodScinHit = data_obj.Branches[m]['P.hod.goodscinhit']
                    pBetanotrk = data_obj.Branches[m]['P.hod.betanotrack']
                    pCaletotnorm = data_obj.Branches[m]['P.cal.etotnorm']
                    pDCntrack = data_obj.Branches[m]['P.dc.ntrack']
                    
                    v.append(pgoodScinHit)   #0
                    v.append(pBetanotrk)     #1  
                    v.append(pCaletotnorm)   #2
                    v.append(pDCntrack)      #3
                    
                    self.var[m] = v
                    
            else:
                self.var = []
                pgoodScinHit = data_obj.Branches['P.hod.goodscinhit']
                pBetanotrk = data_obj.Branches['P.hod.betanotrack']
                pCaletotnorm = data_obj.Branches['P.cal.etotnorm']
                pDCntrack = data_obj.Branches['P.dc.ntrack']
                
                self.var.append(pgoodScinHit)   #0
                self.var.append(pBetanotrk)     #1  
                self.var.append(pCaletotnorm)   #2
                self.var.append(pDCntrack)      #3
        elif self.eff_type == 'HMS':
            if self.many:
                self.var = {}
                for m in self.many:
                    v = []
                    hgoodScinHit = data_obj.Branches[m]['H.hod.goodscinhit']
                    hBetanotrk = data_obj.Branches[m]['H.hod.betanotrack']
                    hCaletotnorm = data_obj.Branches[m]['H.cal.etotnorm']
                    hDCntrack = data_obj.Branches[m]['H.dc.ntrack']
                    
                    v.append(hgoodScinHit)   #0
                    v.append(hBetanotrk)     #1  
                    v.append(hCaletotnorm)   #2
                    v.append(hDCntrack)      #3
                    
                    self.var[m] = v
            else:
                self.var = []
                hgoodScinHit = data_obj.Branches['H.hod.goodscinhit']
                hBetanotrk = data_obj.Branches['H.hod.betanotrack']
                hCaletotnorm = data_obj.Branches['H.cal.etotnorm']
                hDCntrack = data_obj.Branches['H.dc.ntrack']
                
                self.var.append(hgoodScinHit)   #0
                self.var.append(hBetanotrk)     #1  
                self.var.append(hCaletotnorm)   #2
                self.var.append(hDCntrack)      #3
                
        elif self.eff_type == 'LT':
            if self.many:
                self.var = {}
                for m in self.many:
                    v = []
                    pEDTM_tdcTimeRaw = data_obj.Branches[m]['T.coin.pEDTM_tdcTimeRaw']
                    pEDTM_scaler = data_obj.Branches[m]['P.EDTM.scalerCut']
                    scalerTime = data_obj.Branches[m]['P.1MHz.scalerTimeCut']
                    pTRIG6_scaler = data_obj.Branches[m]['P.pTRIG6.scalerCut']
                    
                    v.append(pEDTM_tdcTimeRaw)  #0
                    v.append(pEDTM_scaler)      #1
                    v.append(scalerTime)        #2
                    v.append(pTRIG6_scaler)     #3
                    
                    self.var[m] = v
            else:
                self.var = []
                pEDTM_tdcTimeRaw = data_obj.Branches['T.coin.pEDTM_tdcTimeRaw']
                pEDTM_scaler = data_obj.Branches['P.EDTM.scalerCut']
                scalerTime = data_obj.Branches['P.1MHz.scalerTimeCut']
                pTRIG6_scaler = data_obj.Branches['P.pTRIG6.scalerCut']
                
                self.var.append(pEDTM_tdcTimeRaw)   #0
                self.var.append(pEDTM_scaler)       #1
                self.var.append(scalerTime)         #2
                self.var.append(pTRIG6_scaler)      #3
                
        else:
            print('No efficiency type chosen.')
            
    
    def init_cuts(self,data_obj):
        #initialize cuts
        self.cuts = []
        if self.eff_type == 'SHMS':
            pgoodScinHit_cut = C.VCUT(1.,'pgoodScinHit_cut')
            pBetanotrk_cut = C.WCUT(0.5,1.5,'pBetanotrk_cut')
            pCaletotnorm_cut = C.WCUT(0.6,np.inf,'pCaletotnorm_cut')
            pDCntrack_cut = C.WCUT(0.5,np.inf,'pDCntrack')
            
            self.cuts.append(pgoodScinHit_cut)  #0
            self.cuts.append(pBetanotrk_cut)    #1
            self.cuts.append(pCaletotnorm_cut)  #2
            self.cuts.append(pDCntrack_cut)     #3
            
        elif self.eff_type == 'HMS':
            hgoodScinHit_cut = C.VCUT(1.,'hgoodScinHit_cut')
            hBetanotrk_cut = C.WCUT(0.5,1.5,'hBetanotrk_cut')
            hCaletotnorm_cut = C.WCUT(0.,0.6,'hCaletotnorm_cut')
            hDCntrack_cut = C.WCUT(0.5,np.inf,'hDCntrack')
            
            self.cuts.append(hgoodScinHit_cut)  #0
            self.cuts.append(hBetanotrk_cut)    #1
            self.cuts.append(hCaletotnorm_cut)  #2
            self.cuts.append(hDCntrack_cut)     #3
            
        elif self.eff_type == 'LT':
            TRIG_OFF_cut = C.VCUT(0,name='TRIG_OFF_cut')
            
            self.cuts.append(TRIG_OFF_cut)  #0
        else:
            print('No efficiency type given.')


        #initialize current cut
        self.curr_cut = C.current_cut(data_obj,current=self.curr,
                                      many=self.many)



      
        

#%%

T_sel = ['H.kin.secondary.emiss','H.kin.secondary.pmiss',
         'H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th',
         'P.gtr.ph','H.extcor.xsieve','H.extcor.ysieve',
         'P.extcor.xsieve','P.extcor.ysieve','H.react.z','P.react.z',
         'CTime.epCoinTime_ROC2','H.hod.goodscinhit','H.hod.betanotrack',
         'H.cal.etotnorm','H.cer.npeSum','H.dc.ntrack','P.hod.goodscinhit',
         'P.hod.betanotrack','P.cal.etotnorm','P.hgcer.npeSum',
         'P.ngcer.npeSum','P.dc.ntrack','T.coin.pTRIG1_ROC2_tdcTimeRaw',
         'T.coin.pTRIG2_ROC2_tdcTimeRaw','T.coin.pTRIG6_ROC2_tdcTimeRaw',
         'T.coin.pEDTM_tdcTimeRaw']

TSP_sel = ['evNumber','P.BCM4A.scalerCharge','P.BCM4A.scalerChargeCut',
           'P.BCM4A.scalerCurrent','P.BCM4B.scalerChargeCut',
           'P.BCM4B.scalerCurrent','P.BCM4C.scalerChargeCut',
           'P.BCM4C.scalerCurrent','P.pTRIG1.scalerCut','P.pTRIG2.scalerCut',
           'P.pTRIG6.scalerCut','P.1MHz.scalerTimeCut','P.EDTM.scalerCut',
           'P.S1X.scalerCut']

TSH_sel = ['H.S1X.scalerCut','H.1MHz.scalerTimeCut']

R = 20851
RS = [20840,20841,20846,20851,20858,20861,20868,20869] #heep_coin
# R = [20871,20872]
# ,21065,21066,21067,21068,21069,21070, 
#      21071,21072,21073,21074,21075] #pm_120
# R = [20873,20874,20875,20876,20877,20878,20880,20881,20882,20883,21076,
#      21078,21079,21080,21081,21082,21083,21084,21085,21086,21087,21088,
#      21089,21090,21091,21092,21093,21094,21095,21096,21097,21098,21099,
#      21100,21101,21102] # pm_580
# R = [20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,20896,
#      20897,20898,20899,20900,20901,20902,20903,20904,20905,20907,20908,
#      20909,20910,20911,20912,20913,20914,20915,20916,20921,20922,20923,
#      20924,20925,20926,20927,20928,20929,20930,20931,20932,20933,20934,
#      20935,20936,20937,20938,20939,20940,20941,20942,20943,20944,20945,
#      20949,20950,20951,20953,20954,20955,20956] # pm_800

R = [20840,20841,20846,20851,20858,20861,20868,20869,20871,20872,
     20873,20874,20875,20876,20877,20878,20880,20881,20882,20883]

# R = [20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,
#      20896,20897,20898,20899,20900,20901,20902,20903,20904,20905]

# R = [20907,20908,20909,20910,20911,20912,20913,20914,20915,20916,
#      20921,20922,20923,20924,20925,20926,20927,20928,20929,20930,
#      20931,20932,20933,20934,20935,20936,20937,20938,20939,20940,
#      20941,20942,20943,20944,20945,20949,20950,20951,20953,20954]

# R = [20955,20956,20958,20959,20960,20961,20962,20963,20965,20966,
#      20969,20970,20971,20972,20973,20974,20975,20976,20977,20978,
#      20979,20980,20981,20982,20983,20984,20985,20986,20987,20988,
#      20989,20990,20991,20992,20993,20994,20995,20996,20997,20998]

# R = [20999,21000,21001,21002,21003,21004,21005,21006,21007,21008,
#      21009,21011,21012,21013,21014,21015,21016,21017,21018,21019,
#      21020,21021,21022,21023,21024,21025,21026,21027,21028,21029,
#      21030,21031,21032,21033,21034,21036,21037,21038,21039,21040]

# R = [21041,21042,21043,21044,21045,21046,21047,21048,21065,21066,
#      21067,21068,21069,21070,21071,21072,21073,21074,21075,21076,
#      21077,21078,21079,21080,21081,21082,21083,21084,21085,21086,
#      21087,21088,21089,21090,21091,21092,21093,21094,21095,21096,
#      21097,21098,21099,21100,21101,21102]

rdir = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

RUN = D.DATA_INIT(data_type='deut23_data', run= R, 
                  select_branches= {'T':T_sel,'TSP':TSP_sel,'TSH':TSH_sel},
                  select_trees=['T','TSP','TSH'],
                  ROOTfiles_path=rdir)
#%%
evNumber = RUN.Branches['evNumber'] 
BCM4A_current = RUN.Branches['P.BCM4A.scalerCurrent'] 
BCM4B_current = RUN.Branches['P.BCM4B.scalerCurrent']

currents = [BCM4A_current,BCM4B_current]

max_curr = max(BCM4A_current)

h = B.histo(BCM4A_current,range=(max_curr-5,max_curr+5),bins=50)
h.fit(plot_fit=False)
h.fit(h.mean.value-3*h.sigma.value,h.mean.value+3*h.sigma.value,plot_fit=False)

mean = h.mean.value
sigma = h.sigma.value

curr_cut = C.WCUT(mean-3.5*sigma,np.inf)

BCM4A_cut = curr_cut(BCM4A_current)
BCM4A_cut = np.roll(BCM4A_cut,-1)[:-1] # we ignore the first element in the 
                                        # cut so it can be applied to ranges
                                        # in ranges it is the next number that
                                        # determines if it is a good range or 
                                        # not

# create array of intervals using the event numbers and the same array rolled
# over by 1
intervals = np.vstack((evNumber,np.roll(evNumber,-1)))
u = intervals.T[:-1] # then take the transpose to have an array where each line
                     # is a range, the last range is [last item, first item]
                     # so its not needed.
                     
u = u[BCM4A_cut] # apply the cut to the intervals to only get the good intervals                      


Pm = RUN.Branches['H.kin.secondary.pmiss']
Em = RUN.Branches['H.kin.secondary.emiss']

Pm_curr_cut = []
Em_curr_cut = []
for i in u:
    for j in range(int(i[0])-1,int(i[1])):
        Pm_curr_cut.append(Pm[j])
        Em_curr_cut.append(Em[j])
        
#test current_cut class
c_cut = C.current_cut(RUN)
c_cut_arr = c_cut()


cuts_list = C.acceptance_cuts

cuts_to_apply = []
for cut in cuts_list:
    br = RUN.Branches[C.HCANA_names[cut.name]]
    br_curr_cut = []
    for i in u:
        for j in range(int(i[0])-1,int(i[1])):
            br_curr_cut.append(br[j])
    cut_array = cut(np.array(br_curr_cut))
    cut.stats()
    
    cuts_to_apply.append(cut_array)

all_cuts = cuts_to_apply[0]
for arr in cuts_to_apply:    
    all_cuts = all_cuts & arr
    
cuts_to_apply = []
for cut in cuts_list:
    br = RUN.Branches[C.HCANA_names[cut.name]]
    cut_array = cut(np.array(br))
    cut.stats()
    
    cuts_to_apply.append(cut_array)

all_cuts_nocurrcut = cuts_to_apply[0]
for arr in cuts_to_apply:    
    all_cuts_nocurrcut = all_cuts_nocurrcut & arr    


Pm_c = np.array(Pm_curr_cut)[all_cuts]
Em_c = np.array(Em_curr_cut)[all_cuts]

Em_Pm_h = B.histo2d(Pm_c,Em_c,range=[(-0.025,0.04),(-0.025,0.02)],bins=100)
Em_Pm_nocurrcut_h = B.histo2d(Pm[all_cuts_nocurrcut], Em[all_cuts_nocurrcut],
                     range=[(-0.025,0.04),(-0.025,0.02)],bins=100)

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
           
        proj_values = h2.x_bin_center[bin_range]
        yield_values = counts[:,0]
        yield_errors = counts[:,1]
    
    else:
        counts = []
        for nx in bin_range[:]:
            h = norm*h2.project_y(bins = [nx])
            if plot_proj:
                B.pl.figure()
                h.plot()
                B.pl.vlines([proj_min,proj_max],0,30,
                            linestyles='--',color='black')
                
            counts.append(h.sum(proj_min,proj_max))
        
        counts = np.array(counts)
           
        proj_values = h2.y_bin_center[bin_range]
        yield_values = counts[:,0]
        yield_errors = counts[:,1]
    
    return (proj_values,yield_values,yield_errors)

Pm_values, YIELD, YIELD_E = get_yield(Em_Pm_h, bin_range=np.arange(0,100,1), 
                                      proj_min=-0.02, proj_max=0.02) 

Pm_values_nocurrcut, YIELD_nocurrcut, YIELD_E_nocurrcut = get_yield(\
                            Em_Pm_nocurrcut_h, bin_range=np.arange(0,100,1), 
                            proj_min=-0.02, proj_max=0.02) 

#%% calc hms/shms tracking efficiencies
pgoodScinHit_cut = C.VCUT(1.,'pgoodScinHit_cut')
pBetanotrk_cut = C.WCUT(0.5,1.5,'pBetanotrk_cut')
pCaletotnorm_cut = C.WCUT(0.6,np.inf,'pCaletotnorm_cut')
# pNGCer_npeSum_cut = C.WCUT(0.5,np.inf,'pNGCer_npeSum_cut')
pDCntrack_cut = C.WCUT(0.5,np.inf,'pDCntrack')
hgoodScinHit_cut = C.VCUT(1.,'hgoodScinHit_cut')
hBetanotrk_cut = C.WCUT(0.5,1.5,'hBetanotrk_cut')
hCaletotnorm_cut = C.WCUT(0.,0.6,'hCaletotnorm_cut')
# pNGCer_npeSum_cut = C.WCUT(0.5,np.inf,'pNGCer_npeSum_cut')
hDCntrack_cut = C.WCUT(0.5,np.inf,'hDCntrack')

#apply current cut
curr_cut = C.current_cut(RUN,current='BCM4A',many=True)
curr_cut_arrays = curr_cut()

cuts_list = [C.hms_delta,C.shms_delta]

# cuts_to_apply = {}
# for m in RUN.many:
#     clist = []
#     for cut in cuts_list:
#         br = RUN.Branches[m][C.HCANA_names[cut.name]]
#         cut_array = cut(br)
#         # cut.stats()        
#         clist.append(cut_array)
#     clist.append(curr_cut_arrays[m])
#     cuts_to_apply[m] = clist
    
# all_cuts = {}
# for m in RUN.many:
#     all_cuts[m] = cuts_to_apply[m][0]
#     for arr in cuts_to_apply[m]:    
#         all_cuts[m] = all_cuts[m] & arr     
    
ptrk = []
htrk = []
for m in RUN.many:
    pgoodScinHit = RUN.Branches[m]['P.hod.goodscinhit']
    pBetanotrk = RUN.Branches[m]['P.hod.betanotrack']
    pCaletotnorm = RUN.Branches[m]['P.cal.etotnorm']
    pNGCer_npeSum = RUN.Branches[m]['P.ngcer.npeSum']
    pHGCer_npeSum = RUN.Branches[m]['P.hgcer.npeSum']
    pDCntrack = RUN.Branches[m]['P.dc.ntrack']
    
    pgoodScinHit = pgoodScinHit[curr_cut_arrays[m]]
    pBetanotrk = pBetanotrk[curr_cut_arrays[m]]
    pCaletotnorm = pCaletotnorm[curr_cut_arrays[m]]
    pNGCer_npeSum = pNGCer_npeSum[curr_cut_arrays[m]]
    pHGCer_npeSum = pHGCer_npeSum[curr_cut_arrays[m]]
    pDCntrack = pDCntrack[curr_cut_arrays[m]]
    
    n_pgoodScinHit = pgoodScinHit_cut(pgoodScinHit)
    n_pBetanotrk = pBetanotrk_cut(pBetanotrk)
    n_pCaletotnorm = pCaletotnorm_cut(pCaletotnorm)
    n_pDCntrack = pDCntrack_cut(pDCntrack)
    
    good_elec_should_list = [n_pgoodScinHit,n_pBetanotrk,n_pCaletotnorm]
    
    first = True
    for a in good_elec_should_list:
        if first:
            good_elec_should = a
            first = False
        else:
            good_elec_should = good_elec_should & a
    
    good_elec_did = good_elec_should & n_pDCntrack
    
    ptrk_eff = good_elec_did.sum()/good_elec_should.sum() 
    ptrk_eff_err =  np.sqrt(good_elec_should.sum() - good_elec_did.sum())/\
                                good_elec_should.sum()
    
    print('ptrk_eff = ', ptrk_eff, '$\pm$', ptrk_eff_err)
    ptrk.append((ptrk_eff,ptrk_eff_err))
    
    hgoodScinHit = RUN.Branches[m]['H.hod.goodscinhit']
    hBetanotrk = RUN.Branches[m]['H.hod.betanotrack']
    hCaletotnorm = RUN.Branches[m]['H.cal.etotnorm']
    hCer_npeSum = RUN.Branches[m]['H.cer.npeSum']
    hDCntrack = RUN.Branches[m]['H.dc.ntrack']
    
    hgoodScinHit = hgoodScinHit[curr_cut_arrays[m]]
    hBetanotrk = hBetanotrk[curr_cut_arrays[m]]
    hCaletotnorm = hCaletotnorm[curr_cut_arrays[m]]
    hCer_npeSum = hCer_npeSum [curr_cut_arrays[m]]
    hDCntrack = hDCntrack[curr_cut_arrays[m]]
    
    n_hgoodScinHit = hgoodScinHit_cut(hgoodScinHit)
    n_hBetanotrk = hBetanotrk_cut(hBetanotrk)
    n_hCaletotnorm = hCaletotnorm_cut(hCaletotnorm)
    n_hDCntrack = hDCntrack_cut(hDCntrack)
    
    good_had_should_list = [n_hgoodScinHit,n_hBetanotrk,n_hCaletotnorm]
    
    first = True
    for a in good_had_should_list:
        if first:
            good_had_should = a
            first = False
        else:
            good_had_should = good_had_should & a
    
    good_had_did = good_had_should & n_hDCntrack
    
    htrk_eff = good_had_did.sum()/good_had_should.sum() 
    htrk_eff_err =  np.sqrt(good_had_should.sum() - good_had_did.sum())/\
                                good_had_should.sum()
    
    print('htrk_eff = ', htrk_eff, '$\pm$', htrk_eff_err)
    htrk.append((htrk_eff,htrk_eff_err)) 

# B.pl.figure()
# B.plot_exp(RUN.many,np.array(ptrk)[:,0],np.array(ptrk)[:,1],label='shms_trk') 
# B.plot_exp(RUN.many,np.array(htrk)[:,0],np.array(htrk)[:,1],label='hms_trk')
# xlo,xhi = B.pl.xlim()
# B.pl.hlines(np.mean(np.array(ptrk)[:,0]),xlo,xhi,colors='black',linestyles='--') 
# B.pl.hlines(np.mean(np.array(htrk)[:,0]),xlo,xhi,colors='black',linestyles='--') 
# B.pl.text(20842, 0.9892, '0.9891')
# B.pl.text(20842, 0.9839, '0.9837')
# B.pl.legend()
#%%
#get S1X scaler rates
pS1X_scalerRate = []
hS1X_scalerRate = []
current = []
for m in RUN.many:
    ps1x = RUN.Branches[m]['P.S1X.scalerCut'][-1]
    pstime = RUN.Branches[m]['P.1MHz.scalerTimeCut'][-1]
    
    hs1x = RUN.Branches[m]['H.S1X.scalerCut'][-1]
    hstime = RUN.Branches[m]['H.1MHz.scalerTimeCut'][-1]
    
    ps1xrate = ps1x/pstime/1000. #kHz
    hs1xrate = hs1x/hstime/1000. #kHz
    
    pS1X_scalerRate.append(ps1xrate)
    hS1X_scalerRate.append(hs1xrate)
    
    charge = RUN.Branches[m]['P.BCM4A.scalerChargeCut'][-1]
    
    curr = charge/pstime
    
    current.append(curr)

B.pl.figure()
B.plot_exp(pS1X_scalerRate,np.array(ptrk)[:,0],np.array(ptrk)[:,1],label='shms_trk') 
B.pl.legend()
B.pl.xlabel('S1X_scalerRate')
B.pl.ylabel('tracking efficiency')
B.pl.figure()
B.plot_exp(hS1X_scalerRate,np.array(htrk)[:,0],np.array(htrk)[:,1],label='hms_trk')
B.pl.legend()
B.pl.xlabel('S1X_scalerRate')
B.pl.ylabel('tracking efficiency')
B.pl.figure()
B.plot_exp(current,np.array(ptrk)[:,0],np.array(ptrk)[:,1],label='shms_trk') 
B.plot_exp(current,np.array(htrk)[:,0],np.array(htrk)[:,1],label='hms_trk')    
B.pl.legend()
B.pl.xlabel('BCM4A_current')
B.pl.ylabel('tracking efficiency')

B.pl.figure()
B.plot_exp(current,pS1X_scalerRate) 
B.pl.xlabel('BCM4A_current')
B.pl.ylabel('SHMS S1X scaler rate')
B.pl.figure()
B.plot_exp(current,hS1X_scalerRate)    
B.pl.xlabel('BCM4A_current')
B.pl.ylabel('HMS S1X scaler rate')

#%%
def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    #B.pl.legend(fontsize = fsize)

h_pgoodScinHit = B.histo(pgoodScinHit,range=(-0.5,2.5),bins=3)
B.pl.figure()
h_pgoodScinHit.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(1.0, 0., yhi,color='black',linestyles='--')
pptxify('pgoodScinHit')

h_pBetanotrk = B.histo(pBetanotrk,range=(0.1,2.0),bins=100)
B.pl.figure()
h_pBetanotrk.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines([0.5,1.5], 0., yhi,color='black',linestyles='--')
pptxify('pBetanotrk')

h_pCaletotnorm = B.histo(pCaletotnorm,range=(0.5,1.5),bins=100)
B.pl.figure()
h_pCaletotnorm.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(0.6, 0., yhi,color='black',linestyles='--')
pptxify('pCaletotnorm')

h_pDCntrack = B.histo(pDCntrack,range=(-0.5,5.5),bins=6)
B.pl.figure()
h_pDCntrack.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(0.5, 0., yhi,color='black',linestyles='--')
pptxify('pDCntrack')

h_hgoodScinHit = B.histo(hgoodScinHit,range=(-0.5,2.5),bins=3)
B.pl.figure()
h_hgoodScinHit.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(1.0, 0., yhi,color='black',linestyles='--')
pptxify('hgoodScinHit')

h_hBetanotrk = B.histo(hBetanotrk,range=(0.5,1.5),bins=100)
B.pl.figure()
h_hBetanotrk.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines([0.5,1.5], 0., yhi,color='black',linestyles='--')
pptxify('hBetanotrk')

h_hCaletotnorm = B.histo(hCaletotnorm,range=(0.,1.),bins=100)
B.pl.figure()
h_hCaletotnorm.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines([0.,0.6], 0., yhi,color='black',linestyles='--')
pptxify('hCaletotnorm')

h_hDCntrack = B.histo(hDCntrack,range=(-0.5,5.5),bins=6)
B.pl.figure()
h_hDCntrack.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(0.5, 0., yhi,color='black',linestyles='--')
pptxify('hDCntrack')

#%% calculate live time

pEDTM_tdcTimeRaw = RUN.Branches['T.coin.pEDTM_tdcTimeRaw']
pTRIG6_tdcTimeRaw = RUN.Branches['T.coin.pTRIG6_ROC2_tdcTimeRaw']
pEDTM_scaler = RUN.Branches['P.EDTM.scalerCut']
pTRIG6_scaler = RUN.Branches['P.pTRIG6.scalerCut']
scalerTime = RUN.Branches['P.pTRIG6.scalerCut']

pTRIG6_scalerRateRaw = pTRIG6_scaler/scalerTime
pEDTM_scalerRate = pEDTM_scaler/scalerTime

pTRIG6_scalerRate = pTRIG6_scalerRateRaw - pEDTM_scalerRate

TRIG_OFF_cut = C.VCUT(0,name='TRIG_OFF_cut')

noEDTMcut_array = TRIG_OFF_cut(pEDTM_tdcTimeRaw)
EDTMcut_array = ~TRIG_OFF_cut(pEDTM_tdcTimeRaw)
pTRIG6cut_array = ~TRIG_OFF_cut(pTRIG6_tdcTimeRaw)

#apply current cut
curr_cut = C.current_cut(RUN,current='BCM4A')
curr_cut_arrays = curr_cut()

noEDTMcut_array = noEDTMcut_array[curr_cut_arrays]
EDTMcut_array = EDTMcut_array[curr_cut_arrays]
pTRIG6cut_array = pTRIG6cut_array[curr_cut_arrays]

pEDTM_scalerRate = pEDTM_scalerRate[curr_cut.cut_array]
pTRIG6_scalerRateRaw = pTRIG6_scalerRateRaw[curr_cut.cut_array]
pTRIG6_scalerRate = pTRIG6_scalerRate[curr_cut.cut_array]

tLT_corr_factor = 1 - pTRIG6_scalerRateRaw[-1]*250e-06 +\
    pTRIG6_scalerRate[-1]*250e-06*(1 + pEDTM_scalerRate[-1]/pTRIG6_scalerRateRaw[-1])

tLT = EDTMcut_array.sum()/pEDTM_scaler[curr_cut.cut_array][-1]

tLT_corr = tLT*tLT_corr_factor
  
pTRIG6_noEDTM = noEDTMcut_array & pTRIG6cut_array

cpuLT = pTRIG6_noEDTM.sum()/\
                (pTRIG6_scaler[curr_cut.cut_array][-1] -\
                         pEDTM_scaler[curr_cut.cut_array][-1])   
#%% write charge and current values to db

deutDB_name= 'deuteron_db.db'

hms_trk_eff = get_spec_eff(RUN,eff_type='HMS',many=True)
htrk = hms_trk_eff.calc_eff()

shms_trk_eff = get_spec_eff(RUN,eff_type='SHMS',many=True)
ptrk= shms_trk_eff.calc_eff()

LT_eff = get_spec_eff(RUN,eff_type='LT',many=True)
LT= LT_eff.calc_eff()
 
#%%
run =  RUN.many
for r in run:
    BCM4A_charge = RUN.Branches[r]['P.BCM4A.scalerChargeCut'][-1]/1e3
    BCM4A_current = RUN.Branches[r]['P.BCM4A.scalerChargeCut'][-1]/\
        RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]
    
    BCM4B_charge = RUN.Branches[r]['P.BCM4B.scalerChargeCut'][-1]/1e3
    BCM4B_current = RUN.Branches[r]['P.BCM4B.scalerChargeCut'][-1]/\
        RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]
    
    BCM4C_charge = RUN.Branches[r]['P.BCM4C.scalerChargeCut'][-1]/1e3
    BCM4C_current = RUN.Branches[r]['P.BCM4C.scalerChargeCut'][-1]/\
        RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]
    
    T1_scl_rate =  (RUN.Branches[r]['P.pTRIG1.scalerCut'][-1]-\
                    RUN.Branches[r]['P.EDTM.scalerCut'][-1])/\
                    RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]/1e3
    T2_scl_rate = (RUN.Branches[r]['P.pTRIG2.scalerCut'][-1]-\
                    RUN.Branches[r]['P.EDTM.scalerCut'][-1])/\
                    RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]/1e3
    T6_scl_rate = (RUN.Branches[r]['P.pTRIG6.scalerCut'][-1]-\
                    RUN.Branches[r]['P.EDTM.scalerCut'][-1])/\
                    RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]/1e3
    
    ht = htrk[r]
    pt = ptrk[r]
    lt = LT[r]        
    
    col_names = ['BCM4A_charge','BCM4A_current','BCM4B_charge','BCM4B_current',
                 'BCM4C_charge','BCM4C_current','T1_scl_rates','T2_scl_rates',
                 'T6_scl_rates','tLT','HMS_TrkEff','SHMS_TrkEff']
    col_vals = [f'{BCM4A_charge:.3f}', f'{BCM4A_current:.3f}',
                f'{BCM4B_charge:.3f}', f'{BCM4B_current:.3f}',
                f'{BCM4C_charge:.3f}', f'{BCM4C_current:.3f}',
                f'{T1_scl_rate:.3f}',f'{T2_scl_rate:.3f}',
                f'{T6_scl_rate:.3f}',f'{lt:.4f}',f'{ht:.4f}',f'{pt:.4f}']
    
    # print([f'{nn} = {col_vals[i]}' for i,nn in enumerate(col_names)])
    db.update_row(deutDB_name, 'RUN_LIST_UPDATED',col_names, col_vals, 
                  where=f'run = {r}') 




#%%
               