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
from matplotlib import colormaps as cmp

import numpy as np

# Set default font to times new roman
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14
}
B.pl.rc('font', **font)

#%%
# To get rid of unhelpful groupings (lists with only 1 item, etc)       
def get_list(db_res, index = 0):
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0]

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
                
        elif eff_type == 'cpuLT':
            if many:
                for m in data_obj.many:
                    self.var[m][0] = self.var[m][0][cut[m]]
                    self.var[m][1] = self.var[m][1][cut[m]]
                    self.var[m][2] = self.var[m][2][cut[m]]
                    self.var[m][3] = self.var[m][3][cut[m]]
                    self.var[m][4] = self.var[m][4][scaler_cut[m]]
                    self.var[m][5] = self.var[m][5][scaler_cut[m]]
                    self.var[m][6] = self.var[m][6][scaler_cut[m]]
                    self.var[m][7] = self.var[m][7][scaler_cut[m]]
            else:
                self.var[0] = self.var[0][cut]
                self.var[1] = self.var[1][cut]
                self.var[2] = self.var[2][cut]
                self.var[3] = self.var[3][cut]
                self.var[4] = self.var[4][scaler_cut]
                self.var[5] = self.var[5][scaler_cut]
                self.var[6] = self.var[6][scaler_cut]
                self.var[7] = self.var[7][scaler_cut]
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
                    n_goodstarttime = self.cuts[4](self.var[m][4])
                    n_CerNpeSum = self.cuts[5](self.var[m][5])
                        
                    good_should_list = [n_goodScinHit,n_Betanotrk,
                                        n_Caletotnorm,n_goodstarttime,
                                        n_CerNpeSum]
                                    
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
                n_goodstarttime = self.cuts[4](self.var[4])
                n_CerNpeSum = self.cuts[5](self.var[5])
                    
                good_should_list = [n_goodScinHit,n_Betanotrk,
                                    n_Caletotnorm,n_goodstarttime,
                                    n_CerNpeSum]
                                
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
            
        elif self.eff_type == 'cpuLT':
            if self.many:
                cpuLT_list = {}
                for m in self.many:
                    noEDTMcut_array = self.cuts[0](self.var[m][0])
                    pTRIG1cut_array = ~self.cuts[0](self.var[m][1])
                    pTRIG2cut_array = ~self.cuts[0](self.var[m][2])
                    pTRIG6cut_array = ~self.cuts[0](self.var[m][3])
                    
                    pTRIG1_noEDTM = noEDTMcut_array & pTRIG1cut_array
                    pTRIG2_noEDTM = noEDTMcut_array & pTRIG2cut_array
                    pTRIG6_noEDTM = noEDTMcut_array & pTRIG6cut_array
                    
                    pTRIG1_scl_noEDTM = self.var[m][5][-1] - self.var[m][4][-1]
                    pTRIG2_scl_noEDTM = self.var[m][6][-1] - self.var[m][4][-1]
                    pTRIG6_scl_noEDTM = self.var[m][7][-1] - self.var[m][4][-1]
                    
                    T1cpuLT = (pTRIG1_noEDTM.sum()*self.var[m][8])/\
                                (pTRIG1_scl_noEDTM)
                    T2cpuLT = (pTRIG2_noEDTM.sum()*self.var[m][9])/\
                                (pTRIG2_scl_noEDTM)
                    T6cpuLT = (pTRIG6_noEDTM.sum()*self.var[m][10])/\
                                (pTRIG6_scl_noEDTM)
                    
                    # calculate Binomial error for cpuLT
                    T1cpuLT_err = np.sqrt(pTRIG1_noEDTM.sum()*\
                                    (1-(pTRIG1_noEDTM.sum()/pTRIG1_scl_noEDTM)))*\
                                        self.var[m][8]/pTRIG1_scl_noEDTM
                    T2cpuLT_err = np.sqrt(pTRIG2_noEDTM.sum()*\
                                    (1-(pTRIG2_noEDTM.sum()/pTRIG2_scl_noEDTM)))*\
                                        self.var[m][9]/pTRIG2_scl_noEDTM
                    T6cpuLT_err = np.sqrt(pTRIG6_noEDTM.sum()*\
                                    (1-(pTRIG6_noEDTM.sum()/pTRIG6_scl_noEDTM)))*\
                                        self.var[m][10]/pTRIG6_scl_noEDTM
                
                    if self.var[m][8] == -1:
                        T1cpuLT = 0.
                        T1cpuLT_err = 0.
                    if self.var[m][9] == -1:
                        T2cpuLT = 0.
                        T2cpuLT_err = 0.
                    if self.var[m][10] == -1:
                        T6cpuLT = 0.
                        T6cpuLT_err = 0.    
                
                    cpuLT = [(T1cpuLT,T1cpuLT_err),(T2cpuLT,T2cpuLT_err),
                             (T6cpuLT,T6cpuLT_err)]
                    
                    cpuLT_list[m] = cpuLT
                
                return cpuLT_list
            
            else:
                noEDTMcut_array = self.cuts[0](self.var[0])
                pTRIG1cut_array = ~self.cuts[0](self.var[1])
                pTRIG2cut_array = ~self.cuts[0](self.var[2])
                pTRIG6cut_array = ~self.cuts[0](self.var[3])
                
                pTRIG1_noEDTM = noEDTMcut_array & pTRIG1cut_array
                pTRIG2_noEDTM = noEDTMcut_array & pTRIG2cut_array
                pTRIG6_noEDTM = noEDTMcut_array & pTRIG6cut_array
                
                pTRIG1_scl_noEDTM = self.var[5][-1] - self.var[4][-1]
                pTRIG2_scl_noEDTM = self.var[6][-1] - self.var[4][-1]
                pTRIG6_scl_noEDTM = self.var[7][-1] - self.var[4][-1]
                
                T1cpuLT = (pTRIG1_noEDTM.sum()*self.var[8])/\
                            (pTRIG1_scl_noEDTM)
                T2cpuLT = (pTRIG2_noEDTM.sum()*self.var[9])/\
                            (pTRIG2_scl_noEDTM)
                T6cpuLT = (pTRIG6_noEDTM.sum()*self.var[10])/\
                            (pTRIG6_scl_noEDTM)
                            
                # calculate Binomial error for cpuLT
                T1cpuLT_err = np.sqrt(pTRIG1_noEDTM.sum()*\
                                (1-(pTRIG1_noEDTM.sum()/pTRIG1_scl_noEDTM)))*\
                                    self.var[8]/pTRIG1_scl_noEDTM
                T2cpuLT_err = np.sqrt(pTRIG2_noEDTM.sum()*\
                                (1-(pTRIG2_noEDTM.sum()/pTRIG2_scl_noEDTM)))*\
                                    self.var[9]/pTRIG2_scl_noEDTM
                T6cpuLT_err = np.sqrt(pTRIG6_noEDTM.sum()*\
                                (1-(pTRIG6_noEDTM.sum()/pTRIG6_scl_noEDTM)))*\
                                    self.var[10]/pTRIG6_scl_noEDTM
                                    
                if self.var[8] == -1:
                    T1cpuLT = 0.
                    T1cpuLT_err = 0.
                if self.var[9] == -1:
                    T2cpuLT = 0.
                    T2cpuLT_err = 0.
                if self.var[10] == -1:
                    T6cpuLT = 0.
                    T6cpuLT_err = 0.
            
                cpuLT = [(T1cpuLT,T1cpuLT_err),(T2cpuLT,T2cpuLT_err),
                         (T6cpuLT,T6cpuLT_err)]
                
                return cpuLT
                
                
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
                    pgoodstarttime = data_obj.Branches[m]['P.hod.goodstarttime']
                    pNGCerNpeSum = data_obj.Branches[m]['P.ngcer.npeSum']
                    
                    v.append(pgoodScinHit)   #0
                    v.append(pBetanotrk)     #1  
                    v.append(pCaletotnorm)   #2
                    v.append(pDCntrack)      #3
                    v.append(pgoodstarttime) #4
                    v.append(pNGCerNpeSum)   #5
                    
                    self.var[m] = v
                    
            else:
                self.var = []
                pgoodScinHit = data_obj.Branches['P.hod.goodscinhit']
                pBetanotrk = data_obj.Branches['P.hod.betanotrack']
                pCaletotnorm = data_obj.Branches['P.cal.etotnorm']
                pDCntrack = data_obj.Branches['P.dc.ntrack']
                pgoodstarttime = data_obj.Branches['P.hod.goodstarttime']
                pNGCerNpeSum = data_obj.Branches['P.ngcer.npeSum']

                
                self.var.append(pgoodScinHit)   #0
                self.var.append(pBetanotrk)     #1  
                self.var.append(pCaletotnorm)   #2
                self.var.append(pDCntrack)      #3
                self.var.append(pgoodstarttime) #4
                self.var.append(pNGCerNpeSum)   #5
                
        elif self.eff_type == 'HMS':
            if self.many:
                self.var = {}
                for m in self.many:
                    v = []
                    hgoodScinHit = data_obj.Branches[m]['H.hod.goodscinhit']
                    hBetanotrk = data_obj.Branches[m]['H.hod.betanotrack']
                    hCaletotnorm = data_obj.Branches[m]['H.cal.etotnorm']
                    hDCntrack = data_obj.Branches[m]['H.dc.ntrack']
                    hgoodstarttime = data_obj.Branches[m]['H.hod.goodstarttime']
                    hCerNpeSum = data_obj.Branches[m]['H.cer.npeSum']
                    
                    v.append(hgoodScinHit)   #0
                    v.append(hBetanotrk)     #1  
                    v.append(hCaletotnorm)   #2
                    v.append(hDCntrack)      #3
                    v.append(hgoodstarttime) #4
                    v.append(hCerNpeSum)     #5
                    
                    self.var[m] = v
            else:
                self.var = []
                hgoodScinHit = data_obj.Branches['H.hod.goodscinhit']
                hBetanotrk = data_obj.Branches['H.hod.betanotrack']
                hCaletotnorm = data_obj.Branches['H.cal.etotnorm']
                hDCntrack = data_obj.Branches['H.dc.ntrack']
                hgoodstarttime = data_obj.Branches['H.hod.goodstarttime']
                hCerNpeSum = data_obj.Branches['H.cer.npeSum']

                self.var.append(hgoodScinHit)   #0
                self.var.append(hBetanotrk)     #1  
                self.var.append(hCaletotnorm)   #2
                self.var.append(hDCntrack)      #3
                self.var.append(hgoodstarttime) #4
                self.var.append(hCerNpeSum)     #5
                
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
                    
        elif self.eff_type == 'cpuLT':
            if self.many:
                self.var = {}
                for m in self.many:
                    v = []
                    PS1, PS2, PS6 =\
                        db.retrieve('deuteron_db.db','PS1, PS2, PS6',
                                    'RUN_LIST_UPDATED', 
                                    where = f"run=\'{m}\'")[0]
                    pEDTM_tdcTimeRaw = data_obj.Branches[m]['T.coin.pEDTM_tdcTimeRaw']
                    pTRIG1_tdcTimeRaw = data_obj.Branches[m]['T.coin.pTRIG1_ROC2_tdcTimeRaw']
                    pTRIG2_tdcTimeRaw = data_obj.Branches[m]['T.coin.pTRIG2_ROC2_tdcTimeRaw']
                    pTRIG6_tdcTimeRaw = data_obj.Branches[m]['T.coin.pTRIG6_ROC2_tdcTimeRaw']
                    pEDTM_scaler = data_obj.Branches[m]['P.EDTM.scalerCut']
                    pTRIG1_scaler = data_obj.Branches[m]['P.pTRIG1.scalerCut']
                    pTRIG2_scaler = data_obj.Branches[m]['P.pTRIG2.scalerCut']
                    pTRIG6_scaler = data_obj.Branches[m]['P.pTRIG6.scalerCut']
                    
                    v.append(pEDTM_tdcTimeRaw)   #0
                    v.append(pTRIG1_tdcTimeRaw)  #1
                    v.append(pTRIG2_tdcTimeRaw)  #2
                    v.append(pTRIG6_tdcTimeRaw)  #3
                    v.append(pEDTM_scaler)       #4
                    v.append(pTRIG1_scaler)      #5
                    v.append(pTRIG2_scaler)      #6
                    v.append(pTRIG6_scaler)      #7
                    v.append(PS1)                #8
                    v.append(PS2)                #9
                    v.append(PS6)                #10
                    
                    self.var[m] = v
            else:
                self.var = []
                PS1, PS2, PS6 =\
                    db.retrieve('deuteron_db.db','PS1, PS2, PS6',
                                'RUN_LIST_UPDATED', 
                                where = f"run=\'{data_obj.many}\'")[0]
                pEDTM_tdcTimeRaw = data_obj.Branches['T.coin.pEDTM_tdcTimeRaw']
                pTRIG1_tdcTimeRaw = data_obj.Branches['T.coin.pTRIG1_ROC2_tdcTimeRaw']
                pTRIG2_tdcTimeRaw = data_obj.Branches['T.coin.pTRIG2_ROC2_tdcTimeRaw']
                pTRIG6_tdcTimeRaw = data_obj.Branches['T.coin.pTRIG6_ROC2_tdcTimeRaw']
                pEDTM_scaler = data_obj.Branches['P.EDTM.scalerCut']
                pTRIG1_scaler = data_obj.Branches['P.pTRIG1.scalerCut']
                pTRIG2_scaler = data_obj.Branches['P.pTRIG2.scalerCut']
                pTRIG6_scaler = data_obj.Branches['P.pTRIG6.scalerCut']
                
                self.var.append(pEDTM_tdcTimeRaw)   #0
                self.var.append(pTRIG1_tdcTimeRaw)  #1
                self.var.append(pTRIG2_tdcTimeRaw)  #2
                self.var.append(pTRIG6_tdcTimeRaw)  #3
                self.var.append(pEDTM_scaler)       #4
                self.var.append(pTRIG1_scaler)      #5
                self.var.append(pTRIG2_scaler)      #6
                self.var.append(pTRIG6_scaler)      #7
                self.var.append(PS1)                #8
                self.var.append(PS2)                #9
                self.var.append(PS6)                #10
        else:
            print('No efficiency type chosen.')
            
    
    def init_cuts(self,data_obj):
        #initialize cuts
        self.cuts = []
        if self.eff_type == 'SHMS':
            pgoodScinHit_cut = C.VCUT(1.,'pgoodScinHit_cut')
            pBetanotrk_cut = C.WCUT(0.5,1.5,'pBetanotrk_cut')
            pCaletotnorm_cut = C.WCUT(0.6,np.inf,'pCaletotnorm_cut')
            pDCntrack_cut = C.WCUT(1.,np.inf,'pDCntrack')
            pgoodstarttime_cut = C.VCUT(1.,'pgoodstarttime_cut')
            pNGCerNpeSum_cut = C.WCUT(1.5,np.inf,'hNGCerNpeSum')
            
            self.cuts.append(pgoodScinHit_cut)      #0
            self.cuts.append(pBetanotrk_cut)        #1
            self.cuts.append(pCaletotnorm_cut)      #2
            self.cuts.append(pDCntrack_cut)         #3
            self.cuts.append(pgoodstarttime_cut)    #4
            self.cuts.append(pNGCerNpeSum_cut)      #5     
            
        elif self.eff_type == 'HMS':
            hgoodScinHit_cut = C.VCUT(1.,'hgoodScinHit_cut')
            hBetanotrk_cut = C.WCUT(0.5,1.5,'hBetanotrk_cut')
            hCaletotnorm_cut = C.WCUT(0.01,0.6,'hCaletotnorm_cut')
            hDCntrack_cut = C.WCUT(1.,np.inf,'hDCntrack')
            hgoodstarttime_cut = C.VCUT(1.,'hgoodstarttime_cut')
            hCerNpeSum_cut = C.WCUT(-np.inf,0.5,'hCerNpeSum')
            
            self.cuts.append(hgoodScinHit_cut)      #0
            self.cuts.append(hBetanotrk_cut)        #1
            self.cuts.append(hCaletotnorm_cut)      #2
            self.cuts.append(hDCntrack_cut)         #3
            self.cuts.append(hgoodstarttime_cut)    #4
            self.cuts.append(hCerNpeSum_cut)        #5
            
        elif self.eff_type == 'LT':
            TRIG_OFF_cut = C.VCUT(0,name='TRIG_OFF_cut')
            
            self.cuts.append(TRIG_OFF_cut)  #0
            
        elif self.eff_type == 'cpuLT':
            TRIG_OFF_cut = C.VCUT(0,name='TRIG_OFF_cut')
            
            self.cuts.append(TRIG_OFF_cut)  #0
        else:
            print('No efficiency type given.')


        #initialize current cut
        self.curr_cut = C.current_cut(data_obj,current=self.curr,
                                      many=self.many)



      
        

#%%

# T_sel = ['H.kin.secondary.emiss','H.kin.secondary.pmiss',
#          'H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph','P.gtr.th',
#          'P.gtr.ph','H.extcor.xsieve','H.extcor.ysieve',
#          'P.extcor.xsieve','P.extcor.ysieve','H.react.z','P.react.z',
#          'CTime.epCoinTime_ROC2','H.hod.goodscinhit','H.hod.betanotrack',
#          'H.cal.etotnorm','H.cer.npeSum','H.dc.ntrack','P.hod.goodscinhit',
#          'P.hod.betanotrack','P.cal.etotnorm','P.hgcer.npeSum',
#          'P.ngcer.npeSum','P.dc.ntrack','T.coin.pTRIG1_ROC2_tdcTimeRaw',
#          'T.coin.pTRIG2_ROC2_tdcTimeRaw','T.coin.pTRIG6_ROC2_tdcTimeRaw',
#          'T.coin.pEDTM_tdcTimeRaw']

# TSP_sel = ['evNumber','P.BCM4A.scalerCharge','P.BCM4A.scalerChargeCut',
#            'P.BCM4A.scalerCurrent','P.BCM4B.scalerChargeCut',
#            'P.BCM4B.scalerCurrent','P.BCM4C.scalerChargeCut',
#            'P.BCM4C.scalerCurrent','P.pTRIG1.scalerCut','P.pTRIG2.scalerCut',
#            'P.pTRIG6.scalerCut','P.1MHz.scalerTimeCut','P.EDTM.scalerCut',
#            'P.S1X.scalerCut']

# TSH_sel = ['H.S1X.scalerCut','H.1MHz.scalerTimeCut']

T_sel = ['H.hod.goodscinhit','H.hod.betanotrack',
         'H.cal.etotnorm','H.cer.npeSum','H.dc.ntrack','P.hod.goodscinhit',
         'P.hod.betanotrack','P.cal.etotnorm','P.hgcer.npeSum',
         'P.ngcer.npeSum','P.dc.ntrack','T.coin.pTRIG1_ROC2_tdcTimeRaw',
         'T.coin.pTRIG2_ROC2_tdcTimeRaw','T.coin.pTRIG6_ROC2_tdcTimeRaw',
         'T.coin.pEDTM_tdcTimeRaw','H.hod.goodstarttime','P.hod.goodstarttime']

TSP_sel = ['evNumber','P.BCM1.scalerCurrent','P.BCM2.scalerCurrent',
           'P.BCM4A.scalerCurrent','P.BCM4B.scalerCurrent',
           'P.BCM4C.scalerCurrent','P.BCM1.scalerCharge',
           'P.BCM2.scalerCharge','P.BCM4A.scalerCharge',
           'P.BCM4B.scalerCharge','P.BCM4C.scalerCharge',
           'P.pTRIG1.scalerCut',
           'P.pTRIG2.scalerCut',
           'P.pTRIG6.scalerCut','P.1MHz.scalerTime','P.EDTM.scalerCut',
           'P.S1X.scalerCut']

TSH_sel = ['H.S1X.scalerCut','H.1MHz.scalerTimeCut']

# R = [20840,20841,20846,20851,20858,20861,20868,20869] #heep_coin
# R = [20871,20872]
# ,21065,21066,21067,21068,21069,21070, 
#      21071,21072,21073,21074,21075] #pm_120
# R = [20873,20874,20875,20876,20877,20878,20880,20881,20882,20883,21076]
# R = [21049,21051,21052,21053,21054,21055,21056,21057,21058,21059,21060,
#      21061,21062,21063,21064]

#      21078,21079,21080,21081,21082,21083,21084,21085,21086,21087,21088,
#      21089,21090,21091,21092,21093,21094,21095,21096,21097,21098,21099,
#      21100,21101,21102] # pm_580
# R = [20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,20896,
#      20897,20898,20899,20900,20901,20902,20903,20904,20905,20907,20908,
#      20909,20910,20911,20912,20913,20914,20915,20916,20921,20922,20923,
#      20924,20925,20926,20927,20928,20929,20930,20931,20932,20933,20934,
#      20935,20936,20937,20938,20939,20940,20941,20942,20943,20944,20945,
#      20949,20950,20951,20953,20954,20955,20956] # pm_800

# all runs:
# R = [20840,20841,20846,20851,20858,20861,20868,20869,20871,20872,
#      20873,20874,20875,20876,20877,20878,20880,20881,20882,20883,
#      20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,
#      20896,20897,20898,20899,20900,20901,20902,20903,20904,20905,
#      20907,20908,20909,20910,20911,20912,20913,20914,20915,20916,
#      20921,20922,20923,20924,20925,20926,20927,20928,20929,20930,
#      20931,20932,20933,20934,20935,20936,20937,20938,20939,20940,
#      20941,20942,20943,20944,20945,20949,20950,20951,20953,20954]

# R = [20955,20956,20958,20959,20960,20961,20962,20963,20965,20966,
#      20969,20970,20971,20972,20973,20974,20975,20976,20977,20978,
#      20979,20980,20981,20982,20983,20984,20985,20986,20987,20988,
#      20989,20990,20991,20992,20993,20994,20995,20996,20997,20998,
#      20999,21000,21001,21002,21003,21004,21005,21006,21007,21008,
#      21009,21011,21012,21013,21014,21015,21016,21017,21018,21019,
#      21020,21021,21022,21023,21024,21025,21026,21027,21028,21029,
#      21030,21031,21032,21033,21034,21036,21037,21038,21039,21040]

R = [21041,21042,21043,21044,21045,21046,21047,21048,21065,21066,
     21067,21068,21069,21070,21071,21072,21073,21074,21075,21076,
     21077,21078,21079,21080,21081,21082,21083,21084,21085,21086,
     21087,21088,21089,21090,21091,21092,21093,21094,21095,21096,
     21097,21098,21099,21100,21101,21102]

# 20843,20844,20845,20847,
#      20850,20856,20857,20859,20860,20862,20863,20864,20865,20866,
#      20867,20870,21049,21050,21051,21052,21053,21054,21055,21056,
#      21057,21058,21059,21060,21061,21062,21063,21064]


# R = 21063

rdir = "/media/gvill/Gema's T7/ROOTfiles/pass_3/"

RUN = D.DATA_INIT(data_type='deut23_data', run= R, 
                  select_branches= {'TSP':TSP_sel},
                  select_trees=['TSP'],
                  ROOTfiles_path=rdir)
#%% testing current cut
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

#%% calc hms/shms tracking efficiencies - test
pgoodScinHit_cut = C.VCUT(1.,'pgoodScinHit_cut')
pBetanotrk_cut = C.WCUT(0.5,1.5,'pBetanotrk_cut')
pCaletotnorm_cut = C.WCUT(0.6,np.inf,'pCaletotnorm_cut')
# pNGCer_npeSum_cut = C.WCUT(0.5,np.inf,'pNGCer_npeSum_cut')
pDCntrack_cut = C.WCUT(0.5,np.inf,'pDCntrack')

hgoodStartTime_cut = C.VCUT(1.,'hgoodStartTime_cut')
hgoodScinHit_cut = C.VCUT(1.,'hgoodScinHit_cut')
hBetanotrk_cut = C.WCUT(0.5,1.5,'hBetanotrk_cut')
hCaletotnorm_cut = C.WCUT(0.01,0.6,'hCaletotnorm_cut')
hCer_npeSum_cut = C.WCUT(-np.inf,0.5,'hCer_npeSum_cut')
hDCntrack_cut = C.WCUT(1.0,np.inf,'hDCntrack')

#apply current cut
curr_cut = C.current_cut(RUN,current='BCM4A',many=True)
curr_cut_arrays = curr_cut()

#EDTM cut
cuts_list = [C.noEDTM]

cuts_to_apply = {}
for m in RUN.many:
    clist = []
    for cut in cuts_list:
        br = RUN.Branches[m][C.HCANA_names[cut.name]]
        cut_array = cut(br)
        # cut.stats()        
        clist.append(cut_array)
    clist.append(curr_cut_arrays[m])
    cuts_to_apply[m] = clist
    
all_cuts = {}
for m in RUN.many:
    all_cuts[m] = cuts_to_apply[m][0]
    for arr in cuts_to_apply[m]:    
        all_cuts[m] = all_cuts[m] & arr     
    
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
    
    # print('ptrk_eff = ', ptrk_eff, '$\pm$', ptrk_eff_err)
    ptrk.append((ptrk_eff,ptrk_eff_err))
    
    hgoodStartTime = RUN.Branches[m]['H.hod.goodstarttime']
    hgoodScinHit = RUN.Branches[m]['H.hod.goodscinhit']
    hBetanotrk = RUN.Branches[m]['H.hod.betanotrack']
    hCaletotnorm = RUN.Branches[m]['H.cal.etotnorm']
    hCer_npeSum = RUN.Branches[m]['H.cer.npeSum']
    hDCntrack = RUN.Branches[m]['H.dc.ntrack']
    
    hgoodStartTime = hgoodStartTime[all_cuts[m]]
    hgoodScinHit = hgoodScinHit[all_cuts[m]]
    hBetanotrk = hBetanotrk[all_cuts[m]]
    hCaletotnorm = hCaletotnorm[all_cuts[m]]
    hCer_npeSum = hCer_npeSum [all_cuts[m]]
    hDCntrack = hDCntrack[all_cuts[m]]
    
    n_hgoodStartTime = hgoodStartTime_cut(hgoodStartTime)
    n_hgoodScinHit = hgoodScinHit_cut(hgoodScinHit)
    n_hBetanotrk = hBetanotrk_cut(hBetanotrk)
    n_hCaletotnorm = hCaletotnorm_cut(hCaletotnorm)
    n_hDCntrack = hDCntrack_cut(hDCntrack)
    n_hCerNpeSum = hCer_npeSum_cut(hCer_npeSum) 
    
    # print('n_hgoodScinHit = ',n_hgoodScinHit.sum())
    # print('n_hBetanotrk = ',n_hBetanotrk.sum())
    # print('n_hCaletotnorm = ',n_hCaletotnorm.sum())
    # print('n_hDCntrack = ',n_hDCntrack.sum())
    # print('n_hCerNpeSum = ',n_hCerNpeSum.sum())
    # print('n_hgoodStartTime = ', n_hgoodStartTime.sum())
    
    good_had_should_list = [n_hgoodScinHit,n_hBetanotrk,n_hCaletotnorm,
                            n_hgoodStartTime]
    # good_had_should_list = [n_hgoodScinHit,n_hBetanotrk,n_hgoodStartTime]
    # good_had_should_list = [n_hCaletotnorm]
    
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
#%% test
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

#%% efficiency plots 
def pptxify(t='',x='',y='',fsize=24): 
    B.pl.tick_params('both',labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize':fsize})
    B.pl.xlabel(x, fontdict={'fontsize':fsize})
    B.pl.ylabel(y, fontdict={'fontsize':fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(10,10)
    #B.pl.legend(fontsize = fsize)

# h_pgoodScinHit = B.histo(pgoodScinHit,range=(-0.5,2.5),bins=3)
# B.pl.figure()
# h_pgoodScinHit.plot()
# ylo, yhi = B.pl.ylim()
# B.pl.vlines(1.0, 0., yhi,color='black',linestyles='--')
# pptxify('pgoodScinHit')

# h_pBetanotrk = B.histo(pBetanotrk,range=(0.1,2.0),bins=100)
# B.pl.figure()
# h_pBetanotrk.plot()
# ylo, yhi = B.pl.ylim()
# B.pl.vlines([0.5,1.5], 0., yhi,color='black',linestyles='--')
# pptxify('pBetanotrk')

# h_pCaletotnorm = B.histo(pCaletotnorm,range=(0.5,1.5),bins=100)
# B.pl.figure()
# h_pCaletotnorm.plot()
# ylo, yhi = B.pl.ylim()
# B.pl.vlines(0.6, 0., yhi,color='black',linestyles='--')
# pptxify('pCaletotnorm')

# h_pDCntrack = B.histo(pDCntrack,range=(-0.5,5.5),bins=6)
# B.pl.figure()
# h_pDCntrack.plot()
# ylo, yhi = B.pl.ylim()
# B.pl.vlines(0.5, 0., yhi,color='black',linestyles='--')
# pptxify('pDCntrack')

h_hgoodScinHit = B.histo(hgoodScinHit,range=(-0.5,2.5),bins=3)
# B.pl.figure()
# h_hgoodScinHit.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(1.0, 0., yhi,color='black',linestyles='--')
pptxify('hgoodScinHit')

h_hBetanotrk = B.histo(hBetanotrk,range=(0.5,1.5),bins=100)
# B.pl.figure()
# h_hBetanotrk.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines([0.5,1.5], 0., yhi,color='black',linestyles='--')
pptxify('hBetanotrk')

h_hCaletotnorm = B.histo(hCaletotnorm,range=(0.,1.),bins=100)
# B.pl.figure()
# h_hCaletotnorm.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines([0.,0.6], 0., yhi,color='black',linestyles='--')
pptxify('hCaletotnorm')

h_hDCntrack = B.histo(hDCntrack,range=(-0.5,5.5),bins=6)
# B.pl.figure()
# h_hDCntrack.plot()
ylo, yhi = B.pl.ylim()
B.pl.vlines(0.5, 0., yhi,color='black',linestyles='--')
pptxify('hDCntrack')

h_hCerNpeSum = B.histo(hCer_npeSum,range=(0.5,2.0),bins=100)
B.pl.figure()
h_hCerNpeSum.plot()

#%% calculate live time - test

pEDTM_tdcTimeRaw = RUN.Branches['T.coin.pEDTM_tdcTimeRaw']
pTRIG6_tdcTimeRaw = RUN.Branches['T.coin.pTRIG6_ROC2_tdcTimeRaw']
pTRIG2_tdcTimeRaw = RUN.Branches['T.coin.pTRIG2_ROC2_tdcTimeRaw']
pEDTM_scaler = RUN.Branches['P.EDTM.scalerCut']
pTRIG6_scaler = RUN.Branches['P.pTRIG6.scalerCut']
pTRIG2_scaler = RUN.Branches['P.pTRIG2.scalerCut']
scalerTime = RUN.Branches['P.pTRIG6.scalerCut']

pTRIG6_scalerRateRaw = pTRIG6_scaler/scalerTime
pEDTM_scalerRate = pEDTM_scaler/scalerTime

pTRIG6_scalerRate = pTRIG6_scalerRateRaw - pEDTM_scalerRate

TRIG_OFF_cut = C.VCUT(0,name='TRIG_OFF_cut')

noEDTMcut_array = TRIG_OFF_cut(pEDTM_tdcTimeRaw)
EDTMcut_array = ~TRIG_OFF_cut(pEDTM_tdcTimeRaw)
pTRIG6cut_array = ~TRIG_OFF_cut(pTRIG6_tdcTimeRaw)
pTRIG2cut_array = ~TRIG_OFF_cut(pTRIG2_tdcTimeRaw)

#apply current cut
curr_cut = C.current_cut(RUN,current='BCM4A')
curr_cut_arrays = curr_cut()

noEDTMcut_array = noEDTMcut_array[curr_cut_arrays]
EDTMcut_array = EDTMcut_array[curr_cut_arrays]
pTRIG6cut_array = pTRIG6cut_array[curr_cut_arrays]
pTRIG2cut_array = pTRIG2cut_array[curr_cut_arrays]

pEDTM_scalerRate = pEDTM_scalerRate[curr_cut.cut_array]
pTRIG6_scalerRateRaw = pTRIG6_scalerRateRaw[curr_cut.cut_array]
pTRIG6_scalerRate = pTRIG6_scalerRate[curr_cut.cut_array]

tLT_corr_factor = 1 - pTRIG6_scalerRateRaw[-1]*250e-06 +\
    pTRIG6_scalerRate[-1]*250e-06*(1 + pEDTM_scalerRate[-1]/pTRIG6_scalerRateRaw[-1])

tLT = EDTMcut_array.sum()/pEDTM_scaler[curr_cut.cut_array][-1]

tLT_corr = tLT*tLT_corr_factor
  
pTRIG6_noEDTM = noEDTMcut_array & pTRIG6cut_array
pTRIG2_noEDTM = noEDTMcut_array & pTRIG2cut_array

cpuLT = pTRIG6_noEDTM.sum()/\
                (pTRIG6_scaler[curr_cut.cut_array][-1] -\
                         pEDTM_scaler[curr_cut.cut_array][-1])  

T2cpuLT = pTRIG2_noEDTM.sum()/\
                (pTRIG2_scaler[curr_cut.cut_array][-1] -\
                         pEDTM_scaler[curr_cut.cut_array][-1]) 
print(T2cpuLT)
#%% write charge using new cut
deutDB_name= 'deuteron_db.db'

BCM1_cut = C.current_cut(RUN,current='BCM1',many=True)
BCM2_cut = C.current_cut(RUN,current='BCM2',many=True)
BCM4A_cut = C.current_cut(RUN,current='BCM4A',many=True)
BCM4B_cut = C.current_cut(RUN,current='BCM4B',many=True)
BCM4C_cut = C.current_cut(RUN,current='BCM4C',many=True)

BCM1_cut()
BCM2_cut()
BCM4A_cut()
BCM4B_cut()
BCM4C_cut()

for r in R:

    # get charge with current cut
    charge_BCM1_cut = BCM1_cut.scaler(RUN.Branches[r]['P.BCM1.scalerCharge'],BCM1_cut.cut_array[r])
    charge_BCM2_cut = BCM2_cut.scaler(RUN.Branches[r]['P.BCM2.scalerCharge'],BCM2_cut.cut_array[r])
    charge_BCM4A_cut = BCM4A_cut.scaler(RUN.Branches[r]['P.BCM4A.scalerCharge'],BCM4A_cut.cut_array[r])
    charge_BCM4B_cut = BCM4B_cut.scaler(RUN.Branches[r]['P.BCM4B.scalerCharge'],BCM4B_cut.cut_array[r])
    charge_BCM4C_cut = BCM4C_cut.scaler(RUN.Branches[r]['P.BCM4C.scalerCharge'],BCM4C_cut.cut_array[r])
    
    # get scl time with current cut
    sclTime_BCM1_cut = BCM1_cut.scaler(RUN.Branches[r]['P.1MHz.scalerTime'],BCM1_cut.cut_array[r])
    sclTime_BCM2_cut = BCM2_cut.scaler(RUN.Branches[r]['P.1MHz.scalerTime'],BCM2_cut.cut_array[r])
    sclTime_BCM4A_cut = BCM4A_cut.scaler(RUN.Branches[r]['P.1MHz.scalerTime'],BCM4A_cut.cut_array[r])
    sclTime_BCM4B_cut = BCM4B_cut.scaler(RUN.Branches[r]['P.1MHz.scalerTime'],BCM4B_cut.cut_array[r])
    sclTime_BCM4C_cut = BCM4C_cut.scaler(RUN.Branches[r]['P.1MHz.scalerTime'],BCM4C_cut.cut_array[r])
    
    current_BCM1 = charge_BCM1_cut/sclTime_BCM1_cut
    current_BCM2 = charge_BCM2_cut/sclTime_BCM2_cut
    current_BCM4A = charge_BCM4A_cut/sclTime_BCM4A_cut
    current_BCM4B = charge_BCM4B_cut/sclTime_BCM4B_cut
    current_BCM4C = charge_BCM4C_cut/sclTime_BCM4C_cut
    
    col_names = ['BCM1_charge','BCM1_current','BCM2_charge','BCM2_current',
                 'BCM4A_charge','BCM4A_current','BCM4B_charge','BCM4B_current',
                 'BCM4C_charge','BCM4C_current']
    col_vals = [f'{charge_BCM1_cut/1e3:.3f}',f'{current_BCM1:.3f}',
                f'{charge_BCM2_cut/1e3:.3f}',f'{current_BCM2:.3f}',
                f'{charge_BCM4A_cut/1e3:.3f}',f'{current_BCM4A:.3f}',
                f'{charge_BCM4B_cut/1e3:.3f}',f'{current_BCM4B:.3f}',
                f'{charge_BCM4C_cut/1e3:.3f}',f'{current_BCM4C:.3f}']
    
    # print(r,[f'{nn} = {col_vals[i]}' for i,nn in enumerate(col_names)])
    db.update_row(deutDB_name, 'RUN_LIST_UPDATED',col_names, col_vals, 
                  where=f'run = {r}')
    
#%% write charge and current values to db

deutDB_name= 'deuteron_db.db'

# hms_trk_eff = get_spec_eff(RUN,eff_type='HMS',many=True)
# htrk = hms_trk_eff.calc_eff()
# htrk_err = hms_trk_eff.trk_eff_err

# shms_trk_eff = get_spec_eff(RUN,eff_type='SHMS',many=True)
# ptrk= shms_trk_eff.calc_eff()
# ptrk_err = shms_trk_eff.trk_eff_err

# LT_eff = get_spec_eff(RUN,eff_type='LT',many=True)
# LT= LT_eff.calc_eff()
 
run =  RUN.many
for r in run:
    BCM1_charge = RUN.Branches[r]['P.BCM1.scalerChargeCut'][-1]/1e3
    BCM1_current = RUN.Branches[r]['P.BCM1.scalerChargeCut'][-1]/\
        RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]
    
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
    
    ht = 0.0
    ht_e = 0.0
    # pt = ptrk[r]
    # pt_e = ptrk_err[r]
    # lt = LT[r] 

    ps1x = RUN.Branches[r]['P.S1X.scalerCut'][-1]
    pstime = RUN.Branches[r]['P.1MHz.scalerTimeCut'][-1]
    
    hs1x = 0.0
    hstime = 1.0
    
    ps1xrate = ps1x/pstime/1000. #kHz
    hs1xrate = hs1x/hstime/1000. #kHz       
    
    # col_names = ['BCM4A_charge','BCM4A_current','BCM4B_charge','BCM4B_current',
    #              'BCM4C_charge','BCM4C_current','T1_scl_rates','T2_scl_rates',
    #              'T6_scl_rates','tLT','HMS_TrkEff','HMS_TrkEff_Err',
    #              'SHMS_TrkEff','SHMS_TrkEff_Err','pS1X_scl_rates','hS1X_scl_rates']
    # col_vals = [f'{BCM4A_charge:.3f}', f'{BCM4A_current:.3f}',
    #             f'{BCM4B_charge:.3f}', f'{BCM4B_current:.3f}',
    #             f'{BCM4C_charge:.3f}', f'{BCM4C_current:.3f}',
    #             f'{T1_scl_rate:.3f}',f'{T2_scl_rate:.3f}',
    #             f'{T6_scl_rate:.3f}',f'{lt:.4f}',f'{ht:.4f}',f'{ht_e:.4f}',
    #             f'{pt:.4f}',f'{pt_e:.4f}',f'{ps1xrate:.3f}', f'{hs1xrate:.3f}']
    col_names = ['BCM1_charge','BCM1_current']
    col_vals = [f'{BCM1_charge:.3f}',f'{BCM1_current:.3f}']
    
    # print(r,[f'{nn} = {col_vals[i]}' for i,nn in enumerate(col_names)])
    db.update_row(deutDB_name, 'RUN_LIST_UPDATED',col_names, col_vals, 
                  where=f'run = {r}')
#%%   
deutDB_name= 'deuteron_db.db'

cpuLT_eff = get_spec_eff(RUN,eff_type='cpuLT',many=True)
cpuLT = cpuLT_eff.calc_eff()


for m in RUN.many:
    T1 = cpuLT[m][0][0]
    T2 = cpuLT[m][1][0]
    T6 = cpuLT[m][2][0]
    
    col_names = ['T1_cpuLT','T2_cpuLT','T6_cpuLT']
    col_vals = [f'{T1:.4f}', f'{T2:.4f}', f'{T6:.4f}']
    
    # print([f'{nn} = {col_vals[i]}' for i,nn in enumerate(col_names)])
    db.update_row(deutDB_name, 'RUN_LIST_UPDATED',col_names, col_vals, 
                  where=f'run = {m}') 
 

#%%               
R = [20840,20841,20846,20851,20858,20861,20868,20869,20871,20872,
     20873,20874,20875,20876,20877,20878,20880,20881,20882,20883,
     20886,20887,20888,20889,20890,20891,20892,20893,20894,20895,
     20896,20897,20898,20899,20900,20901,20902,20903,20904,20905,
     20907,20908,20909,20910,20911,20912,20913,20914,20915,20916,
     20921,20922,20923,20924,20925,20926,20927,20928,20929,20930,
     20931,20932,20933,20934,20935,20936,20937,20938,20939,20940,
     20941,20942,20943,20944,20945,20949,20950,20951,20953,20954,
     20955,20956,20958,20959,20960,20961,20962,20963,20965,20966,
     20969,20970,20971,20972,20973,20974,20975,20976,20977,20978,
     20979,20980,20981,20982,20983,20984,20985,20986,20987,20988,
     20989,20990,20991,20992,20993,20994,20995,20996,20997,20998,
     20999,21000,21001,21002,21003,21004,21005,21006,21007,21008,
     21009,21011,21012,21013,21014,21015,21016,21017,21018,21019,
     21020,21021,21022,21023,21024,21025,21026,21027,21028,21029,
     21030,21031,21032,21033,21034,21036,21037,21038,21039,21040,
     21041,21042,21043,21044,21045,21046,21047,21048,21065,21066,
     21067,21068,21069,21070,21071,21072,21073,21074,21075,21076,
     21077,21078,21079,21080,21081,21082,21083,21084,21085,21086,
     21087,21088,21089,21090,21091,21092,21093,21094,21095,21096,
     21097,21098,21099,21100,21101,21102,20843,20844,20845,20847,
     20850,20856,20857,20859,20860,20862,20863,20864,20865,20866,
     20867,20870,21049,21050,21051,21052,21053,21054,21055,21056,
     21057,21058,21059,21060,21061,21062,21063,21064]

hte = []
htee = []
pte = []
ptee = []
lte = []
ps1x = []
hs1x = []
bcm1c=[] 
bcm1q=[]
bcm4ac=[]
bcm4aq=[]
bcm4bc=[]
bcm4bq=[]
bcm4cc=[]
bcm4cq=[]
bcm2c=[]
bcm2q=[]
t6lt = []
t2lt = []
t1lt = []
for run in R:
    h_teff, h_teff_err, p_teff, p_teff_err, lt, ps1x_scl, hs1x_scl,\
    bcm1_c,bcm1_q,bcm4a_c,bcm4a_q,bcm4b_c,bcm4b_q,bcm4c_c,bcm4c_q,bcm2_c,bcm2_q,\
        t6_lt,t2_lt,t1_lt=\
        db.retrieve('deuteron_db.db', 
                    'HMS_TrkEff, HMS_TrkEff_Err, SHMS_TrkEff,'+\
                    ' SHMS_TrkEff_Err, tLT, pS1X_scl_rates, hS1X_scl_rates,'+\
                    ' BCM1_current, BCM1_charge, BCM4A_current, BCM4A_charge,'+\
                    ' BCM4B_current, BCM4B_charge, BCM4C_current, BCM4C_charge,'+\
                    ' BCM2_current, BCM2_charge,T6_cpuLT,T2_cpuLT,T1_cpuLT', 
                    'RUN_LIST_UPDATED', where = f"run=\'{run}\'")[0]
        
    hte.append(h_teff)
    htee.append(h_teff_err)
    pte.append(p_teff)
    ptee.append(p_teff_err)
    lte.append(lt)
    ps1x.append(ps1x_scl)
    hs1x.append(hs1x_scl)
    bcm1c.append(bcm1_c) 
    bcm1q.append(bcm1_q)
    bcm4ac.append(bcm4a_c)
    bcm4aq.append(bcm4a_q)
    bcm4bc.append(bcm4b_c)
    bcm4bq.append(bcm4b_q)
    bcm4cc.append(bcm4c_c)
    bcm4cq.append(bcm4c_q)
    bcm2c.append(bcm2_c) 
    bcm2q.append(bcm2_q)
    t6lt.append(t6_lt)
    t2lt.append(t2_lt)
    t1lt.append(t1_lt)

#%%
R = np.array(R)
hte = np.array(hte)
htee = np.array(htee)
lte = np.array(lte)

B.pl.figure()
B.pl.subplot(2,1,1)
hte_R = B.plot_exp(R[hte >0.],hte[hte >0.],htee[hte >0.],marker='^',label='HMS')
pte_R = B.plot_exp(R,pte,ptee,marker='v',label='SHMS')

hmean = np.mean(hte[hte >0.])
pmean = np.mean(pte)

B.pl.hlines([hmean,pmean,1.0],xmin=20825,xmax=21115,
            color=[cmp['Set1'](1),cmp['Set1'](4),'red'],linestyle='--')
B.pl.text(21100, hmean-0.014, f'{hmean:0.4f}',color=cmp['Set1'](1))
B.pl.text(20830, pmean-0.014, f'{pmean:0.4f}',color=cmp['Set1'](4))

B.pl.xlabel('Run #')
B.pl.ylabel('Tracking Efficiency')
B.pl.legend()
B.pl.ylim(0.90,1.05)

B.pl.subplot(2,1,2)
lte_R = B.plot_exp(R[lte>0.5],lte[lte>0.5],marker='*')
ltmean = np.mean(lte[lte>0.5])
B.pl.hlines([ltmean,1.0],xmin=20825,xmax=21115,
            color=[cmp['Set1'](4),'red'],linestyle='--')
B.pl.text(20830, ltmean-0.014, f'{ltmean:0.4f}',color=cmp['Set1'](4))
B.pl.xlabel('Run #')
B.pl.ylabel('Live Time')
B.pl.ylim(0.90,1.05)
#%%
t6lt = np.array(t6lt)
t2lt = np.array(t2lt)
t1lt = np.array(t1lt)

B.pl.figure()
# B.pl.subplot(3,1,1)
t6_R = B.plot_exp(R[t6lt>0],t6lt[t6lt>0],marker='*',label='T6')
t2_R = B.plot_exp(R[t2lt>0],t2lt[t2lt>0],marker='*',label='T2')
t1_R = B.plot_exp(R[t1lt>0],t1lt[t1lt>0],marker='*',label='T1')

B.pl.figure()
# B.pl.subplot(3,1,1)
t6_R = B.plot_exp(ps1x[t6lt>0],t6lt[t6lt>0],marker='*',label='T6')
t2_R = B.plot_exp(ps1x[t2lt>0],t2lt[t2lt>0],marker='*',label='T2')
t1_R = B.plot_exp(ps1x[t1lt>0],t1lt[t1lt>0],marker='*',label='T1')


#%%
ps1x = np.array(ps1x)
hs1x = np.array(hs1x)

B.pl.figure()
B.pl.subplot(1,2,1)
hte_R = B.plot_exp(hs1x[hte >0.],hte[hte >0.],htee[hte >0.],marker='^',label='HMS')
B.linefit(hs1x[hte >0.] ,hte[hte >0.])
B.pl.xlabel('HMS S1X Scaler Rate [kHz]')
B.pl.ylabel('Tracking Efficiency')

# B.pl.figure()
B.pl.subplot(1,2,2)
pte_R = B.plot_exp(ps1x,pte,ptee,marker='v',label='SHMS')

B.linefit(ps1x ,pte)

B.pl.xlabel('SHMS S1X Scaler Rate [kHz]')
B.pl.ylabel('Tracking Efficiency')

#%%
B.pl.figure()
hs1x_R = B.plot_exp(R,hs1x,marker='^',label='HMS')
B.pl.title('HMS S1X rates vs. Run')

B.pl.figure()
ps1x_R = B.plot_exp(R,ps1x,marker='v',label='SHMS')
B.pl.title('SHMS S1X rates vs. Run')

#%%
R = np.array(R)
bcm1q = np.array(bcm1q)
bcm2q = np.array(bcm2q)
bcm4aq = np.array(bcm4aq)
bcm4cq = np.array(bcm4cq)
bcm1c = np.array(bcm1c)
bcm2c = np.array(bcm2c)
bcm4ac = np.array(bcm4ac)
bcm4cc = np.array(bcm4cc)

B.pl.figure()
B.plot_exp(R,bcm1c,marker='^',color=cmp['Set1'](0))
B.plot_exp(R,bcm4ac,marker='^',color=cmp['Set1'](1))
B.plot_exp(R,bcm4bc,marker='^',color=cmp['Set1'](2))
B.plot_exp(R,bcm4cc,marker='^',color=cmp['Set1'](3))

B.pl.figure()
B.plot_exp(R,bcm1q,marker='^',color=cmp['Set1'](0))
B.plot_exp(R,bcm4aq,marker='^',color=cmp['Set1'](1))
B.plot_exp(R,bcm4bq,marker='^',color=cmp['Set1'](2))

#BCM1 ratios
B.pl.figure()
B.plot_exp(R,bcm1q/bcm2q,marker='^',color=cmp['Set1'](0),label='BCM1/BCM2')
B.plot_exp(R,bcm1q/bcm4aq,marker='^',color=cmp['Set1'](1),label='BCM1/BCM4A')
B.plot_exp(R,bcm1q/bcm4cq,marker='^',color=cmp['Set1'](2),label='BCM1/BCM4C')
B.pl.title('BCM1 Ratios')
B.pl.legend()