#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 10:26:33 2025

@author: gvill

Target Boiling Studies for SHMS

Determining yields as a function of current for the luminosity scan runs.

"""

import data_init as D
import database_operations as db
import cut_handler as C
import LT.box as B

import numpy as np
# import matplotlib.lines as mlines
from matplotlib import colormaps as cmp

#Set default font to times new roman
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14
}
B.pl.rc('font', **font)

#%% runs for the different targets

C12_runs = [21059,21060,21061,21062,21063] # 15,30,40,55,70 muA

LD2_runs = [21049,21050,21051,21052,21053,  # 10,15,25,35,40
            21054,21055,21056,21057,21058]  # 40,45,55,65,70 muA

#%% hcana variables needed for this study

TSP_branches = ['evNumber','P.BCM1.scalerCurrent','P.BCM2.scalerCurrent',
                'P.BCM4A.scalerCurrent',
                'P.BCM4B.scalerCurrent','P.BCM4C.scalerCurrent',
                'P.1MHz.scalerTime',
                'P.S1X.scaler','P.S1X.scalerCut','P.S1XS1Y.scaler',
                'P.EDTM.scaler','P.pTRIG1.scaler','P.pTRIG2.scaler',
                'P.pTRIG3.scalerCut',
                'P.pTRIG6.scaler','P.pTRIG1.scaler','P.pTRIG2.scaler',
                'P.EDTM.scalerCut',
                'P.pEL_REAL.scaler','P.pEL_REAL.scalerCut',
                'P.BCM4A.scalerCharge','P.BCM4B.scalerCharge',
                'P.BCM4C.scalerCharge','P.BCM1.scalerCharge',
                'P.BCM2.scalerCharge']

T_branches = ['P.ngcer.npeSum','P.hgcer.npeSum','P.cal.etotnorm',
              'P.cal.etottracknorm','P.gtr.dp','P.kin.primary.Q2',
              'H.react.z','P.react.z','P.extcor.xsieve','P.extcor.ysieve',
              'T.coin.pEDTM_tdcTimeRaw','P.hod.goodscinhit','P.hod.betanotrack',
              'P.dc.ntrack','P.hod.goodstarttime','P.ngcer.npeSum',
              'T.coin.pTRIG1_ROC2_tdcTimeRaw',
              'T.coin.pTRIG2_ROC2_tdcTimeRaw',
              'T.coin.pTRIG6_ROC2_tdcTimeRaw']

#%% load rootfiles
ROOTfiles_DIR = "/media/gvill/Gema's T7/ROOTfiles/xem2BCMcalib/"

# C12_data = D.DATA_INIT(data_type='deut23_data',run=C12_runs,
#                          select_branches={'T':T_branches,'TSP':TSP_branches},
#                          select_trees=['T','TSP'],
#                          ROOTfiles_path=ROOTfiles_DIR) 

LD2_data = D.DATA_INIT(data_type='deut23_data',run=LD2_runs,
                         select_branches={'T':T_branches,'TSP':TSP_branches},
                         select_trees=['T','TSP'],
                         ROOTfiles_path=ROOTfiles_DIR) 
#%% calculate trk eff and LT using class
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
        
        if eff_type == 'LT':
            if many:
                for m in data_obj.many:
                    self.var[m][0] = self.var[m][0][cut[m]]
                    self.var[m][1] = self.curr_cut.scaler(self.var[m][1],self.curr_cut.cut_array[m])
                    self.var[m][2] = self.curr_cut.scaler(self.var[m][2],self.curr_cut.cut_array[m])
                    self.var[m][3] = self.curr_cut.scaler(self.var[m][3],self.curr_cut.cut_array[m])
            else:
                self.var[0] = self.var[0][cut]
                self.var[1] = self.curr_cut.scaler(self.var[1],self.curr_cut.cut_array)
                self.var[2] = self.curr_cut.scaler(self.var[2],self.curr_cut.cut_array)
                self.var[3] = self.curr_cut.scaler(self.var[3],self.curr_cut.cut_array)
                
        elif eff_type == 'cpuLT':
            if many:
                for m in data_obj.many:
                    self.var[m][0] = self.var[m][0][cut[m]]
                    self.var[m][1] = self.var[m][1][cut[m]]
                    self.var[m][2] = self.var[m][2][cut[m]]
                    self.var[m][3] = self.var[m][3][cut[m]]
                    self.var[m][4] = self.curr_cut.scaler(self.var[m][4],self.curr_cut.cut_array[m])
                    self.var[m][5] = self.curr_cut.scaler(self.var[m][5],self.curr_cut.cut_array[m])
                    self.var[m][6] = self.curr_cut.scaler(self.var[m][6],self.curr_cut.cut_array[m])
                    self.var[m][7] = self.curr_cut.scaler(self.var[m][7],self.curr_cut.cut_array[m])
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

                    # tLT_corr_factor = 1 - (pTRIG6_scalerRate[-1] +\
                    #                        pEDTM_scalerRate[-1])*250e-06 +\
                    #                         pTRIG6_scalerRate[-1]*250e-06*\
                    #                         (1 + pEDTM_scalerRate[-1]/\
                    #                          (pTRIG6_scalerRate[-1]+\
                    #                           pEDTM_scalerRate[-1]))

                    self.tLT[m] = EDTMcut_array.sum()/self.var[m][1]
                    
                    if self.tLT[m] > 1.0:
                        self.tLT[m] = 1.0

                    # tLT_corr[m] = self.tLT[m]*tLT_corr_factor
                      
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
    
                self.tLT = EDTMcut_array.sum()/self.var[1]
    
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
                    
                    pTRIG1_scl_noEDTM = self.var[m][5] - self.var[m][4]
                    pTRIG2_scl_noEDTM = self.var[m][6] - self.var[m][4]
                    pTRIG6_scl_noEDTM = self.var[m][7] - self.var[m][4]
                    
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
                        
                    if T1cpuLT > 1.0:
                        T1cpuLT = 1.0
                    if T2cpuLT > 1.0:
                        T2cpuLT = 1.0
                    if T6cpuLT > 1.0:
                        T6cpuLT = 1.0    
                        
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
                    pEDTM_scaler = data_obj.Branches[m]['P.EDTM.scaler']
                    scalerTime = data_obj.Branches[m]['P.1MHz.scalerTime']
                    pTRIG6_scaler = data_obj.Branches[m]['P.pTRIG6.scaler']
                    
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
                    pEDTM_scaler = data_obj.Branches[m]['P.EDTM.scaler']
                    pTRIG1_scaler = data_obj.Branches[m]['P.pTRIG1.scaler']
                    pTRIG2_scaler = data_obj.Branches[m]['P.pTRIG2.scaler']
                    pTRIG6_scaler = data_obj.Branches[m]['P.pTRIG6.scaler']
                    
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
                                      many=self.many,cut_lim='tgt_boil')

# C12_shms_trk_eff = get_spec_eff(C12_data,eff_type='SHMS',curr='BCM1',many=True)
# C12_trk_eff = C12_shms_trk_eff.calc_eff()

LD2_shms_trk_eff = get_spec_eff(LD2_data,eff_type='SHMS',curr='BCM1',many=True)
LD2_trk_eff = LD2_shms_trk_eff.calc_eff()

# C12_cpuLT = get_spec_eff(C12_data,eff_type='cpuLT',curr='BCM1',many=True)
# C12_allLT = C12_cpuLT.calc_eff()
# C12_LT = {}
# for r in C12_allLT:
#     C12_LT[r] = C12_allLT[r][1][0]

LD2_cpuLT = get_spec_eff(LD2_data,eff_type='cpuLT',curr='BCM1',many=True)
LD2_allLT = LD2_cpuLT.calc_eff()
LD2_LT = {}
for r in LD2_allLT:
    LD2_LT[r] = LD2_allLT[r][1][0]

#%% LD2 get PS for each run
    
LD2_PS = {}

for r in LD2_runs:
    PS=\
        db.retrieve('deuteron_db.db','PS2','RUN_LIST_UPDATED',
                    where = f"run=\'{r}\'")[0][0]

    LD2_PS[r] = PS

#%% C12 initialize current cuts for each set
C12_BCM1_cut = C.current_cut(C12_data,current='BCM1',many=True,cut_lim='tgt_boil')
# C12_BCM2_cut = C.current_cut(C12_data,current='BCM2',many=True,cut_lim='tgt_boil')
# C12_BCM4A_cut = C.current_cut(C12_data,current='BCM4A',many=True,cut_lim='tgt_boil')
# C12_BCM4B_cut = C.current_cut(C12_data,current='BCM4B',many=True,cut_lim='tgt_boil')
# C12_BCM4C_cut = C.current_cut(C12_data,current='BCM4C',many=True,cut_lim='tgt_boil')

C12_BCM1_cut_arrays = C12_BCM1_cut()
# C12_BCM2_cut_arrays = C12_BCM2_cut()
# C12_BCM4A_cut_arrays = C12_BCM4A_cut()
# C12_BCM4B_cut_arrays = C12_BCM4B_cut()
# C12_BCM4C_cut_arrays = C12_BCM4C_cut()

C12_1char_cut = {}
C12_2char_cut = {}
C12_Achar_cut = {}
C12_Bchar_cut = {}
C12_Cchar_cut = {}

C12_pTRIG2scl_cut = {}

C12_1S1XS1Yscl_cut = {}
C12_2S1XS1Yscl_cut = {}
C12_AS1XS1Yscl_cut = {}
C12_BS1XS1Yscl_cut = {}
C12_CS1XS1Yscl_cut = {}

C12_EDTMscl_cut = {}
C12_sclTime_cut = {}
C12_avg_current = []

for r in C12_runs:
    # C12_Achar_cut[r] = C12_BCM4A_cut.scaler(C12_data.Branches[r]['P.BCM4A.scalerCharge'],C12_BCM4A_cut.cut_array[r])
    # C12_pTRIG2scl_cut[r] = C12_BCM4A_cut.scaler(C12_data.Branches[r]['P.pTRIG2.scaler'],C12_BCM4A_cut.cut_array[r])
    # C12_AS1XS1Yscl_cut[r] = C12_BCM4A_cut.scaler(C12_data.Branches[r]['P.S1XS1Y.scaler'],C12_BCM4A_cut.cut_array[r])
    # C12_EDTMscl_cut[r] = C12_BCM4A_cut.scaler(C12_data.Branches[r]['P.EDTM.scaler'],C12_BCM4A_cut.cut_array[r])
    # C12_sclTime_cut[r] = C12_BCM4A_cut.scaler(C12_data.Branches[r]['P.1MHz.scalerTime'],C12_BCM4A_cut.cut_array[r])   
    # C12_avg_current.append(C12_Achar_cut[r]/C12_sclTime_cut[r])
    
    C12_1char_cut[r] = C12_BCM1_cut.scaler(C12_data.Branches[r]['P.BCM1.scalerCharge'],C12_BCM1_cut.cut_array[r])
    C12_pTRIG2scl_cut[r] = C12_BCM1_cut.scaler(C12_data.Branches[r]['P.pTRIG2.scaler'],C12_BCM1_cut.cut_array[r])
    C12_1S1XS1Yscl_cut[r] = C12_BCM1_cut.scaler(C12_data.Branches[r]['P.S1XS1Y.scaler'],C12_BCM1_cut.cut_array[r])
    C12_EDTMscl_cut[r] = C12_BCM1_cut.scaler(C12_data.Branches[r]['P.EDTM.scaler'],C12_BCM1_cut.cut_array[r])
    C12_sclTime_cut[r] = C12_BCM1_cut.scaler(C12_data.Branches[r]['P.1MHz.scalerTime'],C12_BCM1_cut.cut_array[r])    
    C12_avg_current.append(C12_1char_cut[r]/C12_sclTime_cut[r])
    
    # C12_2char_cut[r] = C12_BCM2_cut.scaler(C12_data.Branches[r]['P.BCM2.scalerCharge'],C12_BCM2_cut.cut_array[r])
    # C12_2S1XS1Yscl_cut[r] = C12_BCM2_cut.scaler(C12_data.Branches[r]['P.S1XS1Y.scaler'],C12_BCM2_cut.cut_array[r])

    # C12_Bchar_cut[r] = C12_BCM4B_cut.scaler(C12_data.Branches[r]['P.BCM4B.scalerCharge'],C12_BCM4B_cut.cut_array[r])
    # C12_BS1XS1Yscl_cut[r] = C12_BCM4B_cut.scaler(C12_data.Branches[r]['P.S1XS1Y.scaler'],C12_BCM4B_cut.cut_array[r])

    # C12_Cchar_cut[r] = C12_BCM4C_cut.scaler(C12_data.Branches[r]['P.BCM4C.scalerCharge'],C12_BCM4C_cut.cut_array[r])
    # C12_CS1XS1Yscl_cut[r] = C12_BCM4C_cut.scaler(C12_data.Branches[r]['P.S1XS1Y.scaler'],C12_BCM4C_cut.cut_array[r])

#%% LD2 initialize current cuts for each set
LD2_BCM1_cut = C.current_cut(LD2_data,current='BCM1',many=True,cut_lim='tgt_boil')
# LD2_BCM2_cut = C.current_cut(LD2_data,current='BCM2',many=True,cut_lim='tgt_boil')
# LD2_BCM4A_cut = C.current_cut(LD2_data,current='BCM4A',many=True,cut_lim='tgt_boil')
# LD2_BCM4B_cut = C.current_cut(LD2_data,current='BCM4B',many=True,cut_lim='tgt_boil')
# LD2_BCM4C_cut = C.current_cut(LD2_data,current='BCM4C',many=True,cut_lim='tgt_boil')

LD2_BCM1_cut_arrays = LD2_BCM1_cut()
# LD2_BCM2_cut_arrays = LD2_BCM2_cut()
# LD2_BCM4A_cut_arrays = LD2_BCM4A_cut()
# LD2_BCM4B_cut_arrays = LD2_BCM4B_cut()
# LD2_BCM4C_cut_arrays = LD2_BCM4C_cut()

LD2_1char_cut = {}
LD2_2char_cut = {}
LD2_Achar_cut = {}
LD2_Bchar_cut = {}
LD2_Cchar_cut = {}

LD2_pTRIG2scl_cut = {}

LD2_1S1XS1Yscl_cut = {}
LD2_2S1XS1Yscl_cut = {}
LD2_AS1XS1Yscl_cut = {}
LD2_BS1XS1Yscl_cut = {}
LD2_CS1XS1Yscl_cut = {}

LD2_EDTMscl_cut = {}
LD2_sclTime_cut = {}
LD2_avg_current = []

for r in LD2_runs:
    # LD2_Achar_cut[r] = LD2_BCM4A_cut.scaler(LD2_data.Branches[r]['P.BCM4A.scalerCharge'],LD2_BCM4A_cut.cut_array[r])
    # LD2_pTRIG2scl_cut[r] = LD2_BCM4A_cut.scaler(LD2_data.Branches[r]['P.pTRIG2.scaler'],LD2_BCM4A_cut.cut_array[r])
    # LD2_AS1XS1Yscl_cut[r] = LD2_BCM4A_cut.scaler(LD2_data.Branches[r]['P.S1XS1Y.scaler'],LD2_BCM4A_cut.cut_array[r])
    # LD2_EDTMscl_cut[r] = LD2_BCM4A_cut.scaler(LD2_data.Branches[r]['P.EDTM.scaler'],LD2_BCM4A_cut.cut_array[r])
    # LD2_sclTime_cut[r] = LD2_BCM4A_cut.scaler(LD2_data.Branches[r]['P.1MHz.scalerTime'],LD2_BCM4A_cut.cut_array[r])   
    # LD2_avg_current.append(LD2_Achar_cut[r]/LD2_sclTime_cut[r])
    
    LD2_1char_cut[r] = LD2_BCM1_cut.scaler(LD2_data.Branches[r]['P.BCM1.scalerCharge'],LD2_BCM1_cut.cut_array[r])
    LD2_pTRIG2scl_cut[r] = LD2_BCM1_cut.scaler(LD2_data.Branches[r]['P.pTRIG2.scaler'],LD2_BCM1_cut.cut_array[r])
    LD2_1S1XS1Yscl_cut[r] = LD2_BCM1_cut.scaler(LD2_data.Branches[r]['P.S1XS1Y.scaler'],LD2_BCM1_cut.cut_array[r])
    LD2_EDTMscl_cut[r] = LD2_BCM1_cut.scaler(LD2_data.Branches[r]['P.EDTM.scaler'],LD2_BCM1_cut.cut_array[r])
    LD2_sclTime_cut[r] = LD2_BCM1_cut.scaler(LD2_data.Branches[r]['P.1MHz.scalerTime'],LD2_BCM1_cut.cut_array[r])    
    LD2_avg_current.append(LD2_1char_cut[r]/LD2_sclTime_cut[r])
    
    # LD2_2char_cut[r] = LD2_BCM2_cut.scaler(LD2_data.Branches[r]['P.BCM2.scalerCharge'],LD2_BCM2_cut.cut_array[r])
    # LD2_2S1XS1Yscl_cut[r] = LD2_BCM2_cut.scaler(LD2_data.Branches[r]['P.S1XS1Y.scaler'],LD2_BCM2_cut.cut_array[r])

    # LD2_Bchar_cut[r] = LD2_BCM4B_cut.scaler(LD2_data.Branches[r]['P.BCM4B.scalerCharge'],LD2_BCM4B_cut.cut_array[r])
    # LD2_BS1XS1Yscl_cut[r] = LD2_BCM4B_cut.scaler(LD2_data.Branches[r]['P.S1XS1Y.scaler'],LD2_BCM4B_cut.cut_array[r])

    # LD2_Cchar_cut[r] = LD2_BCM4C_cut.scaler(LD2_data.Branches[r]['P.BCM4C.scalerCharge'],LD2_BCM4C_cut.cut_array[r])
    # LD2_CS1XS1Yscl_cut[r] = LD2_BCM4C_cut.scaler(LD2_data.Branches[r]['P.S1XS1Y.scaler'],LD2_BCM4C_cut.cut_array[r])

#%% C12 define cuts for trk and no trk yields
pNGCer_cut = C.WCUT(3.0,np.inf,'pNGCer_cut') 
pdelta_cut = C.WCUT(-10,22,'pdelta_cut')
# pztar_cut = C.WCUT(-1,1,'pztar_cut')

C12_all_cuts_notrk = {}
C12_all_cuts_trk = {}
for r in C12_runs:
    #make ztar cut around mean
    pztar = C12_data.Branches[r]['P.react.z']
    
    hc = B.histo(pztar,range=(-5.,5.),bins=100)
    hc.fit(plot_fit=False)
    
    cmin = hc.mean.value - 1.5*hc.sigma.value
    cmax = hc.mean.value + 1.5*hc.sigma.value

    pztar_cut = C.WCUT(cmin,cmax,'pztar_cut')
    
    
    ngcerNpeSum = C12_data.Branches[r]['P.ngcer.npeSum']
    pcalEtotNorm = C12_data.Branches[r]['P.cal.etotnorm']
    pEDTM = C12_data.Branches[r]['T.coin.pEDTM_tdcTimeRaw']
    pcalEtotTrkNorm = C12_data.Branches[r]['P.cal.etottracknorm']
    pdp = C12_data.Branches[r]['P.gtr.dp']

    
    c_noedtm = C.noEDTM(pEDTM)
    c_ngcer = pNGCer_cut(ngcerNpeSum)
    c_caletotnorm = C.shms_calPID(pcalEtotNorm)
    c_caletottracknorm = C.shms_calPID(pcalEtotTrkNorm)
    c_pdelta = pdelta_cut(pdp)
    c_ztar = pztar_cut(pztar)
    
       
    cuts_all_notrk = c_noedtm & c_ngcer & c_caletotnorm & c_ztar
    cuts_all_trk = c_noedtm & c_pdelta & c_caletottracknorm & c_ztar
    
    C12_all_cuts_notrk[r] = cuts_all_notrk
    C12_all_cuts_trk[r] = cuts_all_trk
    
    # h = B.histo(pztar,range=(-5,5),bins=100)
    # hc = B.histo(pztar[c_ztar],range=(-5,5),bins=100)
    
    # B.pl.figure()
    # h.plot(filled=False)
    # hc.plot()
    
# =============================================================================
#     if r == 21063:
#         h10 = B.histo(ngcerNpeSum,range=(0.5,20),bins=100)
#         h20 = B.histo(pcalEtotNorm,range=(0.1,2),bins=100)
#         h30 = B.histo(pEDTM,range=(0,1000),bins=100)
#         h40 = B.histo(pcalEtotTrkNorm,range=(0.1,2),bins=100)
#         h50 = B.histo(pdp,range=(-15,25),bins=100)
#         
#         h1 = B.histo(ngcerNpeSum[c_ngcer],range=(0.5,20),bins=100,
#                      title=f'NGCER NPE SUM: RUN {r}',xlabel='',ylabel='')
#         h2 = B.histo(pcalEtotNorm[c_caletotnorm],range=(0.1,2),bins=100,
#                      title='Calorimeter $E_{dep}/p$:'+f' RUN {r}',xlabel='',ylabel='')
#         h3 = B.histo(pEDTM[c_noedtm],range=(0,1000),bins=100,
#                      title='')
#         h4 = B.histo(pcalEtotTrkNorm[c_caletottracknorm],range=(0.1,2),bins=100,
#                      title='Calorimeter Tracking $E_{dep}/p$:'+f' RUN {r}',
#                      xlabel='',ylabel='')
#         h5 = B.histo(pdp[c_pdelta],range=(-15,25),bins=100, 
#                      title=f'SHMS Delta: RUN {r}',xlabel='',ylabel='')
#         # h6 = B.histo(pztar[c_ztar],range=(-5,5),bins=100)
#         hlist_nocut = [h10,h20,h30,h40,h50]
#         hlist = [h1, h2, h3, h4, h5]
#         # for i in range(len(hlist)):
#         #     B.pl.figure()
#         #     hlist_nocut[i].plot(filled=False,color=cmp['Paired'].colors[1])
#         #     hlist[i].plot(color=cmp['Paired'].colors[0])
#             
#         h20 = B.histo(pcalEtotNorm,range=(0.1,2),bins=100)
#         h50 = B.histo(pdp,range=(-15,25),bins=100)
#         
#         
#         h2 = B.histo(pcalEtotNorm[cuts_all_notrk],range=(0.1,2),bins=100,
#                      title='Calorimeter $E_{dep}/p$:'+f' RUN {r}',xlabel='',ylabel='')
#         h5 = B.histo(pdp[cuts_all_trk],range=(-15,25),bins=100, 
#                      title=f'SHMS Delta: RUN {r}',xlabel='',ylabel='')
# 
#         B.pl.figure()
#         h20.plot(filled=False,color=cmp['Paired'].colors[1],label='No Cuts')
#         h2.plot(color=cmp['Paired'].colors[0],label='Cuts')
#         B.pl.vlines([0.7,1.3],ymin=0,ymax=80000,linestyles='--',colors='r')
#         B.pl.text(1.35,70000,'EDTM > 0\nNGCER NPE > 3.0\n0.7 > Cal E/P > 1.3')
#         B.pl.legend()
#         
#         B.pl.figure()
#         h50.plot(filled=False,color=cmp['Paired'].colors[1],label='No Cuts')
#         h5.plot(color=cmp['Paired'].colors[0],label='Cuts')
#         B.pl.vlines([-10,22],ymin=0,ymax=28000,linestyles='--',colors='r')
#         B.pl.text(-5,20000,'EDTM > 0\n0.7 > Cal Tracking E/P > 1.3\n-10 > SHMS Delta > 22')
#         B.pl.legend()
# =============================================================================

#%% C12 yields: scaler, no trk, trk

C12_scl_yield_BCM4A = []
C12_scl_yield_BCM4B = []
C12_scl_yield_BCM4C = []
C12_scl_yield_BCM1 = []
C12_scl_yield_BCM2 = []
C12_scl_yield_BCM4A_err = []
C12_scl_yield_BCM4B_err = []
C12_scl_yield_BCM4C_err = []
C12_scl_yield_BCM1_err = []
C12_scl_yield_BCM2_err = []

C12_notrk_yield_BCM4A = []
C12_notrk_yield_BCM4B = []
C12_notrk_yield_BCM4C = []
C12_notrk_yield_BCM4A_err = []
C12_notrk_yield_BCM4B_err = []
C12_notrk_yield_BCM4C_err = []

C12_trk_yield_BCM4A = []
C12_trk_yield_BCM4B = []
C12_trk_yield_BCM4C = []
C12_trk_yield_BCM4A_err = []
C12_trk_yield_BCM4B_err = []
C12_trk_yield_BCM4C_err = []

for r in C12_runs:
    # declare and combine cuts for easy use
    # BCM4A_curr_cut = C12_BCM4A_cut.cut_array[r]
    # BCM4B_curr_cut = C12_BCM4B_cut.cut_array[r]
    # BCM4C_curr_cut = C12_BCM4C_cut.cut_array[r]
    
    cuts_all_1_notrk = C12_BCM1_cut_arrays[r] & C12_all_cuts_notrk[r]
    # cuts_all_4A_notrk = C12_BCM4A_cut_arrays[r] & C12_all_cuts_notrk[r]
    # cuts_all_4B_notrk = C12_BCM4B_cut_arrays[r] & C12_all_cuts_notrk[r] 
    # cuts_all_4C_notrk = C12_BCM4C_cut_arrays[r] & C12_all_cuts_notrk[r]
    
    cuts_all_1_trk = C12_BCM1_cut_arrays[r] & C12_all_cuts_trk[r]
    # cuts_all_4A_trk = C12_BCM4A_cut_arrays[r] & C12_all_cuts_trk[r]
    # cuts_all_4B_trk = C12_BCM4B_cut_arrays[r] & C12_all_cuts_trk[r] 
    # cuts_all_4C_trk = C12_BCM4C_cut_arrays[r] & C12_all_cuts_trk[r]
    
    # get variables used to calculate yields
    # y_scl = C12_data.Branches[r]['P.S1XS1Y.scaler']
    # y_scl = C12_data.Branches[r]['P.pTRIG2.scaler']
    # y_scl = newpTRIG2_scl[r]
    # y_scl = C12_pTRIG2scl_cut[r]
    y_scl = C12_1S1XS1Yscl_cut[r]
    # y_scl = C12_AS1XS1Yscl_cut[r]
    y_notrk = C12_data.Branches[r]['P.cal.etotnorm']
    y_trk = C12_data.Branches[r]['P.gtr.dp'] 
    # EDTM_scl = C12_data.Branches[r]['P.EDTM.scaler']
    EDTM_scl = C12_EDTMscl_cut[r]
    charge = C12_1char_cut[r]/1e3 # muC to mC
    # charge = C12_Achar_cut[r]/1e3 # muC to mC

 
    # subtract EDTM events & normalize by charge    
    BCM4A_scl_yield = (y_scl-EDTM_scl)/charge
    BCM4A_scl_yield_err = np.sqrt(y_scl)/charge
    print(f'Run {r} \nBCM4A scaler yield = ', BCM4A_scl_yield)
    
    # BCM4B_scl_yield = (C12_BS1XS1Yscl_cut[r]-EDTM_scl)/(C12_Bchar_cut[r]/1e3)
    # BCM4B_scl_yield_err = np.sqrt(C12_BS1XS1Yscl_cut[r])/(C12_Bchar_cut[r]/1e3)
    # print('BCM4B scaler yield = ',BCM4B_scl_yield)
    
    # BCM4C_scl_yield = (C12_CS1XS1Yscl_cut[r]-EDTM_scl)/(C12_Cchar_cut[r]/1e3)
    # BCM4C_scl_yield_err = np.sqrt(C12_CS1XS1Yscl_cut[r])/(C12_Cchar_cut[r]/1e3)
    # print('BCM4C scaler yield = ',BCM4C_scl_yield)
    
    # BCM1_scl_yield = (C12_1S1XS1Yscl_cut[r]-EDTM_scl)/(C12_1char_cut[r]/1e3)
    # BCM1_scl_yield_err = np.sqrt(C12_1S1XS1Yscl_cut[r])/(C12_1char_cut[r]/1e3)
    # print('BCM1 scaler yield = ',BCM1_scl_yield)
    
    # BCM2_scl_yield = (C12_2S1XS1Yscl_cut[r]-EDTM_scl)/(C12_2char_cut[r]/1e3)
    # BCM2_scl_yield_err = np.sqrt(C12_2S1XS1Yscl_cut[r])/(C12_2char_cut[r]/1e3)
    # print('BCM2 scaler yield = ',BCM2_scl_yield)
    
    norm = 1/(charge*C12_LT[r]) 
    BCM4A_notrk_h = B.histo(y_notrk[cuts_all_1_notrk],range=(0,2),bins=100,
                                 title=f'Run {r} No Track Yield')
    BCM4A_notrk_yield,BCM4A_notrk_yield_err = BCM4A_notrk_h.sum() 
    BCM4A_notrk_yield = norm*BCM4A_notrk_yield
    BCM4A_notrk_yield_err = norm*BCM4A_notrk_yield_err
    print('notrk yield = ', BCM4A_notrk_yield)
    
    # B.pl.figure()
    # BCM4A_notrk_h.plot()
    
    norm = 1/(charge*C12_LT[r]*C12_trk_eff[r]) 
    BCM4A_trk_h = B.histo(y_trk[cuts_all_1_trk],range=(-15,15),bins=100,
                               title=f'Run {r} Yield w/ Tracking')
    # BCM4A_trk_h = B.histo(y_trk[C12_BCM4A_cut_arrays[r]],range=(-10,22),bins=100,
    #                            title=f'Run {r} Yield w/ Tracking')
    BCM4A_trk_yield,BCM4A_trk_yield_err = BCM4A_trk_h.sum()
    # normalize yield
    print('trk yield no norm = ', BCM4A_trk_yield)
    BCM4A_trk_yield = norm*BCM4A_trk_yield
    BCM4A_trk_yield_err = norm*BCM4A_trk_yield_err
    print('trk norm = ', norm)
    
    print('trk yield = ', BCM4A_trk_yield,'\n')
    
    # B.pl.figure()
    # BCM4A_trk_h.plot()
    
    C12_scl_yield_BCM4A.append(BCM4A_scl_yield)
    # C12_scl_yield_BCM4B.append(BCM4B_scl_yield)
    # C12_scl_yield_BCM4C.append(BCM4C_scl_yield) 
    # C12_scl_yield_BCM1.append(BCM1_scl_yield)
    # C12_scl_yield_BCM2.append(BCM2_scl_yield) 
    C12_scl_yield_BCM4A_err.append(BCM4A_scl_yield_err)
    # C12_scl_yield_BCM4B_err.append(BCM4B_scl_yield_err)
    # C12_scl_yield_BCM4C_err.append(BCM4C_scl_yield_err)
    # C12_scl_yield_BCM1_err.append(BCM1_scl_yield_err)
    # C12_scl_yield_BCM2_err.append(BCM2_scl_yield_err) 
    
    C12_notrk_yield_BCM4A.append(BCM4A_notrk_yield)
    C12_notrk_yield_BCM4A_err.append(BCM4A_notrk_yield_err) 
    
    C12_trk_yield_BCM4A.append(BCM4A_trk_yield) 
    C12_trk_yield_BCM4A_err.append(BCM4A_trk_yield_err)

# B.pl.figure(figsize=(15,15))
# B.pl.suptitle('C12 Normalized Scaler Yield')
# B.pl.subplot(2,3,1)
# B.plot_exp(C12_avg_current,C12_scl_yield_BCM4A,C12_scl_yield_BCM4A_err)
# B.pl.title('BCM4A')
# B.pl.ylabel('Normalized Scaler Yield: Y/Q')
# B.pl.xlabel('Avg. Current [$\mu$A]')
# # B.pl.figure()
# B.pl.subplot(2,3,2)
# B.plot_exp(C12_avg_current,C12_scl_yield_BCM4B,C12_scl_yield_BCM4B_err)
# B.pl.title('BCM4B')
# # B.pl.figure()
# B.pl.subplot(2,3,3)
# B.plot_exp(C12_avg_current,C12_scl_yield_BCM4C,C12_scl_yield_BCM4C_err)
# B.pl.title('BCM4C')
# B.pl.subplot(2,3,4)
# B.plot_exp(C12_avg_current,C12_scl_yield_BCM1,C12_scl_yield_BCM1_err)
# B.pl.title('BCM1')
# B.pl.subplot(2,3,5)
# B.plot_exp(C12_avg_current,C12_scl_yield_BCM2,C12_scl_yield_BCM2_err)
# B.pl.title('BCM2')

plot_currname = 'BCM1'
plot_current = np.array(C12_avg_current)
plot_yscl = np.array(C12_scl_yield_BCM4A)
plot_ysclerr = np.array(C12_scl_yield_BCM4A_err)
plot_yntrk = C12_notrk_yield_BCM4A
plot_yntrkerr = C12_notrk_yield_BCM4A_err
plot_ytrk = C12_trk_yield_BCM4A
plot_ytrkerr = C12_trk_yield_BCM4A_err

B.pl.figure()
B.pl.suptitle('C12 Normalized Scaler Yield')
B.plot_exp(plot_current,plot_yscl,plot_ysclerr,
           color='black')
for i, label in enumerate(C12_runs):
    B.pl.annotate(label, (plot_current[i], plot_yscl[i]),
                  textcoords="offset points", xytext=(0,10), ha='left')
    
linefit = B.linefit(plot_current,plot_yscl,plot_ysclerr,plot_fit=False)
linefit.plot(label = f'Slope = {linefit.slope:0.1f}$\pm${linefit.sigma_s:0.1f}\n'+\
                     f'Offset = {linefit.offset:4.1f}$\pm${linefit.sigma_o:4.1f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Scaler Yield: Y/Q')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.ylim(2.7e6,2.9e6)
# fit_legend_line(linefit)
B.pl.legend()

B.pl.figure()
B.pl.suptitle('C12 Fractional Normalized Scaler Yield')
B.plot_exp(plot_current,plot_yscl/plot_yscl[0],
           plot_ysclerr/plot_yscl[0],label='C12_data')
linefit = B.linefit(plot_current,
                    plot_yscl/plot_yscl[0],
                    plot_ysclerr/plot_yscl[0],plot_fit=False)
linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Scaler Yield: Y/Q')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.ylim(0.9,1.1)
# fit_legend_line(linefit)
B.pl.legend()


B.pl.figure()
B.plot_exp(plot_current,plot_yntrk,plot_yntrkerr)
B.pl.suptitle('C12 Normalized No Tracking Yield')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized No Tracking Yield: Y/(Q*tLT)')
B.pl.xlabel('Avg. Current [$\mu$A]')
for i, label in enumerate(C12_runs):
    B.pl.annotate(label, (plot_current[i], plot_yntrk[i]),
                  textcoords="offset points", xytext=(0,10), ha='left')
# linefit = B.linefit(np.array(C12_Acurrent),
#                     np.array(C12_notrk_yield_BCM4A),
#                     np.array(C12_notrk_yield_BCM4A_err),plot_fit=False)
# linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
#                      f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
#                      '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
# B.pl.legend()
    
B.pl.figure()
B.pl.suptitle('C12 Fractional Normalized No Tracking Yield')
B.plot_exp(plot_current,plot_yntrk/plot_yntrk[0],
           plot_yntrkerr/plot_yntrk[0],label='C12_data_notrk')
B.pl.ylim(0.9,1.1)
linefit = B.linefit(plot_current,
                    plot_yntrk/plot_yntrk[0],
                    plot_yntrkerr/plot_yntrk[0],plot_fit=False)
linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized No Tracking Yield: Y/(Q*tLT)')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.legend()

B.pl.figure()
B.plot_exp(plot_current,plot_ytrk,plot_ytrkerr)
B.pl.suptitle('C12 Normalized Tracking Yield')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Tracking Yield: Y/(Q*tLT*TrkEff)')
B.pl.xlabel('Avg. Current [$\mu$A]')
for i, label in enumerate(C12_runs):
    B.pl.annotate(label, (plot_current[i], plot_ytrk[i]),
                  textcoords="offset points", xytext=(0,10), ha='left')
# linefit = B.linefit(np.array(C12_Acurrent),
#                     np.array(C12_trk_yield_BCM4A),
#                     np.array(C12_trk_yield_BCM4A_err),plot_fit=False)
# linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
#                      f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
#                      '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
# B.pl.legend()
    
B.pl.figure()
B.pl.suptitle('C12 Fractional Normalized Tracking Yield')
B.plot_exp(plot_current,plot_ytrk/plot_ytrk[0],
           plot_ytrkerr/plot_ytrk[0],label='C12_data_trk')
B.pl.ylim(0.9,1.1)
linefit = B.linefit(plot_current,
                    plot_ytrk/plot_ytrk[0],
                    plot_ytrkerr/plot_ytrk[0],plot_fit=False)
linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Tracking Yield: Y/(Q*tLT*TrkEff)')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.legend()

#%% C12 NORMALIZED YIELD PLOT
B.pl.figure(figsize=(8,10),layout='constrained')
# B.pl.suptitle('C12 Normalized Yield')

B.pl.subplot(3,1,1)
B.plot_exp(plot_current,plot_yscl,plot_ysclerr,color=cmp['Paired'].colors[5],
          markersize=10,capsize=8)
for i, label in enumerate(C12_runs):
    B.pl.annotate(label, (plot_current[i], plot_yscl[i]),
                  textcoords="offset points", xytext=(0,-15), ha='left')
B.pl.ylim(2.7e6,2.9e6)
B.pl.xlim(5,75)
B.pl.title(plot_currname)
B.pl.text(8,2.72e6,'C12 Scaler Yield',fontsize=24)

B.pl.subplot(3,1,2)
B.plot_exp(plot_current,plot_yntrk,plot_yntrkerr,color=cmp['Paired'].colors[5],
           markersize=10,capsize=8)
B.pl.ylim(11000,14000)
B.pl.xlim(5,75)
B.pl.ylabel('Normalized Yield')
B.pl.text(8,11500,'C12 No Tracking Yield',fontsize=24)

B.pl.subplot(3,1,3)
B.plot_exp(plot_current,plot_ytrk,plot_ytrkerr,color=cmp['Paired'].colors[5],
           markersize=10,capsize=8)
B.pl.ylim(10000,13000)
B.pl.xlim(5,75)
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.text(8,10500,'C12 Tracking Yield',fontsize=24)

B.pl.savefig('target boiling/plots/C12_norm_yield_all.png')


#%% LD2 define cuts for trk and no trk yields
pNGCer_cut = C.WCUT(3.0,np.inf,'pNGCer_cut') 
pdelta_cut = C.WCUT(-10,22,'pdelta_cut')

LD2_all_cuts_notrk = {}
LD2_all_cuts_trk = {}
for r in LD2_runs:
    ngcerNpeSum = LD2_data.Branches[r]['P.ngcer.npeSum']
    pcalEtotNorm = LD2_data.Branches[r]['P.cal.etotnorm']
    pEDTM = LD2_data.Branches[r]['T.coin.pEDTM_tdcTimeRaw']
    pcalEtotTrkNorm = LD2_data.Branches[r]['P.cal.etottracknorm']
    pdp = LD2_data.Branches[r]['P.gtr.dp']
    
    c_noedtm = C.noEDTM(pEDTM)
    c_ngcer = pNGCer_cut(ngcerNpeSum)
    c_caletotnorm = C.shms_calPID(pcalEtotNorm)
    c_caletottracknorm = C.shms_calPID(pcalEtotTrkNorm)
    c_pdelta = pdelta_cut(pdp)
    
       
    cuts_all_notrk = c_noedtm & c_ngcer & c_caletotnorm
    cuts_all_trk = c_noedtm & c_pdelta & c_caletottracknorm
    
    LD2_all_cuts_notrk[r] = cuts_all_notrk
    LD2_all_cuts_trk[r] = cuts_all_trk
#%% LD2 yields: scaler, no trk, trk

LD2_scl_yield_BCM4A = []
LD2_scl_yield_BCM4B = []
LD2_scl_yield_BCM4C = []
LD2_scl_yield_BCM4A_err = []
LD2_scl_yield_BCM4B_err = []
LD2_scl_yield_BCM4C_err = []

LD2_notrk_yield_BCM4A = []
LD2_notrk_yield_BCM4B = []
LD2_notrk_yield_BCM4C = []
LD2_notrk_yield_BCM4A_err = []
LD2_notrk_yield_BCM4B_err = []
LD2_notrk_yield_BCM4C_err = []

LD2_trk_yield_BCM4A = []
LD2_trk_yield_BCM4B = []
LD2_trk_yield_BCM4C = []
LD2_trk_yield_BCM4A_err = []
LD2_trk_yield_BCM4B_err = []
LD2_trk_yield_BCM4C_err = []

for r in LD2_runs:
    # declare and combine cuts for easy use
    # BCM4A_curr_cut = LD2_BCM4A_cut.cut_array[r]
    # BCM4B_curr_cut = LD2_BCM4B_cut.cut_array[r]
    # BCM4C_curr_cut = LD2_BCM4C_cut.cut_array[r]
    
    cuts_all_1_notrk = LD2_BCM1_cut_arrays[r] & LD2_all_cuts_notrk[r]
    # cuts_all_4A_notrk = LD2_BCM4A_cut_arrays[r] & LD2_all_cuts_notrk[r]
    # cuts_all_4B_notrk = LD2_BCM4B_cut_arrays[r] & LD2_all_cuts_notrk[r] 
    # cuts_all_4C_notrk = LD2_BCM4C_cut_arrays[r] & LD2_all_cuts_notrk[r]
    
    cuts_all_1_trk = LD2_BCM1_cut_arrays[r] & LD2_all_cuts_trk[r]
    # cuts_all_4A_trk = LD2_BCM4A_cut_arrays[r] & LD2_all_cuts_trk[r]
    # cuts_all_4B_trk = LD2_BCM4B_cut_arrays[r] & LD2_all_cuts_trk[r] 
    # cuts_all_4C_trk = LD2_BCM4C_cut_arrays[r] & LD2_all_cuts_trk[r]
    
    # get variables used to calculate yields
    # y_scl = LD2_data.Branches[r]['P.S1XS1Y.scaler']
    # y_scl = LD2_data.Branches[r]['P.pTRIG2.scaler']
    # y_scl = newpTRIG2_scl[r]
    y_scl = LD2_1S1XS1Yscl_cut[r]
    y_notrk = LD2_data.Branches[r]['P.cal.etotnorm']
    y_trk = LD2_data.Branches[r]['P.gtr.dp'] 
    # EDTM_scl = LD2_data.Branches[r]['P.EDTM.scaler']
    EDTM_scl = LD2_EDTMscl_cut[r]
    charge = LD2_1char_cut[r]/1e3
 
    # subtract EDTM events      
    # BCM4A_scl_yield = y_scl[BCM4A_curr_cut][-1] - EDTM_scl[BCM4A_curr_cut][-1]
    # normalize by charge
    # BCM4A_scl_yield = BCM4A_scl_yield/LD2_Acharge[r]
    BCM4A_scl_yield = (y_scl-EDTM_scl)/charge
    # BCM4A_scl_yield = y_scl[-1]/newLD2_Acharge[r]
    # BCM4A_scl_yield_err = np.sqrt(y_scl[BCM4A_curr_cut][-1])/LD2_Acharge[r]
    BCM4A_scl_yield_err = np.sqrt(y_scl)/charge
    print(f'Run {r} \nscaler yield = ', BCM4A_scl_yield)
    
    # BCM4B_scl_yield = y_scl[BCM4B_curr_cut][-1] - EDTM_scl[BCM4B_curr_cut][-1]
    # BCM4B_scl_yield = BCM4B_scl_yield/LD2_Bcharge[r]
    # BCM4B_scl_yield_err = np.sqrt(y_scl[BCM4B_curr_cut][-1])/LD2_Bcharge[r]
    # print(BCM4B_scl_yield)
    
    # BCM4C_scl_yield = y_scl[BCM4C_curr_cut][-1] - EDTM_scl[BCM4C_curr_cut][-1]
    # BCM4C_scl_yield = BCM4C_scl_yield/LD2_Ccharge[r]
    # BCM4C_scl_yield_err = np.sqrt(y_scl[BCM4C_curr_cut][-1])/LD2_Ccharge[r]
    # print(BCM4C_scl_yield)
    
    norm = LD2_PS[r]/(charge*LD2_LT[r]) 
    BCM4A_notrk_h = B.histo(y_notrk[cuts_all_1_notrk],range=(0,2),bins=100,
                                 title=f'Run {r} No Track Yield')
    BCM4A_notrk_yield,BCM4A_notrk_yield_err = BCM4A_notrk_h.sum() 
    BCM4A_notrk_yield = norm*BCM4A_notrk_yield
    BCM4A_notrk_yield_err = norm*BCM4A_notrk_yield_err
    print('notrk yield = ', BCM4A_notrk_yield)
    
    # B.pl.figure()
    # BCM4A_notrk_h.plot()
    
    norm = LD2_PS[r]/(charge*LD2_LT[r]*LD2_trk_eff[r]) 
    BCM4A_trk_h = B.histo(y_trk[cuts_all_1_trk],range=(-15,15),bins=100,
                               title=f'Run {r} Yield w/ Tracking')
    # BCM4A_trk_h = B.histo(y_trk[LD2_BCM4A_cut_arrays[r]],range=(-10,22),bins=100,
    #                            title=f'Run {r} Yield w/ Tracking')
    BCM4A_trk_yield,BCM4A_trk_yield_err = BCM4A_trk_h.sum()
    # normalize yield
    print('trk yield no norm = ', BCM4A_trk_yield)
    BCM4A_trk_yield = norm*BCM4A_trk_yield
    BCM4A_trk_yield_err = norm*BCM4A_trk_yield_err
    print('trk norm = ', norm)
    
    print('trk yield = ', BCM4A_trk_yield,'\n')
    
    # B.pl.figure()
    # BCM4A_trk_h.plot()
    
    LD2_scl_yield_BCM4A.append(BCM4A_scl_yield)
    # LD2_scl_yield_BCM4B.append(BCM4B_scl_yield)
    # LD2_scl_yield_BCM4C.append(BCM4C_scl_yield)   
    LD2_scl_yield_BCM4A_err.append(BCM4A_scl_yield_err)
    # LD2_scl_yield_BCM4B_err.append(BCM4B_scl_yield_err)
    # LD2_scl_yield_BCM4C_err.append(BCM4C_scl_yield_err) 
    
    LD2_notrk_yield_BCM4A.append(BCM4A_notrk_yield)
    LD2_notrk_yield_BCM4A_err.append(BCM4A_notrk_yield_err) 
    
    LD2_trk_yield_BCM4A.append(BCM4A_trk_yield) 
    LD2_trk_yield_BCM4A_err.append(BCM4A_trk_yield_err)

# B.pl.figure()
# B.pl.suptitle('LD2 Normalized Scaler Yield')
# B.pl.subplot(1,3,1)
# B.plot_exp(LD2_Acurrent,LD2_scl_yield_BCM4A,LD2_scl_yield_BCM4A_err)
# B.pl.title('BCM4A')
# B.pl.ylabel('Normalized Scaler Yield: Y/Q')
# B.pl.figure()
# B.pl.subplot(1,3,2)
# B.plot_exp(LD2_Bcurrent,LD2_scl_yield_BCM4B,LD2_scl_yield_BCM4B_err)
# B.pl.title('BCM4B')
# B.pl.xlabel('Avg. Current [$\mu$A]')
# B.pl.figure()
# B.pl.subplot(1,3,3)
# B.plot_exp(LD2_Ccurrent,LD2_scl_yield_BCM4C,LD2_scl_yield_BCM4C_err)
# B.pl.title('BCM4C')
# plot_current = np.array(LD2_Acurrent)/LD2_Acurrent[0]
plot_currname = 'BCM1'
plot_current = np.array(LD2_avg_current)
plot_yscl = LD2_scl_yield_BCM4A
plot_ysclerr = LD2_scl_yield_BCM4A_err
plot_yntrk = LD2_notrk_yield_BCM4A
plot_yntrkerr = LD2_notrk_yield_BCM4A_err
plot_ytrk = LD2_trk_yield_BCM4A
plot_ytrkerr = LD2_trk_yield_BCM4A_err

B.pl.figure()
B.pl.suptitle('LD2 Normalized Scaler Yield')
B.plot_exp(plot_current,plot_yscl,plot_ysclerr,label='LD2_data')
for i, label in enumerate(LD2_runs):
    B.pl.annotate(label, (plot_current[i], plot_yscl[i]),
                  textcoords="offset points", xytext=(0,10), ha='left')
# linefit = B.linefit(plot_current,plot_yscl,plot_ysclerr,plot_fit=False)
# linefit.plot(label = f'Slope = {linefit.slope:0.1f}$\pm${linefit.sigma_s:0.1f}\n'+\
#                      f'Offset = {linefit.offset:4.1f}$\pm${linefit.sigma_o:4.1f}\n'+\
#                      '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Scaler Yield: Y/Q')
B.pl.xlabel('Avg. Current [$\mu$A]')
# fit_legend_line(linefit)
# B.pl.legend()

B.pl.figure()
B.pl.suptitle('LD2 Fractional Normalized Scaler Yield')
B.plot_exp(plot_current,plot_yscl/plot_yscl[0],
           plot_ysclerr/plot_yscl[0],label='LD2_data')
linefit_scl = B.linefit(plot_current,
                    plot_yscl/plot_yscl[0],
                    plot_ysclerr/plot_yscl[0],plot_fit=False)
linefit_scl.plot(label = f'Slope = {linefit_scl.slope:0.2e}$\pm${linefit_scl.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_scl.offset:0.2f}$\pm${linefit_scl.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_scl.chi_red:0.2f}')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Scaler Yield: Y/Q')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.ylim(0.9,1.1)
# fit_legend_line(linefit)
B.pl.legend()


B.pl.figure()
B.plot_exp(plot_current,plot_yntrk,plot_yntrkerr)
B.pl.suptitle('LD2 Normalized No Tracking Yield')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized No Tracking Yield: Y/(Q*tLT)')
B.pl.xlabel('Avg. Current [$\mu$A]')
for i, label in enumerate(LD2_runs):
    B.pl.annotate(label, (plot_current[i], plot_yntrk[i]),
                  textcoords="offset points", xytext=(0,10), ha='left')
# linefit = B.linefit(np.array(LD2_Acurrent),
#                     np.array(LD2_notrk_yield_BCM4A),
#                     np.array(LD2_notrk_yield_BCM4A_err),plot_fit=False)
# linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
#                      f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
#                      '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
# B.pl.legend()
    
B.pl.figure()
B.pl.suptitle('LD2 Fractional Normalized No Tracking Yield')
B.plot_exp(plot_current,plot_yntrk/plot_yntrk[0],
           plot_yntrkerr/plot_yntrk[0],label='LD2_data_notrk')
B.pl.ylim(0.9,1.1)
linefit_notrk = B.linefit(plot_current,
                    plot_yntrk/plot_yntrk[0],
                    plot_yntrkerr/plot_yntrk[0],plot_fit=False)
linefit_notrk.plot(label = f'Slope = {linefit_notrk.slope:0.2e}$\pm${linefit_notrk.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_notrk.offset:0.2f}$\pm${linefit_notrk.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_notrk.chi_red:0.2f}')
B.pl.ylabel('Normalized No Tracking Yield: Y/(Q*tLT)')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.legend()

B.pl.figure()
B.plot_exp(plot_current,plot_ytrk,plot_ytrkerr)
B.pl.suptitle('LD2 Normalized Tracking Yield')
B.pl.title(plot_currname)
B.pl.ylabel('Normalized Tracking Yield: Y/(Q*tLT*TrkEff)')
B.pl.xlabel('Avg. Current [$\mu$A]')
for i, label in enumerate(LD2_runs):
    B.pl.annotate(label, (plot_current[i], plot_ytrk[i]),
                  textcoords="offset points", xytext=(0,10), ha='left')
# linefit = B.linefit(np.array(LD2_Acurrent),
#                     np.array(LD2_trk_yield_BCM4A),
#                     np.array(LD2_trk_yield_BCM4A_err),plot_fit=False)
# linefit.plot(label = f'Slope = {linefit.slope:0.2e}$\pm${linefit.sigma_s:0.2e}\n'+\
#                      f'Offset = {linefit.offset:0.2f}$\pm${linefit.sigma_o:0.2f}\n'+\
#                      '$\chi_{red}^2$' + f' = {linefit.chi_red:0.2f}')
# B.pl.legend()
    
B.pl.figure()
B.pl.suptitle('LD2 Fractional Normalized Tracking Yield')
B.plot_exp(plot_current,plot_ytrk/plot_ytrk[0],
           plot_ytrkerr/plot_ytrk[0],label='LD2_data_trk')
B.pl.ylim(0.9,1.1)
linefit_trk = B.linefit(plot_current,
                    plot_ytrk/plot_ytrk[0],
                    plot_ytrkerr/plot_ytrk[0],plot_fit=False)
linefit_trk.plot(label = f'Slope = {linefit_trk.slope:0.2e}$\pm${linefit_trk.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_trk.offset:0.2f}$\pm${linefit_trk.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_trk.chi_red:0.2f}')
B.pl.ylabel('Normalized Tracking Yield: Y/(Q*tLT*TrkEff)')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.legend()

#%% LD2 NORMALIZED YIELD PLOT
B.pl.figure(figsize=(8,10),layout='constrained')
# B.pl.suptitle('LD2 Normalized Yield')

B.pl.subplot(3,1,1)
B.plot_exp(plot_current,plot_yscl,plot_ysclerr,color=cmp['Paired'].colors[5],
          markersize=10,capsize=8)
for i, label in enumerate(LD2_runs):
    B.pl.annotate(label, (plot_current[i], plot_yscl[i]),
                  textcoords="offset points", xytext=(0,5), ha='left')
B.pl.ylim(9.2e6,9.9e6)
B.pl.xlim(5,75)
B.pl.title(plot_currname)
B.pl.text(8,9.25e6,'LD2 Scaler Yield',fontsize=24)

B.pl.subplot(3,1,2)
B.plot_exp(plot_current,plot_yntrk,plot_yntrkerr,color=cmp['Paired'].colors[5],
           markersize=10,capsize=8)
# B.pl.ylim(20300,22000)
B.pl.xlim(5,75)
B.pl.ylabel('Normalized Yield')
B.pl.text(8,61500,'LD2 No Tracking Yield',fontsize=24)

B.pl.subplot(3,1,3)
B.plot_exp(plot_current,plot_ytrk,plot_ytrkerr,color=cmp['Paired'].colors[5],
           markersize=10,capsize=8)
# B.pl.ylim(13000,14500)
B.pl.xlim(5,75)
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.text(8,40150,'LD2 Tracking Yield',fontsize=24)

# B.pl.savefig('target boiling/plots/LD2_norm_yield_all.png')

#%% plot frac norm yields 
B.pl.figure(figsize=(20,9),layout='constrained')
B.pl.suptitle('Fractional Normalized Yield')

B.pl.subplot(1,3,1)
B.plot_exp(np.array(C12_avg_current),C12_trk_yield_BCM4A/C12_trk_yield_BCM4A[0],
           C12_trk_yield_BCM4A_err/C12_trk_yield_BCM4A[0],label='Carbon-12')
linefit_trk = B.linefit(np.array(C12_avg_current),C12_trk_yield_BCM4A/C12_trk_yield_BCM4A[0],
           C12_trk_yield_BCM4A_err/C12_trk_yield_BCM4A[0],plot_fit=False)
linefit_trk.plot(label = f'Slope = {linefit_trk.slope:0.2e}$\pm${linefit_trk.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_trk.offset:0.2f}$\pm${linefit_trk.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_trk.chi_red:0.2f}')

B.plot_exp(np.array(LD2_avg_current),LD2_trk_yield_BCM4A/LD2_trk_yield_BCM4A[0],
           LD2_trk_yield_BCM4A_err/LD2_trk_yield_BCM4A[0],label='$LD_2$')

linefit_trk = B.linefit(np.array(LD2_avg_current),LD2_trk_yield_BCM4A/LD2_trk_yield_BCM4A[0],
           LD2_trk_yield_BCM4A_err/LD2_trk_yield_BCM4A[0],plot_fit=False)
linefit_trk.plot(label = f'Slope = {linefit_trk.slope:0.2e}$\pm${linefit_trk.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_trk.offset:0.2f}$\pm${linefit_trk.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_trk.chi_red:0.2f}')
B.pl.ylabel('Normalized Tracking Yield: Y/(Q*tLT*TrkEff)')
B.pl.legend()
B.pl.ylim(0.94,1.1)
B.pl.xlim(0,80)
B.pl.text(5,1.03,'Fit Function: $\\frac{Y_{norm}}{Y_0} = \\frac{m}{Y_0}I_{avg} +1$',
          fontsize=18)

B.pl.text(5,0.95,'Tracking Yield', fontsize=24)

B.pl.subplot(1,3,2)
B.plot_exp(np.array(C12_avg_current),C12_notrk_yield_BCM4A/C12_notrk_yield_BCM4A[0],
           C12_notrk_yield_BCM4A_err/C12_notrk_yield_BCM4A[0],label='Carbon-12')
linefit_notrk = B.linefit(np.array(C12_avg_current),C12_notrk_yield_BCM4A/C12_notrk_yield_BCM4A[0],
           C12_notrk_yield_BCM4A_err/C12_notrk_yield_BCM4A[0],plot_fit=False)
linefit_notrk.plot(label = f'Slope = {linefit_notrk.slope:0.2e}$\pm${linefit_notrk.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_notrk.offset:0.2f}$\pm${linefit_notrk.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_notrk.chi_red:0.2f}')

B.plot_exp(np.array(LD2_avg_current),LD2_notrk_yield_BCM4A/LD2_notrk_yield_BCM4A[0],
           LD2_notrk_yield_BCM4A_err/LD2_notrk_yield_BCM4A[0],label='$LD_2$')

linefit_notrk = B.linefit(np.array(LD2_avg_current),LD2_notrk_yield_BCM4A/LD2_notrk_yield_BCM4A[0],
           LD2_notrk_yield_BCM4A_err/LD2_notrk_yield_BCM4A[0],plot_fit=False)
linefit_notrk.plot(label = f'Slope = {linefit_notrk.slope:0.2e}$\pm${linefit_notrk.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_notrk.offset:0.2f}$\pm${linefit_notrk.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_notrk.chi_red:0.2f}')
B.pl.legend()
B.pl.xlim(0,80)
B.pl.ylim(0.94,1.1)
B.pl.ylabel('Normalized No Tracking Yield: Y/(Q*tLT)')
B.pl.xlabel('Avg. Current [$\mu$A]')
B.pl.text(5,0.95,'No Tracking Yield', fontsize=24)

B.pl.subplot(1,3,3)
B.plot_exp(np.array(C12_avg_current),C12_scl_yield_BCM4A/C12_scl_yield_BCM4A[0],
           C12_scl_yield_BCM4A_err/C12_scl_yield_BCM4A[0],label='Carbon-12')
linefit_scl = B.linefit(np.array(C12_avg_current),C12_scl_yield_BCM4A/C12_scl_yield_BCM4A[0],
           C12_scl_yield_BCM4A_err/C12_scl_yield_BCM4A[0],plot_fit=False)
linefit_scl.plot(label = f'Slope = {linefit_scl.slope:0.2e}$\pm${linefit_scl.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_scl.offset:0.2f}$\pm${linefit_scl.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_scl.chi_red:0.2f}')

B.plot_exp(np.array(LD2_avg_current),LD2_scl_yield_BCM4A/LD2_scl_yield_BCM4A[0],
           LD2_scl_yield_BCM4A_err/LD2_scl_yield_BCM4A[0],label='$LD_2$')

linefit_scl = B.linefit(np.array(LD2_avg_current),LD2_scl_yield_BCM4A/LD2_scl_yield_BCM4A[0],
           LD2_scl_yield_BCM4A_err/LD2_scl_yield_BCM4A[0],plot_fit=False)
linefit_scl.plot(label = f'Slope = {linefit_scl.slope:0.2e}$\pm${linefit_scl.sigma_s:0.2e}\n'+\
                     f'Offset = {linefit_scl.offset:0.2f}$\pm${linefit_scl.sigma_o:0.2f}\n'+\
                     '$\chi_{red}^2$' + f' = {linefit_scl.chi_red:0.2f}')
B.pl.legend()
B.pl.xlim(0,80)
B.pl.ylim(0.94,1.1)
B.pl.ylabel('Normalized Scaler Yield: Y/Q')
B.pl.text(5,0.95,'Scaler Yield', fontsize=24)
#%% write target boiling correction factors into db
runlist = D.get_list(db.retrieve('deuteron_db.db', 'run', 'RUN_LIST_UPDATED',
                                 where = 'target = \'LD2\' AND kin_study = \'deep\''))
tgtBoil_corr = []
tgtBoil_corr_err = []
## TBD: need to check BCM calib for uncertainty in charge
dI = 0.02 # from Carlos' thesis estimate in charge unc
m0 = linefit_trk.slope
dm0 = linefit_trk.sigma_s
for r in runlist:    
    I_avg = D.get_list(db.retrieve('deuteron_db.db', 'BCM4A_current', 
                                   'RUN_LIST_UPDATED', where = f'run=\'{r}\''))
    e_tb = 1 + m0*I_avg
    de_tb = np.sqrt(I_avg*I_avg*dm0*dm0 + m0*m0*dI*dI)
    
    
    tgtBoil_corr.append(e_tb)
    tgtBoil_corr_err.append(de_tb)
    db.update_row('../deuteron_db/deuteron_db.db', 'RUN_LIST_UPDATED',
                  ['tgtBoil_corr','tgtBoil_corr_err'],
                  [f'{e_tb:.4f}',f'{de_tb:.4f}'], 
                  where=f'run = {r}')
    
#%% write target boiling correction factors for heep based on carlos' slope
runlist = D.get_list(db.retrieve('deuteron_db.db', 'run', 'RUN_LIST_UPDATED',
                                 where = 'target = \'LH2\' AND kin_study = \'heep_coin\''))
tgtBoil_corr = []
tgtBoil_corr_err = []
## TBD: need to check BCM calib for uncertainty in charge
dI = 0.02 # from Carlos' thesis estimate in charge unc
m0 = -6.34e-4
dm0 = 6.2e-5
for r in runlist:    
    I_avg = D.get_list(db.retrieve('deuteron_db.db', 'BCM4A_current', 
                                   'RUN_LIST_UPDATED', where = f'run=\'{r}\''))
    e_tb = 1 + m0*I_avg
    de_tb = np.sqrt(I_avg*I_avg*dm0*dm0 + m0*m0*dI*dI)
    
    
    tgtBoil_corr.append(e_tb)
    tgtBoil_corr_err.append(de_tb)
    # print(e_tb,de_tb)
    db.update_row('../deuteron_db/deuteron_db.db', 'RUN_LIST_UPDATED',
                  ['tgtBoil_corr','tgtBoil_corr_err'],
                  [f'{e_tb:.4f}',f'{de_tb:.4f}'], 
                  where=f'run = {r}')

#%%
data =np.array(db.retrieve('deuteron_db.db', 'run,tgtBoil_corr,tgtBoil_corr_err', 'RUN_LIST_UPDATED',
                                 where = 'target = \'LD2\' AND kin_study = \'deep\''))
runs = data[:,0]
tgtb = data[:,1]
tgtbe = data[:,2]

#%%
B.pl.figure()
p = B.plot_exp(runs,tgtb,tgtbe)
tgtmean = np.mean(tgtb)
B.pl.hlines([tgtmean],xmin=20825,xmax=21115,
            color=['red'],linestyle='--')
B.pl.text(20830, tgtmean-0.004, f'{tgtmean:0.4f}',color='r')
B.pl.ylim(0.95,1.0)
B.pl.xlabel('Run #')
B.pl.ylabel('Target Boiling Factor')

    
#%% plot current cuts

for r in C12_runs:
    curr1 = C12_BCM1_cut.current[r]
    cmin1 = C12_BCM1_cut.c_cut_min[r]
    cmax1 = C12_BCM1_cut.c_cut_max[r]
    
    Acurr = C12_BCM4A_cut.current[r]
    Acmin = C12_BCM4A_cut.c_cut_min[r]
    Acmax = C12_BCM4A_cut.c_cut_max[r]
    
    Bcurr = C12_BCM4B_cut.current[r]
    Bcmin = C12_BCM4B_cut.c_cut_min[r]
    Bcmax = C12_BCM4B_cut.c_cut_max[r]
    
    Ccurr = C12_BCM4C_cut.current[r]
    Ccmin = C12_BCM4C_cut.c_cut_min[r]
    Ccmax = C12_BCM4C_cut.c_cut_max[r]
    
    h1 = B.histo(curr1,range=(0,70),bins=100,title='BCM1',xlabel='',ylabel='')
    hA = B.histo(Acurr,range=(0,70),bins=100,title='BCM4A',xlabel='',ylabel='')
    hB = B.histo(Bcurr,range=(0,70),bins=100,title='BCM4B',xlabel='',ylabel='')
    hC = B.histo(Ccurr,range=(0,70),bins=100,title='BCM4C',xlabel='',ylabel='')
    
    hcurrs = [h1,hA,hB,hC]
    hmins = [cmin1,Acmin,Bcmin,Ccmin]
    hmaxs = [cmax1,Acmax,Bcmax,Ccmax]
    
    B.pl.figure(figsize=(15,5))
    B.pl.suptitle(f'BCM Scaler Current\nRun {r}')
    for i in range(4):  
        B.pl.subplot(1,4,i+1)
        hcurrs[i].plot()
        B.pl.vlines([hmins[i],hmaxs[i]],ymin=0,ymax=max(hcurrs[i].bin_content),
                    linestyles='--',color='black')
    B.pl.savefig(f'target boiling/plots/current cuts/C12_run_{r}')
        