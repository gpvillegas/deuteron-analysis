#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 09:58:08 2025

@author: gvill

hcana_dict is a function that returns a list of desired hcana leaf variables
categorized by type. The available types are:
    

"""
def hcana_dict(b=None):
    
    ###################
    # T TREE BRANCHES #
    ###################
    
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
    Phodo_ana_branches = ['P.hod.beta']
    
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
    Hhodo_ana_branches = ['H.hod.beta']
    
    # DRIFT CHAMBER branches
    Pdc_branches = ['P.dc.x_fp','P.dc.y_fp','P.dc.xp_fp','P.dc.yp_fp'] 
    
    Hdc_branches = ['H.dc.x_fp','H.dc.y_fp','H.dc.xp_fp','H.dc.yp_fp']
    
    # CALORIMETER branches
    Pcal_branches = ['P.cal.etottracknorm', 'P.cal.etotnorm']
    
    Hcal_branches = ['H.cal.etottracknorm', 'H.cal.etotnorm']
    
    # CERENKOV branches
    Phgcer_branches = ['P.hgcer.npeSum']
    
    Pngcer_branches = ['P.ngcer.npeSum']
    
    Hcer_branches = ['H.cer.npeSum']
    
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
    
    kin_ana_branches = ['H.kin.secondary.Erecoil','H.kin.secondary.Mrecoil',
                        'H.kin.secondary.emiss','H.kin.secondary.emiss_nuc',
                        'H.kin.secondary.pmiss','H.kin.secondary.pmiss_x',
                        'H.kin.secondary.pmiss_y','H.kin.secondary.pmiss_z',
                        'H.kin.secondary.tb','H.kin.secondary.th_bq',
                        'H.kin.secondary.th_xq','H.kin.secondary.tx',
                        'H.kin.secondary.xangle','H.kin.secondary.Prec_x',
                        'H.kin.secondary.Prec_y','H.kin.secondary.Prec_z',
                        'H.kin.secondary.ph_bq','H.kin.secondary.ph_xq',
                        'P.kin.primary.scat_ang_rad',
                        'P.kin.primary.Q2','P.kin.primary.W',
                        'P.kin.primary.x_bj','P.kin.primary.nu',
                        'P.kin.primary.q3m','P.kin.primary.th_q']
    
    # GOLDEN TRACK branches
    
    gtr_branches= ['H.dc.gtrack_nsp','H.gtr.beta','H.gtr.dp','H.gtr.index',
                   'H.gtr.ok','H.gtr.p','H.gtr.ph','H.gtr.px','H.gtr.py',
                   'H.gtr.pz','H.gtr.th','H.gtr.x','H.gtr.y','P.dc.gtrack_nsp',
                   'P.gtr.beta','P.gtr.dp','P.gtr.index','P.gtr.ok','P.gtr.p',
                   'P.gtr.ph','P.gtr.px','P.gtr.py','P.gtr.pz','P.gtr.th',
                   'P.gtr.x','P.gtr.y']
    
    gtr_ana_branches = ['P.gtr.dp','P.gtr.p','P.gtr.th','P.gtr.ph','P.gtr.y',
                        'H.gtr.dp','H.gtr.p','H.gtr.y','H.gtr.th','H.gtr.ph']
    
    # REACTION VERTEX branches
    
    react_branches= ['H.react.x','H.react.y','H.react.z',
                     'P.react.x','P.react.y','P.react.z']
    
    # sieve branches
    sieve_branches = ['H.extcor.xsieve', 'H.extcor.ysieve',
                      'P.extcor.xsieve', 'P.extcor.ysieve']
    
    # coin time branches
    ctime_branches = ['CTime.epCoinTime_ROC2']
    
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
    
    
    #####################
    # TSP TREE BRANCHES #
    #####################
    
    curr_cut_branches = ['evNumber','P.BCM4A.scalerCharge',
                         'P.BCM4A.scalerChargeCut','P.BCM4A.scalerCurrent',
                         'P.BCM4B.scalerCharge','P.BCM4B.scalerCurrent',
                         'P.BCM4C.scalerCharge','P.BCM4C.scalerCurrent',
                         'P.1MHz.scalerTime']
    
    
    T_branch_types= {'phodo':Phodo_branches,
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
                     'optics':optics_branches,
                     'kin_ana':kin_ana_branches,
                     'phodo_ana':Phodo_ana_branches,
                     'hhodo_ana':Hhodo_ana_branches,
                     'gtr_ana':gtr_ana_branches,
                     'sieve':sieve_branches,
                     'ctime':ctime_branches
                     }

    TSP_branch_types = {'curr_cut':curr_cut_branches}
    
    analysis_group = {'T' : T_branch_types['kin_ana']+\
                            T_branch_types['phodo_ana']+\
                            T_branch_types['hhodo_ana']+\
                            T_branch_types['hcer']+\
                            T_branch_types['pngcer']+\
                            T_branch_types['phgcer']+\
                            T_branch_types['hcal']+\
                            T_branch_types['pcal']+\
                            T_branch_types['gtr_ana']+\
                            T_branch_types['react']+\
                            T_branch_types['sieve']+\
                            T_branch_types['ctime']+\
                            T_branch_types['hdc']+\
                            T_branch_types['pdc'],
                      'TSP' : TSP_branch_types['curr_cut']}
    
    all_T_branches = []
    for br in T_branch_types:
        all_T_branches = all_T_branches + T_branch_types[br]

    if b is list:
        branches_sel = []
        for name in b:
            branches_sel = branches_sel + T_branch_types[name]        
        return branches_sel
    elif b == 'analysis':
        return analysis_group
    else: 
        return all_T_branches
