#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 12:05:41 2025

@author: gvill
"""
import data_init as D
import hcana_dict as hc
import cut_handler as C
import LT.box as B

import numpy as np

#%%
class analysis_histos:
    def __init__(self,dataInit_obj,cuts,data_type = 'deut23_data', 
                 histo_type = 'all',react_type = 'deep', curr_type = 'BCM4A'): 
        d = dataInit_obj.Branches
        
        #define constants
        MD = 1.87561    # deuteron mass 
        MN = 0.939566   # neutron mass
        MP = 0.938272   # proton mass
        
        rtd = 180./np.pi
        
        if data_type == 'deut23_data':
            # detector
            charge = D.get_charge_norm(dataInit_obj.many)
            current = d[f'P.{curr_type}.scalerCurrent_cut']
            cTime = d['CTime.epCoinTime_ROC2_cut']
            
            hBeta = d['H.hod.beta'][cuts]
            hCer = d['H.cer.npeSum'][cuts]
            hCal = d['H.cal.etotnorm'][cuts]
            
            pBeta = d['P.hod.beta'][cuts]
            pNGCer = d['P.ngcer.npeSum'][cuts]
            pHGCer = d['P.hgcer.npeSum'][cuts]
            pCal = d['P.cal.etottracknorm'][cuts]
            
            # primary kinematics
            Q2 = d['P.kin.primary.Q2'][cuts]
            omega = d['P.kin.primary.nu'][cuts]
            W = d['P.kin.primary.W'][cuts]
            xbj = d['P.kin.primary.x_bj'][cuts]
            kf = d['P.gtr.p'][cuts]
            thq = d['P.kin.primary.th_q'][cuts]*rtd
            qvec = d['P.kin.primary.q3m'][cuts]
            the = d['P.kin.primary.scat_ang_rad'][cuts]*rtd
            
            # calculate invariant mass squared
            W2 = W**2
            
            # secondary kinematics
            Em = d['H.kin.secondary.emiss'][cuts]
            Em_nuc = d['H.kin.secondary.emiss_nuc'][cuts]
            Pm = d['H.kin.secondary.pmiss'][cuts]
            Pmx_lab = d['H.kin.secondary.Prec_x'][cuts]
            Pmy_lab = d['H.kin.secondary.Prec_y'][cuts]
            Pmz_lab = d['H.kin.secondary.Prec_z'][cuts]
            Pmx_q = d['H.kin.secondary.pmiss_x'][cuts]
            Pmy_q = d['H.kin.secondary.pmiss_y'][cuts]
            Pmz_q = d['H.kin.secondary.pmiss_z'][cuts]
            Kp = d['H.kin.secondary.tx'][cuts]
            Kn = d['H.kin.secondary.tb'][cuts]
            pf = d['H.gtr.p'][cuts]
            thp = d['H.kin.secondary.xangle'][cuts]*rtd - the
            th_pq = d['H.kin.secondary.th_xq'][cuts]*rtd
            th_nq = d['H.kin.secondary.th_bq'][cuts]*rtd
            cosph_pq = np.cos(d['H.kin.secondary.ph_xq'][cuts])
            cosph_nq = np.cos(d['H.kin.secondary.ph_bq'][cuts])
            sinph_pq = np.sin(d['H.kin.secondary.ph_xq'][cuts])
            sinph_nq = np.sin(d['H.kin.secondary.ph_bq'][cuts])
            
            # calculate proton and neutron energy
            Ep = Kp + MP
            En = Kn + MN
            
            # calculate missing mass and missing mass squared
            if react_type == 'deep':
                MM2 = (omega + MD - Ep)**2 - Pm**2
            elif react_type == 'heep':
                MM2 = Em**2 - Pm**2
            
            MM = np.sqrt(MM2)
            
            # light front variables
            PmPar = Pm*np.cos(th_nq)
            PmPerp = np.sqrt(Pm**2 - PmPar**2)
            
            # LF neutron momentum fraction
            alpha_n = 2*(En - PmPar)/MD
            
            # LF momentum fraction of struck nucleon (proton)
            alpha = 2. - alpha_n  
            
            # target reconstruction (lab)
            htar_x = d['H.react.x'][cuts]
            htar_y = d['H.react.y'][cuts]
            htar_z = d['H.react.z'][cuts]
            
            ptar_x = d['P.react.x'][cuts]
            ptar_y = d['P.react.y'][cuts]
            ptar_z = d['P.react.z'][cuts]
            
            ztar_diff = htar_z - ptar_z
            
            # hadron arm reconstructed quantities
            hytar = d['H.gtr.y'][cuts]
            hxptar = d['H.gtr.th'][cuts]
            hyptar = d['H.gtr.ph'][cuts]
            hdp = d['H.gtr.dp'][cuts]
            
            # electron arm reconstructed quantities
            pytar = d['P.gtr.y'][cuts]
            pxptar = d['P.gtr.th'][cuts]
            pyptar = d['P.gtr.ph'][cuts]
            pdp = d['P.gtr.dp'][cuts]
            
            # hadron arm focal plane quantities
            hxfp = d['H.dc.x_fp'][cuts]
            hyfp = d['H.dc.y_fp'][cuts]
            hxpfp = d['H.dc.xp_fp'][cuts]
            hypfp = d['H.dc.yp_fp'][cuts]
            
            # electron arm focal plane quantities
            pxfp = d['P.dc.x_fp'][cuts]
            pyfp = d['P.dc.y_fp'][cuts]
            pxpfp = d['P.dc.xp_fp'][cuts]
            pypfp = d['P.dc.yp_fp'][cuts]
            
            # HMS collimator quantities
            hXColl = d['H.extcor.xsieve']
            hYColl = d['H.extcor.ysieve']
            hXColl_cut = d['H.extcor.xsieve_cut']
            hYColl_cut = d['H.extcor.ysieve_cut']
            
            # SHMS collimator quantities
            pXColl = d['P.extcor.xsieve']
            pYColl = d['P.extcor.ysieve']
            pXColl_cut = d['P.extcor.xsieve_cut']
            pYColl_cut = d['P.extcor.ysieve_cut']
            
            
        self.histos_1D = {'Charge':[charge,(-0.5,2.5),3],
                          f'{curr_type} Scaler Current Cut':
                              [current,(0,70),100],
                          'Coincidence Time Corrected':[cTime,(-3,3),100],
                          'HMS Beta':[hBeta,(0,2),100],
                          'HMS Cerenkov Npe Sum':[hCer,(0,10),100],
                          'HMS Calorimeter EtotNorm':[hCal,(0,2),100],
                          'SHMS Beta':[pBeta,(0,2),100],
                          'SHMS NG Cerenkov Npe Sum':[pNGCer,(0,50),100],
                          'SHMS HG Cerenkov Npe Sum':[pHGCer,(0,20),100],
                          'SHMS Calorimeter EtotTrackNorm':[pCal,(0,2),100],
                          'Q2':[Q2,(0,7),100],
                          'Energy Transfer, $\omega$':[omega,(0,4),100],
                          'Invariant Mass, W':[W,(0,2.5),100],
                          'Invariant Mass Squared, W2':[W2,(0,2.5),100],
                          'Bjorken x':[xbj,(0,2),100],
                          'Final Electron Momentum, kf':[kf,(2,10),100],
                          'q vector angle (Lab)':[thq,(0,180),100],
                          'q vector magnitude':[qvec,(0,7),100],
                          'Electron Scattering Angle (Lab), $\\theta_e$':
                              [the,(0,180),100],
                          'Missing Energy: From Conservation':
                              [Em,(-0.1,1.0),100],
                          'Nuclear Missing Energy: Binding Energy':
                              [Em_nuc,(-0.1,1.0),100],
                          'Missing Mass':[MM,(0,2.),100],
                          'Missing Mass Squared':[MM2,(0,2.),100],
                          'Missing Momentum':[Pm,(-0.1,1.5),100],
                          'Pmiss X (Lab)':[Pmx_lab,(-1,1),100],
                          'Pmiss Y (Lab)':[Pmy_lab,(-1,1),100],
                          'Pmiss Z (Lab)':[Pmz_lab,(-1,1),100],
                          'Pmiss X (q)':[Pmx_q,(-1,1),100],
                          'Pmiss Y (q)':[Pmy_q,(-1,1),100],
                          'Pmiss Z (q)':[Pmz_q,(-1,1),100],
                          'Final Proton Energy':[Ep,(0,6),100],
                          'Final Proton Kinetic Energy':[Kp,(0,4),100],
                          'Final Proton Momentum':[pf,(0,4),100],
                          'Proton Scattering Angle (Lab), $\\theta_p$':
                              [thp,(0,180),100],
                          'Proton Scattering Angle (q), $\\theta_{pq}$':
                              [th_pq,(0,180),100],
                          'Final Neutron Energy':[En,(0,2),100],
                          'Final Neutron Kinetic Energy':[Kn,(0,2),100],
                          'Neutron Recoil Angle (q), $\\theta_{nq}$':[th_nq,(0,180),100],
                          'cos($\phi_{pq}$)':[cosph_pq,(0,1),100],
                          'cos($\phi_{nq}$)':[cosph_nq,(0,1),100],
                          'sin($\phi_{pq}$)':[sinph_pq,(0,1),100],
                          'sin($\phi_{nq}$)':[sinph_nq,(-1.5,1.5),100],
                          'Pmiss Parallel (q)':[PmPar,(-0.5,0.5),100],
                          'Pmiss Transverse (q)':[PmPerp,(-0.5,0.5),100],
                          'LF Neutron Momentum Fraction, $\\alpha_n$':
                              [alpha_n,(0,2.5),100],
                          'LF Proton Momentum Fraction, $\\alpha$':
                              [alpha,(0,2.5),100],
                          'HMS X Target':[htar_x,(-0.1,0.1),100],
                          'HMS Y Target':[htar_y,(-0.1,0.1),100],
                          'HMS Z Target':[htar_z,(-10,10),100],
                          'SHMS X Target':[ptar_x,(-0.1,0.1),100],
                          'SHMS Y Target':[ptar_y,(-0.1,0.1),100],
                          'SHMS Z Target':[ptar_z,(-10,10),100],
                          'Ztar Difference: HMS - SHMS':[ztar_diff,(-5,5),50],
                          'HMS $Y_{tar}$':[hytar,(-5,5),100],
                          'HMS $X\'_{tar}$':[hxptar,(-0.1,0.1),100],
                          'HMS $Y\'_{tar}$':[hyptar,(-0.1,0.1),100],
                          'HMS Momentum Acceptance, $\delta$':[hdp,(-10,10),100],
                          'SHMS $Y_{tar}$':[pytar,(-2,2),100],
                          'SHMS $X\'_{tar}$':[pxptar,(-0.1,0.1),100],
                          'SHMS $Y\'_{tar}$':[pyptar,(-0.1,0.1),100],
                          'SHMS Momentum Acceptance, $\delta$':[pdp,(-10,22),100],
                          'HMS $X_{fp}$':[hxfp,(-40,40),100],
                          'HMS $Y_{fp}$':[hyfp,(-20,20),100],
                          'HMS $X\'_{fp}$':[hxpfp,(-0.1,0.1),100],
                          'HMS $Y\'_{fp}$':[hypfp,(-0.1,0.1),100],
                          'SHMS $X_{fp}$':[pxfp,(-6,4),50],
                          'SHMS $Y_{fp}$':[pyfp,(-6,4),50],
                          'SHMS $X\'_{fp}$':[pxpfp,(-0.1,0.1),100],
                          'SHMS $Y\'_{fp}$':[pypfp,(-0.1,0.1),100],
                          'HMS $X_{coll}$ no Cuts':[hXColl,(-20,20),100],
                          'HMS $Y_{coll}$ no Cuts':[hYColl,(-15,15),100],
                          'HMS $X_{coll}$ Cut':[hXColl_cut,(-20,20),100],
                          'HMS $Y_{coll}$ Cut':[hYColl_cut,(-15,15),100],
                          'SHMS $X_{coll}$ no Cuts':[pXColl,(-20,20),100],
                          'SHMS $Y_{coll}$ no Cuts':[pYColl,(-20,20),100],
                          'SHMS $X_{coll}$ Cut':[pXColl_cut,(-20,20),100],
                          'SHMS $Y_{coll}$ Cut':[pYColl_cut,(-20,20),100],
                          }
        
        self.histos_2D = {'HMS XColl vs YColl':
                          [hYColl_cut,hXColl_cut,(-15,15),100,(-20,20),100],
                          'SHMS XColl vs YColl':
                          [pYColl_cut,pXColl_cut,(-20,20),100,(-20,20),100],
                          'HMS $X_{fp}$ vs $Y_{fp}$':
                          [hyfp,hxfp,(-20,20),100,(-40,40),100],
                          'SHMS $X_{fp}$ vs $Y_{fp}$':
                          [pyfp,pxfp,(-6,4),50,(-6,4),50],
                          '$E_{miss}$ vs $P_{miss}$':
                          [Pm,Em,(-0.1,1.5),100,(-0.1,0.2),100],
                          'Nuclear E_{miss} vs P_{miss}':
                          [Pm,Em_nuc,(-0.1,1.5),100,(-0.1,0.2),100],
                          'Q2 vs W':
                          [W,Q2,(0,2.5),100,(0,7),100],
                          'Q2 vs Pm':
                          [Pm,Q2,(-0.1,1.5),100,(0,7),100], 
                          'Pm vs $\\theta_{nq}$':
                          [Pm,th_nq,(-0.1,1.5),64,(0,180),60],
                          'Pt vs Alpha':
                          [alpha,PmPerp,(-0.5,0.5),100,(0,2.5),100]
                          }    

#%%

r = 20871
branches_sel = hc.hcana_dict('analysis')
t_sel = ['T','TSP']
DATA_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

RUN = D.DATA_INIT(data_type='deut23_data',run=r,
                         select_branches=branches_sel,
                         select_trees=t_sel,
                         ROOTfiles_path= DATA_DIR)

#%%
cuts_list = C.event_selection_cuts

hcoll_cut = C.coll_cut(RUN,spec='HMS')
hcoll_cut_arrays = hcoll_cut()

pcoll_cut = C.coll_cut(RUN,spec='SHMS')
pcoll_cut_arrays = pcoll_cut()

curr_type ='BCM4A'
curr_cut = C.current_cut(RUN,current=curr_type)
curr_cut_arrays = curr_cut()

zt_cut = C.ztar_cut(RUN)
zt_cut_arrays = zt_cut()

ct_cut = C.CTime_cut(RUN)
ct_cut_arrays = ct_cut() 

cuts_to_apply_list = []
for cut in cuts_list:
    br = RUN.Branches[C.HCANA_names[cut.name]]
    cut_array = cut(br)
    cut.stats()
    
    cuts_to_apply_list.append(cut_array)
    
# add special cut arrays to cuts_to_apply_list
cuts_to_apply_list.append(hcoll_cut_arrays)
cuts_to_apply_list.append(pcoll_cut_arrays)
cuts_to_apply_list.append(curr_cut_arrays)
cuts_to_apply_list.append(zt_cut_arrays)
cuts_to_apply_list.append(ct_cut_arrays)
    

all_cuts = cuts_to_apply_list[0]
for arr in cuts_to_apply_list:    
    all_cuts = all_cuts & arr  

# now apply the cuts and store the cut arrays back in the DATA_INIT obj
  
RUN.Branches.update({'H.extcor.xsieve_cut':hcoll_cut.xcoll[hcoll_cut_arrays],
                     'H.extcor.ysieve_cut':hcoll_cut.ycoll[hcoll_cut_arrays],
                     'P.extcor.xsieve_cut':pcoll_cut.xcoll[pcoll_cut_arrays],
                     'P.extcor.ysieve_cut':pcoll_cut.ycoll[pcoll_cut_arrays], 
                     'CTime.epCoinTime_ROC2_cut':ct_cut.CTime_corr,
                     f'P.{curr_type}.scalerCurrent_cut':
                         curr_cut.current[curr_cut.cut_array]
                     })
#%%
my_histos = analysis_histos(RUN, all_cuts)

for h in my_histos.histos_1D:
    print(h)
    if h == 'Charge':
        continue
        his = B.histo([],range=my_histos.histos_1D[h][1],
                      bins=my_histos.histos_1D[h][2])
        his.bin_content[1] = my_histos.histos_1D[h][0]
        
        B.pl.figure()
        his.plot()
        B.pl.title(h)
        B.pl.xlabel('')
        B.pl.ylabel('')
    else:
        # continue
        his = B.histo(my_histos.histos_1D[h][0],
                      range=my_histos.histos_1D[h][1],
                      bins=my_histos.histos_1D[h][2])
        B.pl.figure()
        his.plot()
        B.pl.title(h)
        B.pl.xlabel('')
        B.pl.ylabel('')