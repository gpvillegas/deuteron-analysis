#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:59:11 2026

@author: gvill
"""

import data_init as D
import LT.box as B

selection = ['H.gtr.th','H.gtr.ph','P.gtr.th','P.gtr.ph',]
br_sel_SIMC = ['h_xptar','h_yptar','e_xptar','e_yptar',
               'h_xptar_v','h_yptar_v','e_xptar_v','e_yptar_v',
               'Weight','Normfac','probabs']

DATA_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/pass_3/"
SIMC_baseDIR = "/media/gvill/Gema's T7/ROOTfiles/worksim/deep_wide_ang_acc/"

data = D.DATA_INIT(data_type='deut23_data', run = 20871,
                   select_branches={'T':selection},
                   ROOTfiles_path=DATA_baseDIR)
simc_norad = D.DATA_INIT(data_type='SIMC', kin_study='deep', setting='pm_120',
                   select_branches={'SNT':br_sel_SIMC},simc_type='jmlfsi_norad',
                   SIMC_ROOTfiles_path=SIMC_baseDIR)
simc_rad = D.DATA_INIT(data_type='SIMC', kin_study='deep', setting='pm_120',
                   select_branches={'SNT':br_sel_SIMC},simc_type='jmlfsi_rad',
                   SIMC_ROOTfiles_path=SIMC_baseDIR)
#%%
data_norm = 1/(D.get_charge_norm(20871,curr='BCM1')*D.get_eff_norm(20871))
simc_norad_weight = D.calc_weights(simc_norad.Branches)
simc_rad_weight = D.calc_weights(simc_rad.Branches)

#%%
# hxptar
data_hxptar = data_norm*B.histo(data.Branches['H.gtr.th'],range=(-0.15,0.15),bins=100)
simc_norad_hxptar = B.histo(simc_norad.Branches['h_xptar'],range=(-0.15,0.15),bins=100,
                      weights=simc_norad_weight)
simc_norad_hxptar_v = B.histo(simc_norad.Branches['h_xptar_v'],range=(-0.15,0.15),bins=100,
                        weights=simc_norad_weight)
simc_rad_hxptar = B.histo(simc_rad.Branches['h_xptar'],range=(-0.15,0.15),bins=100,
                      weights=simc_rad_weight)
simc_rad_hxptar_v = B.histo(simc_rad.Branches['h_xptar_v'],range=(-0.15,0.15),bins=100,
                        weights=simc_rad_weight)

# hyptar
data_hyptar = data_norm*B.histo(data.Branches['H.gtr.ph'],range=(-0.05,0.05),bins=100)
simc_norad_hyptar = B.histo(simc_norad.Branches['h_yptar'],range=(-0.05,0.05),bins=100,
                      weights=simc_norad_weight)
simc_norad_hyptar_v = B.histo(simc_norad.Branches['h_yptar_v'],range=(-0.05,0.05),bins=100,
                        weights=simc_norad_weight)
simc_rad_hyptar = B.histo(simc_rad.Branches['h_yptar'],range=(-0.05,0.05),bins=100,
                      weights=simc_rad_weight)
simc_rad_hyptar_v = B.histo(simc_rad.Branches['h_yptar_v'],range=(-0.05,0.05),bins=100,
                        weights=simc_rad_weight)

# pxptar
data_pxptar = data_norm*B.histo(data.Branches['P.gtr.th'],range=(-0.06,0.06),bins=100)
simc_norad_pxptar = B.histo(simc_norad.Branches['e_xptar'],range=(-0.06,0.06),bins=100,
                      weights=simc_norad_weight)
simc_norad_pxptar_v = B.histo(simc_norad.Branches['e_xptar_v'],range=(-0.06,0.06),bins=100,
                        weights=simc_norad_weight)
simc_rad_pxptar = B.histo(simc_rad.Branches['e_xptar'],range=(-0.06,0.06),bins=100,
                      weights=simc_rad_weight)
simc_rad_pxptar_v = B.histo(simc_rad.Branches['e_xptar_v'],range=(-0.06,0.06),bins=100,
                        weights=simc_rad_weight)

# pyptar
data_pyptar = data_norm*B.histo(data.Branches['P.gtr.ph'],range=(-0.04,0.04),bins=100)
simc_norad_pyptar = B.histo(simc_norad.Branches['e_yptar'],range=(-0.04,0.04),bins=100,
                      weights=simc_norad_weight)
simc_norad_pyptar_v = B.histo(simc_norad.Branches['e_yptar_v'],range=(-0.04,0.04),bins=100,
                        weights=simc_norad_weight)
simc_rad_pyptar = B.histo(simc_rad.Branches['e_yptar'],range=(-0.04,0.04),bins=100,
                      weights=simc_rad_weight)
simc_rad_pyptar_v = B.histo(simc_rad.Branches['e_yptar_v'],range=(-0.04,0.04),bins=100,
                        weights=simc_rad_weight)

#%%
B.pl.figure(figsize=(15,15))
B.pl.subplot(2,2,1)
data_hxptar.plot(label='data')
simc_norad_hxptar.plot(filled=False,color='r',label='simc_norad')
simc_norad_hxptar_v.plot(filled=False,color='orange',label='vertex_simc_norad')
simc_rad_hxptar.plot(filled=False,color='blue',label='simc_rad')
simc_rad_hxptar_v.plot(filled=False,color='cyan',label='vertex_simc_rad')

B.pl.title('h_xptar')
B.pl.xlabel('')
B.pl.ylabel('')
B.pl.ylim(0,35)
B.pl.legend()

B.pl.subplot(2,2,2)
data_hyptar.plot(label='data')
simc_norad_hyptar.plot(filled=False,color='r',label='simc_norad')
simc_norad_hyptar_v.plot(filled=False,color='orange',label='vertex_simc_norad')
simc_rad_hyptar.plot(filled=False,color='blue',label='simc_rad')
simc_rad_hyptar_v.plot(filled=False,color='cyan',label='vertex_simc_rad')

B.pl.title('h_yptar')
B.pl.xlabel('')
B.pl.ylabel('')
B.pl.ylim(0,38)
B.pl.legend()

B.pl.subplot(2,2,3)
data_pxptar.plot(label='data')
simc_norad_pxptar.plot(filled=False,color='r',label='simc_norad')
simc_norad_pxptar_v.plot(filled=False,color='orange',label='vertex_simc_norad')
simc_rad_pxptar.plot(filled=False,color='blue',label='simc_rad')
simc_rad_pxptar_v.plot(filled=False,color='cyan',label='vertex_simc_rad')

B.pl.title('p_xptar')
B.pl.xlabel('')
B.pl.ylabel('')
B.pl.ylim(0,40)
B.pl.legend()

B.pl.subplot(2,2,4)
data_pyptar.plot(label='data')
simc_norad_pyptar.plot(filled=False,color='r',label='simc_norad')
simc_norad_pyptar_v.plot(filled=False,color='orange',label='vertex_simc_norad')
simc_rad_pyptar.plot(filled=False,color='blue',label='simc_rad')
simc_rad_pyptar_v.plot(filled=False,color='cyan',label='vertex_simc_rad')

B.pl.title('p_yptar')
B.pl.xlabel('')
B.pl.ylabel('')
B.pl.ylim(0,45)
B.pl.legend()