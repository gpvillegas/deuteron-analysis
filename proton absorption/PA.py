#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:48:32 2026

@author: gvill
"""

import data_init as D
import cut_handler as C
import numpy as np
import LT.box as B
import database_operations as db
from matplotlib import colormaps as cmp

#Set default font to times new roman
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14
}
B.pl.rc('font', **font)

T_branches = ['P.gtr.dp','H.gtr.dp','H.kin.secondary.emiss','P.gtr.th',
              'P.gtr.ph','H.gtr.th','H.gtr.ph','P.react.z','P.kin.primary.W',
              'P.cal.etottracknorm']

coin_run = 20851
singles_run = 20856
Al_singles_run = 20857 

ROOTfiles_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_2/"

coinData = D.DATA_INIT(data_type='deut23_data',run=coin_run,
                         select_branches={'T':T_branches},
                         ROOTfiles_path=ROOTfiles_DIR)
singData = D.DATA_INIT(data_type='deut23_data',run=singles_run,
                         select_branches={'T':T_branches+['T.coin.hTRIG1_ROC1_tdcTimeRaw']},
                         ROOTfiles_path=ROOTfiles_DIR)
AlData = D.DATA_INIT(data_type='deut23_data',run=Al_singles_run,
                         select_branches={'T':T_branches},
                         ROOTfiles_path=ROOTfiles_DIR)

#%% determine delta acceptance region for elastic events
# create missing energy cut
Em_cut = C.WCUT(-0.03,0.03,name='elastics_cut')
Em = coinData.Branches['H.kin.secondary.emiss']

c_Em = Em_cut(Em) # cut boolean array

# missing energy histograms
em_h = B.histo(Em,range=(-0.03,0.1),bins=100)
em_cut_h = B.histo(Em[c_Em],range=(-0.03,0.1),bins=100)

## plot missing energy
B.pl.figure()
em_h.plot(filled=False)
em_cut_h.plot()
B.pl.vlines(0.03,ymin=0,ymax=4300,linestyles='--',colors='black')
B.pl.title('Missing Energy')
B.pl.xlabel('[GeV]')
B.pl.ylabel('')

# delta acceptance histograms
hdelta = coinData.Branches['H.gtr.dp']
pdelta = coinData.Branches['P.gtr.dp']

delta_h = B.histo2d(hdelta[c_Em],pdelta[c_Em],range=[(-10,10),(-10,10)],bins=100,
                    title='Delta Acceptance SHMS vs. HMS',
                    xlabel='HMS $\delta$ [%]',ylabel='SHMS $\delta$ [%]')

# delta cut limits
hdelta_cut_lims = (-8,8)
pdelta_cut_lims = (-2.5,2.5)

# plot delta accptance
B.pl.figure()
delta_h.plot(colormap=cmp['viridis'])
B.pl.vlines([hdelta_cut_lims[0],hdelta_cut_lims[1]],ymin=-10,ymax=10,linestyles='--',colors='r')
B.pl.hlines([pdelta_cut_lims[0],pdelta_cut_lims[1]],xmin=-10,xmax=10,linestyles='--',colors='b')

#%% determine initial angular acceptance region within delta region for elastics

# create delta cuts
hdelta_cut = C.WCUT(hdelta_cut_lims[0],hdelta_cut_lims[1],name='hdelta_cut')
pdelta_cut = C.WCUT(pdelta_cut_lims[0],pdelta_cut_lims[1],name='pdelta_cut')

c_hdp = hdelta_cut(hdelta)
c_pdp = pdelta_cut(pdelta)

cuts = c_Em & c_hdp & c_pdp # combine cuts so far

# angular acceptance histograms
hxptar = coinData.Branches['H.gtr.th']
hyptar = coinData.Branches['H.gtr.ph']
pxptar = coinData.Branches['P.gtr.th']
pyptar = coinData.Branches['P.gtr.ph']

hms_ang_acc_h = B.histo2d(hxptar[cuts],hyptar[cuts],
                          range=[(-0.08,0.08),(-0.08,0.08)],bins=100,
                          title='HMS Angular Acceptance',
                          xlabel='$X\'_{tar}$ [rad]',
                          ylabel='$Y\'_{tar}$ [rad]')
shms_ang_acc_h = B.histo2d(pxptar[cuts],pyptar[cuts],
                          range=[(-0.03,0.03),(-0.03,0.03)],bins=100,
                          title='SHMS Angular Acceptance',
                          xlabel='$X\'_{tar}$ [rad]',
                          ylabel='$Y\'_{tar}$ [rad]')

# angular acceptance cut limits
hxptar_cut_lims = (-0.07,0.07)
hyptar_cut_lims = (-0.025,0.025)
pxptar_cut_lims = (-0.025,0.025)
pyptar_cut_lims = (-0.01,0.01)

# plot angular acceptance
B.pl.figure()
hms_ang_acc_h.plot(colormap=cmp['viridis'])
B.pl.vlines([hxptar_cut_lims[0],hxptar_cut_lims[1]],ymin=-0.08,ymax=0.08,linestyles='--',colors='r')
B.pl.hlines([hyptar_cut_lims[0],hyptar_cut_lims[1]],xmin=-0.08,xmax=0.08,linestyles='--',colors='b')

B.pl.figure()
shms_ang_acc_h.plot(colormap=cmp['viridis'])
B.pl.vlines([pxptar_cut_lims[0],pxptar_cut_lims[1]],ymin=-0.03,ymax=0.03,linestyles='--',colors='r')
B.pl.hlines([pyptar_cut_lims[0],pyptar_cut_lims[1]],xmin=-0.03,xmax=0.03,linestyles='--',colors='b')

#%% determine z target cut on shms using prev cuts

# create shms delta, xptar, yptar cuts for Al dummy run
pdelta_Al = AlData.Branches['P.gtr.dp']
pxptar_Al = AlData.Branches['P.gtr.th']
pyptar_Al = AlData.Branches['P.gtr.ph']

pxptar_cut = C.WCUT(pxptar_cut_lims[0],pxptar_cut_lims[1],name='pxptar_cut')
pyptar_cut = C.WCUT(pyptar_cut_lims[0],pyptar_cut_lims[1],name='pyptar_cut')

c_pdp = pdelta_cut(pdelta_Al)
c_pxptar = pxptar_cut(pxptar_Al)
c_pyptar = pyptar_cut(pyptar_Al)

cuts = c_pdp & c_pxptar & c_pyptar      # combine all cuts here

# ztar histogram from Al dummy data
pztar = AlData.Branches['P.react.z']

pztar_h = B.histo(pztar[cuts],range=(-10,10),bins=100,
                  title='SHMS $Z_{tar}$'+f' Al Dummy Run {Al_singles_run}',
                  xlabel='$Z_{tar}$ [cm]',
                  ylabel='')

# ztar cut limits
pztar_cut_lims = (-2.5,2.5)

# plot ztar
B.pl.figure()
pztar_h.plot()
B.pl.vlines([pztar_cut_lims[0],pztar_cut_lims[1]],ymin=0,ymax=310,linestyles='--',colors='r')

#%% look at singles angular acceptance for region determined form coin data

# recreate delta, xptar, yptar, ztar cuts for singles run
pdelta_s = singData.Branches['P.gtr.dp']
pxptar_s = singData.Branches['P.gtr.th']
pyptar_s = singData.Branches['P.gtr.ph']
pztar_s = singData.Branches['P.react.z']

pztar_cut = C.WCUT(pztar_cut_lims[0],pztar_cut_lims[1],name='pztar_cut')

c_pdp = pdelta_cut(pdelta_s)
c_pxptar = pxptar_cut(pxptar_s)
c_pyptar = pyptar_cut(pyptar_s)
c_pztar = pztar_cut(pztar_s)

cuts = c_pdp & c_pxptar & c_pyptar & c_pztar   # combine cuts here


# SHMS singles angular acceptance histogram
shms_ang_acc_s_h = B.histo2d(pxptar_s[cuts],pyptar_s[cuts],
                          range=[(-0.03,0.03),(-0.03,0.03)],bins=100,
                          title='SHMS Angular Acceptance (Coin Cuts)',
                          xlabel='$X\'_{tar}$ [rad]',
                          ylabel='$Y\'_{tar}$ [rad]')

# plot angular acceptance
B.pl.figure()
shms_ang_acc_s_h.plot(colormap=cmp['viridis'])


#%% prepare invariant mass cuts

# create cuts needed for should/ did cuts
W_s = singData.Branches['P.kin.primary.W']
hdp_s = singData.Branches['H.gtr.dp']
hTRIG1_s = singData.Branches['T.coin.hTRIG1_ROC1_tdcTimeRaw']
pcal_s = singData.Branches['P.cal.etottracknorm']

W_cut = C.WCUT(0.9,1.0,name='invariant_mass_cut')
TRIG_OFF_cut = C.VCUT(0,name='trigger_off_cut')
pcal_cut = C.WCUT(0.7,1.3,name='pcal_cut')

c_W = W_cut(W_s)
c_hdp = hdelta_cut(hdp_s)
c_hTRIG1 = ~TRIG_OFF_cut(hTRIG1_s)
c_pcal = pcal_cut(pcal_s)

# e_should: number of electrons detected in SHMS with elastics kinematics and 
#   acceptance that SHOULD have a proton coincident on HMS side 
e_should = c_pdp & c_pxptar & c_pyptar & c_pztar & c_W & c_pcal

# e_did: number of events that were actually detected on the HMS with a proton
#   in coincidence withing the HMS elastic acceptance region
e_did = e_should & c_hdp & c_hTRIG1

#%% get invariant mass did/should

# get tracking efficiencies to normalize yields
pTrkEff, hTrkEff, tLT, T1cpuLT = db.retrieve('deuteron_db.db','SHMS_TrkEff, HMS_TrkEff,tLT,'+\
                                    'T1_cpuLT', 
                               'RUN_LIST_UPDATED', 
                               where = f"run=\'{singles_run}\'")[0]
pTrkEff_coin, hTrkEff_coin = db.retrieve('deuteron_db.db','SHMS_TrkEff, HMS_TrkEff', 
                               'RUN_LIST_UPDATED', 
                               where = f"run=\'{coin_run}\'")[0]

# invariant mass histogram: no cuts
W_h = B.histo(W_s,range=(0.8,1.2),bins=80)

# invariant mass histogram: W cut
W_cut_h = B.histo(W_s[c_W],range=(0.8,1.2),bins=80)

# invariant mass histogram: e_did
W_did_h = B.histo(W_s[e_did],range=(0.8,1.2),bins=30)
W_did, W_did_e = W_did_h.sum()

# invariant mass histogram: e_should
W_should_h = B.histo(W_s[e_should],range=(0.8,1.2),bins=30)
W_should, W_should_e = W_should_h.sum()

#%% get proton absorption factor

# proton absorption factor, normalized by tracking eff
e_pAbs = W_did/W_should

print('Proton Absorption Factor: ',e_pAbs)

# plot invariant mass
B.pl.figure()
W_h.plot(filled=False)
W_cut_h.plot()
B.pl.vlines([0.9,1.0],0,1000,linestyles='--',colors='r')
B.pl.title('Invariant Mass')
B.pl.xlabel('W [GeV]')
B.pl.ylabel('')

B.pl.figure()
W_should_h.plot(label='W_should')
W_did_h.plot(label='W_did')
B.pl.title('H(e,e\'p) Singles Invariant Mass')
B.pl.xlabel('W [GeV]')
B.pl.ylabel('')
B.pl.text(0.99,100,'Proton Absorption:\n   $\\frac{e_{did}}{e_{should}}$'+\
          f' = {e_pAbs:0.4f}',fontsize=14)
B.pl.legend()

#%% get SHMS transmittance across ang acceptance

# xptar_did/xptar_should histograms
xptar_did_h = B.histo(pxptar_s[e_did],range=(-0.05,0.05),bins=15)
xptar_should_h = B.histo(pxptar_s[e_should],range=(-0.05,0.05),bins=15)

# yptar_did/yptar_should histograms
yptar_did_h = B.histo(pyptar_s[e_did],range=(-0.04,0.04),bins=15)
yptar_should_h = B.histo(pyptar_s[e_should],range=(-0.04,0.04),bins=15)

# get rid of zeroes in denominatyor
def set_zero(inpt,set_zero_to=1.):
    is_zero = np.where(inpt == 0.)
    for i in is_zero[0]:
        inpt[i] = set_zero_to
    print(f'Set {is_zero[0].size} inf to {set_zero_to}\n')
    return

set_zero(xptar_should_h.bin_content)
set_zero(yptar_should_h.bin_content)

# xptar, yptar ratios
xptar_ratio = xptar_did_h/xptar_should_h
yptar_ratio = yptar_did_h/yptar_should_h

new_pxptar_cut_lims = (-0.017,0.017)
new_pyptar_cut_lims = (-0.007,0.007)

# plot xptar_did/xptar_should
B.pl.figure()
xptar_should_h.plot(label='$X\'_{tar}$ Should')
xptar_did_h.plot(label='$X\'_{tar}$ Did')
B.pl.title('')
B.pl.xlabel('SHMS $X\'_{tar}$ [rad]')
B.pl.ylabel('')
B.pl.legend()

B.pl.figure()
xptar_ratio.plot_exp(ignore_zeros=True)
B.pl.hlines(0.9,-0.05,0.05,linestyles='--',colors='r')
B.pl.text(-0.03,0.92,'90%',color='r')
B.pl.ylim(0,1.2)
B.pl.vlines([new_pxptar_cut_lims[0],new_pxptar_cut_lims[1]],
            ymin=0,ymax=1.2,linestyles='--',colors='b')
B.pl.title('$X\'_{tar,did}/X\'_{tar,should}$')
B.pl.xlabel('SHMS $X\'_{tar}$ [rad]')
B.pl.ylabel('Proton Absorption Factor Ratio')

# plot yptar_did/yptar_should
B.pl.figure()
yptar_should_h.plot(label='$Y\'_{tar}$ Should')
yptar_did_h.plot(label='$Y\'_{tar}$ Did')
B.pl.title('')
B.pl.xlabel('SHMS $Y\'_{tar}$ [rad]')
B.pl.ylabel('')
B.pl.legend()

B.pl.figure()
yptar_ratio.plot_exp(ignore_zeros=True)
B.pl.hlines(0.9,-0.05,0.05,linestyles='--',colors='r')
B.pl.text(-0.03,0.92,'90%',color='r')
B.pl.ylim(0,1.2)
B.pl.vlines([new_pyptar_cut_lims[0],new_pyptar_cut_lims[1]],
            ymin=0,ymax=1.2,linestyles='--',colors='b')
B.pl.title('$Y\'_{tar,did}/Y\'_{tar,should}$')
B.pl.xlabel('SHMS $Y\'_{tar}$ [rad]')
B.pl.ylabel('Proton Absorption Factor Ratio')

#%% update cuts with new pxptar/pyptar cuts and plot invariant mass again

new_pxptar_cut = C.WCUT(new_pxptar_cut_lims[0],new_pxptar_cut_lims[1],
                        name='new_pxptar_cut')
new_pyptar_cut = C.WCUT(new_pyptar_cut_lims[0],new_pyptar_cut_lims[1],
                        name='new_pyptar_cut')

c_pxptar = new_pxptar_cut(pxptar_s)
c_pyptar = new_pyptar_cut(pyptar_s)

# e_should: number of electrons detected in SHMS with elastics kinematics and 
#   acceptance that SHOULD have a proton coincident on HMS side 
e_should = c_pdp & c_pxptar & c_pyptar & c_pztar & c_pcal & c_W

e_should_noW = c_pdp & c_pxptar & c_pyptar & c_pztar & c_pcal

# e_did: number of events that were actually detected on the HMS with a proton
#   in coincidence withing the HMS elastic acceptance region
e_did = e_should & c_hdp & c_hTRIG1

e_did_noW = e_should_noW & c_hdp & c_hTRIG1

# invariant mass histogram: e_did
W_did_h = B.histo(W_s[e_did],range=(0.8,1.2),bins=30)
W_did, W_did_e = W_did_h.sum()

W_did_nocut_h = B.histo(W_s[e_did_noW],range=(0.8,1.2),bins=30)

# invariant mass histogram: e_should
W_should_h = B.histo(W_s[e_should],range=(0.8,1.2),bins=30)
W_should, W_should_e = W_should_h.sum()

W_should_nocut_h = B.histo(W_s[e_should_noW],range=(0.8,1.2),bins=30)

# proton absorption factor
e_pAbs = (W_did/hTrkEff_coin)/W_should
e_pAbs_err = np.sqrt(W_should-(W_did/hTrkEff_coin))/W_should

print('Proton Absorption Factor: ',e_pAbs,'$\pm$',e_pAbs_err,
      ' \nor ',(1-e_pAbs)*100,'$\pm$',e_pAbs_err*100,'%')

# plot invariant mass
B.pl.figure()
W_should_nocut_h.plot(label='W_should',color=cmp['Paired'].colors[4])
W_did_nocut_h.plot(label='W_did',color=cmp['Paired'].colors[0])
B.pl.title('H(e,e\'p) Singles Invariant Mass')
B.pl.xlabel('W [GeV]')
B.pl.ylabel('')
B.pl.text(1.02,80,'Proton Transmission:\n$\\frac{e_{did}}{e_{should}}$'+\
          f' = {e_pAbs:0.3f}$\pm${e_pAbs_err:0.3f}',fontsize=14)
B.pl.legend()
B.pl.vlines([0.9,1.0],0,100,linestyles='--',colors='r')