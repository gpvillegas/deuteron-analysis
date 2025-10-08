#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:36:33 2025

@author: gvill

Plot reference times
"""

import data_init as D
import database_operations as db
import LT.box as B

import re
import numpy as np
import matplotlib.pyplot as plt

# %% helper functions


def pptxify(t='', x='', y='', fsize=24, fig_width=10, fig_length=10):
    B.pl.tick_params('both', labelsize='x-large')
    B.pl.title(t, fontdict={'fontsize': fsize})
    B.pl.xlabel(x, fontdict={'fontsize': fsize})
    B.pl.ylabel(y, fontdict={'fontsize': fsize})
    fig = B.pl.gcf()
    fig.set_size_inches(fig_width, fig_length)
    # B.pl.legend(fontsize = fsize)


# %% constants
chan_per_ns = 0.09766

# %%
deut_db = 'deuteron_db.db'

hDCREF_branches = D.get_list(db.retrieve(deut_db, 'Branches', 'TTree_Branches',
                                         where='Branches like \'T.coin.hDCREF__tdc%\''))
pDCREF_branches = D.get_list(db.retrieve(deut_db, 'Branches', 'TTree_Branches',
                                         where='Branches like \'T.coin.pDCREF__tdc%\''))

hDCREF_tdcTimeRaw_branches = D.get_list(db.retrieve(deut_db, 'Branches', 'TTree_Branches',
                                                    where='Branches like \'T.coin.hDCREF__tdcTimeRaw%\''))
pDCREF_tdcTimeRaw_branches = D.get_list(db.retrieve(deut_db, 'Branches', 'TTree_Branches',
                                                    where='Branches like \'T.coin.pDCREF__tdcTimeRaw%\''))

br_sel = hDCREF_branches + pDCREF_branches +\
    ['T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw',
     'T.coin.hFADC_TREF_ROC1_adcMultiplicity',
     'T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw',
     'T.coin.pFADC_TREF_ROC2_adcMultiplicity',
     'T.coin.hT1_tdcTimeRaw', 'T.coin.hT1_tdcMultiplicity',
     'T.coin.hT2_tdcTimeRaw', 'T.coin.hT2_tdcMultiplicity',
     'T.coin.pT1_tdcTimeRaw', 'T.coin.pT1_tdcMultiplicity',
     'T.coin.pT2_tdcTimeRaw', 'T.coin.pT2_tdcMultiplicity']

# %% load root files
# heep = D.DATA_INIT(data_type='deut23_data', kin_study='heep_coin',
#                    select_branches=br_sel)

run = D.DATA_INIT(data_type='deut23_data', run=20851,
                  select_branches=br_sel)

# %%
DCREF_tdcTimeRaw_pat = re.compile(".DCREF._tdcTimeRaw")
DCREF_tdcTimeRaw_list = [m for m in br_sel if DCREF_tdcTimeRaw_pat.search(m)]

TimeRaw_histos = {}
for v in DCREF_tdcTimeRaw_list:
    t = run.Branches[v]
    TimeRaw_histos[v] = B.histo(t, range=(0, 19000), bins=100)

HODO_tdcTimeRaw_pat = re.compile(".T._tdcTimeRaw")
HODO_tdcTimeRaw_list = [m for m in br_sel if HODO_tdcTimeRaw_pat.search(m)]

for v in HODO_tdcTimeRaw_list:
    t = run.Branches[v]
    TimeRaw_histos[v] = B.histo(t, range=(0, 6000), bins=100)

adcTimeRaw_pat = re.compile("._adcPulseTimeRaw")
adcTimeRaw_list = [m for m in br_sel if adcTimeRaw_pat.search(m)]

for v in adcTimeRaw_list:
    t = run.Branches[v]
    TimeRaw_histos[v] = B.histo(t, range=(0, 7000), bins=100)


# %% create multiplicity cut == 1
pattern = re.compile("._adcMultiplicity")

adcMult_list = [m for m in br_sel if pattern.search(m)]

pattern = re.compile("._tdcMultiplicity")

tdcMult_list = [m for m in br_sel if pattern.search(m)]

Mult_list = adcMult_list + tdcMult_list

Mult_cut_dict = {}
for cut in Mult_list:
    p = re.compile(cut.strip('Multiplicity') + '.')
    match = [m for m in TimeRaw_histos.keys() if p.search(m)]
    Mult_cut_dict[match[0]] = {cut: 1}

mult_cut = {}
for k in Mult_cut_dict:
    mult_cut[k] = run.CUT(cut_dict=Mult_cut_dict[k],
                          where_branches=run.Branches)

# %%

TimeRaw_multcut_histos = {}
for v in DCREF_tdcTimeRaw_list:
    t = run.Branches[v]
    t_cut = t[mult_cut[v]]
    TimeRaw_multcut_histos[v] = B.histo(t_cut, range=(0, 19000), bins=100)

for v in HODO_tdcTimeRaw_list:
    t = run.Branches[v]
    t_cut = t[mult_cut[v]]
    TimeRaw_multcut_histos[v] = B.histo(t_cut, range=(0, 6000), bins=100)

for v in adcTimeRaw_list:
    t = run.Branches[v]
    t_cut = t[mult_cut[v]]
    TimeRaw_multcut_histos[v] = B.histo(t_cut, range=(0, 7000), bins=100)

# %% plot ref time plots

for v in TimeRaw_histos:
    B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    B.pl.legend()
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()

# %%
hdc_tdcrefcut = 16000.
hhodo_tdcrefcut = 3500.
hhodo_adcrefcut = 5000.
hcer_adcrefcut = 5000.
hcal_adcrefcut = 5000.

pdc_tdcrefcut = 15500.
phodo_tdcrefcut = 4700.
phodo_adcrefcut = 5600.
pngcer_adcrefcut = 5600.
phgcer_adcrefcut = 5600.
paero_adcrefcut = 5600.
pcal_adcrefcut = 5600.

# %% hDC
for v in hDCREF_tdcTimeRaw_branches:
    # B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    (y1, y2) = plt.ylim()

    B.pl.vlines(hdc_tdcrefcut, y1, y2, colors='black', linestyles='--')
    B.pl.legend()
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()
# %%
for v in pDCREF_tdcTimeRaw_branches:
    # B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    (y1, y2) = plt.ylim()

    B.pl.vlines(pdc_tdcrefcut, y1, y2, colors='black', linestyles='--')
    B.pl.legend()
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()

# %%
for v in ['T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw']:
    # B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    (y1, y2) = plt.ylim()

    B.pl.vlines(hhodo_adcrefcut, y1, y2, colors='black', linestyles='--')
    B.pl.legend()
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()

# %%
for v in ['T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw']:
    # B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    (y1, y2) = plt.ylim()

    B.pl.vlines(phodo_adcrefcut, y1, y2, colors='black', linestyles='--')
    B.pl.legend()
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()
# %%
for v in ['T.coin.hT2_tdcTimeRaw']:
    # B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    (y1, y2) = plt.ylim()

    B.pl.vlines(hhodo_tdcrefcut, y1, y2, colors='black', linestyles='--')
    B.pl.legend()
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()

# %%

for v in ['T.coin.pT1_tdcTimeRaw', 'T.coin.pT2_tdcTimeRaw']:
    # B.pl.figure()
    h = TimeRaw_histos[v]
    h_multcut = TimeRaw_multcut_histos[v]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')

    (y1, y2) = plt.ylim()

    B.pl.vlines(phodo_tdcrefcut, y1, y2, colors='black', linestyles='--')
    B.pl.legend(fontsize=14, loc='upper left')
    B.pl.yscale('log')
    pptxify(f'Reference Time: {v}',
            '[channel]',
            fsize=14, fig_width=7, fig_length=5)

    B.pl.show()

# %%
plot_same = ['T.coin.pDCREF1_tdcTimeRaw', 'T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw',
             'T.coin.pT1_tdcTimeRaw', 'T.coin.hDCREF1_tdcTimeRaw',
             'T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw',
             'T.coin.hT2_tdcTimeRaw']
vlines = [pdc_tdcrefcut, phodo_adcrefcut, phodo_tdcrefcut, hdc_tdcrefcut,
          hhodo_adcrefcut, hhodo_tdcrefcut]

fig = B.pl.figure(layout='constrained')
# fig.suptitle(f'SHMS Drift time vs. Wire number\n Run {r}', fontsize=16)
rows = 2
cols = 3
i = 1

for pl in plot_same:
    ax = plt.subplot(rows, cols, i)
    h = TimeRaw_histos[pl]
    h_multcut = TimeRaw_multcut_histos[pl]

    h.plot(hatch='......', facecolor='white', edgecolor='#7b8fd4', label='raw')
    # h_multcut.plot(color='#de425b',alpha=0.8)
    h_multcut.plot(hatch='......', facecolor='white', edgecolor='#de425b',
                   label='m=1')
    (y1, y2) = plt.ylim()

    B.pl.vlines(vlines[i-1], y1, y2, colors='black', linestyles='--')
    B.pl.legend(fontsize=14, loc='upper left')
    B.pl.yscale('log')

    i += 1
    ax.set_xlabel('[channel]', fontdict={'fontsize': 14})
    ax.set_ylabel('')
    ax.set_title(f'{v}', fontdict={'fontsize': 14})

    fig = B.pl.gcf()
    fig.set_size_inches(18, 20)
