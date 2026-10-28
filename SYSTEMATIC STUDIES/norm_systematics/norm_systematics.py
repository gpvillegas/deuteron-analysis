#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 12:04:24 2026

@author: gvill

NORMALIZATON SYSTEMATICS: This code calculates the systematic uncertainty due
to each normalization factor on the cross section. 3 errors are calculated: 
a run-by-run error added in quadrature at the end, an error on the regular 
average, and an error on the weighted average of a set per setting.

The normalization factors are: 
    1. f_tgtBoil --> Target Boiling Correction
    2. f_pTr --> Proton Transmission Correction
    3. f_pTrk --> SHMS Tracking Efficiency
    4. f_hTrk --> HMS Tracking Efficiency
    5. f_tLT --> EDTM Total Live Time
    6. f_Qtot --> Total accumulated charge
    
The total systematic uncertainty on the cross section due to these factors 
would be:
    dsig^2 = \sum[i=6 (dsig/df_i)^2 * df_i^2 ]
where sig is the experimental cross section (normalized) 
    and dsig/df_i = sig/f_i

So we end up with a simpler form:
    dsig^2 = sig^2 * \sum[i=6 (df_i/f_i)^2 ]

Explicitly this is:
    dsig^2 = sig^2 *( [df/f]_tgtBoil^2 + [df/f]_pTr^2 + [df/f]_pTrk^2 + \
         [df/f]_hTrk^2 + [df/f]_tLT^2 + [df/f]_Qtot^2)

So the systematic error squared is the sum of the squares of the relative errors 
of all efficiency and correction factors times the cross section squared.

The uncertainties of each factor are located:
    1. tgtBoil: written in DB per run
    2. pTr: one overall = 0.953 +/- 0.011
    3. pTrk: written in DB per run
    4. hTrk: written in DB per run
    5. tLT: written in DB per run
    6. Qtot: TBD, one overall, estimate dI/I ~ 2%, charge uncertainty has been
        determined based on BCM calibration
"""

import numpy as np
import database_operations as db
import LT.box as B
from matplotlib import colormaps as cmp

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 18
}
#%%

# helper functions
def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr
      
def get_list(db_res, index = 0):
    # To get rid of unhelpful groupings (lists with only 1 item, etc) 
    l = [k[index] for k in db_res]
    if len(l) > 1:
        return l
    else:
        return l[0]
    
def get_run_nums(pm_set,data_set):
    r = get_list(db.retrieve('deuteron_db.db','run','RUN_LIST_UPDATED', 
                          where = f"setting=\'{pm_set}\'"+\
                              f"AND data_set=\'{data_set}\'"))
    return r
    
def load_xsec_file(pm_set,data_set):
    filename = f'yield_n_xsec/bin_centering/acos_phi/{pm_set}_jml_bc_corr_{data_set}.txt'
    f = B.get_file(filename)
    xsec = convert2NaN(f['fsiRC_dataXsec_fsibc_corr'], value=-1)
    xsec = convert2NaN(xsec,value=np.inf)
    dxsec = convert2NaN(f['fsiRC_dataXsec_fsibc_corr_err'], value=-1)
    dxsec = convert2NaN(dxsec,value=np.inf)
    return (f, xsec, dxsec)
    
class norm_syst:
    def __init__(self,Xsec,run_dict={},name=''):
        self.runs = run_dict
        self.Xsec = Xsec
        self.name = name
        self.error_name = name+'_Err'
        
        self.__get_syst_variable() # method defines self.syst_var and 
                                    # self.syst_var_err
        self.__get_relative_error() # method to calculate relative error
                                        # we use f = 1/syst_var
                                        # df = -(1/syst_var**2)*syst_var_err
                                        # relative error = df/f
                    
    def calc_syst(self,run_by_run=False,avg=False,wavg=False):
        if run_by_run:
            self.__get_syst_error_rbr()
        elif avg:
            self.__get_syst_error_avg()
        elif wavg:
            self.__get_syst_error_wavg()
           
    def __get_syst_error_rbr(self):
        if self.name == 'pTr':
            self.syst_err = {}
            for s in self.Xsec: # loop over settings
                se = {}
                for d in self.Xsec[s]: # loop over data set
                    se[d] = self.Xsec[s][d]**2 * self.relative_error**2
                self.syst_err[s] = se  
        else:
            self.syst_err = {}
            self.syst_err_run = {}
            for s in self.runs: # loop over settings
                se = {}
                se_run = {}
                for d in self.runs[s]: # loop over data set
                    e = []
                    for i in range(len(self.runs[s][d])): # loop over runs 
                                                            # (1 rel_err_factor per run)
                        e.append(self.Xsec[s][d]**2 * self.relative_error[s][d][i]**2)
                    se_run[d] = np.array(e)
                    se[d] = np.linalg.norm(e, axis=0)
                self.syst_err_run[s] = se_run
                self.syst_err[s] = se  
    
    def __get_syst_error_avg(self):
        if self.name == 'pTr':
            self.syst_err = {}
            for s in self.Xsec: # loop over settings
                se = {}
                for d in self.Xsec[s]: # loop over data set
                    se[d] = self.Xsec[s][d]**2 * self.relative_error**2
                self.syst_err[s] = se  
        else:
            self.relative_error_avg = {}
            self.syst_err = {}
            for s in self.runs: # loop over settings
                se = {}
                avg = {}
                for d in self.runs[s]: # loop over data set
                    a = np.average(self.relative_error[s][d]) # average over data_set
                    e = self.Xsec[s][d]**2 * a**2
                    
                    se[d] = np.array(e)
                    avg[d] = a
                self.syst_err[s] = se
                self.relative_error_avg[s] = avg
                
    def __get_syst_error_wavg(self):
        if self.name == 'pTr':
            self.syst_err = {}
            for s in self.Xsec: # loop over settings
                se = {}
                for d in self.Xsec[s]: # loop over data set
                    se[d] = self.Xsec[s][d]**2 * self.relative_error**2
                self.syst_err[s] = se  
        else:
            self.relative_error_wavg = {}
            self.weights = {}
            self.syst_err = {}
            for s in self.runs: # loop over settings
                se = {}
                avg = {}
                wei = {}
                for d in self.runs[s]: # loop over data set
                    # weights
                    w = 1/self.df[s][d]**2
                    w[w == np.inf] = 0.0
                    if any(w != 0.0):
                        # get weighted average of f = 1/syst_var
                        f = np.average(self.f[s][d],weights=w) # average over data_set
                        df = np.sqrt(1./np.sum(w))
                        
                        a = df/f 
                        e = self.Xsec[s][d]**2 * a**2
                    else:
                        a = 0.0 
                        e = self.Xsec[s][d]**2 * a**2
                        
                    se[d] = np.array(e)
                    avg[d] = a
                    wei[d] = w
                    
                self.syst_err[s] = se
                self.relative_error_wavg[s] = avg
                self.weights[s] = wei
                
    def __get_syst_variable(self):
        if self.name == 'pTr':
            self.syst_var = 0.9530
            self.syst_var_err = 0.0110
        # elif self.name == 'Qtot':
        #     self.syst_var = 0.0
        #     self.syst_var_err = 0.0
        else:
            self.syst_var = {}
            self.syst_var_err = {}
            for s in self.runs: # loop over settings
                a1 = {}
                a2 = {}
                for d in self.runs[s]: # loop over data set
                    a11 = []
                    a22 = []
                    for r in self.runs[s][d]: # loop over run number
                        t, t_err =\
                            db.retrieve('deuteron_db.db',
                                        f'{self.name}, {self.error_name}',
                                        'RUN_LIST_UPDATED', where = f"run=\'{r}\'")[0]
                        a11.append(t)   # systematic variable
                        a22.append(t_err) # systematic variable error
                    a1[d] = np.array(a11)
                    a2[d] = np.array(a22)
                self.syst_var[s] = a1
                self.syst_var_err[s] = a2
                
    def __get_relative_error(self):
        if self.name == 'pTr':
            self.f = 1/self.syst_var
            self.df = -(1/self.syst_var**2)*self.syst_var_err
            self.relative_error = np.abs(self.df/self.f)
        # elif self.name == 'Qtot':
        #     self.relative_error = 0.02 # hard coded until figured out... 
        else:
            self.f = {}
            self.df = {}
            self.relative_error = {}
            for s in self.runs:
                f = {}
                df = {}
                rel_err = {}
                for d in self.runs[s]:
                    f[d] = 1/self.syst_var[s][d]
                    df[d] = -(1/self.syst_var[s][d]**2)*self.syst_var_err[s][d]
                    rel_err[d] = np.abs(df[d]/f[d])
                self.f[s] = f
                self.df[s] = df
                self.relative_error[s] = rel_err
            

# get deep run numbers from runlist 
r_120_1 = get_run_nums('pm_120','set_1')
r_120_2 = get_run_nums('pm_120','set_2')
r_580_1 = get_run_nums('pm_580','set_1')
r_580_2 = get_run_nums('pm_580','set_2')
r_800_1 = get_run_nums('pm_800','set_1')
r_900_1 = get_run_nums('pm_900','set_1')

deep_runs = {'pm_120':{'set_1':r_120_1,'set_2':r_120_2},
             'pm_580':{'set_1':r_580_1,'set_2':r_580_2},
             'pm_800':{'set_1':r_800_1},
             'pm_900':{'set_1':r_900_1}}

# get xsec files
fPm120_1, dataXsec_120_1, dataXsec_err_120_1 = load_xsec_file('pm_120','set_1')
fPm120_2, dataXsec_120_2, dataXsec_err_120_2 = load_xsec_file('pm_120','set_2')
fPm580_1, dataXsec_580_1, dataXsec_err_580_1 = load_xsec_file('pm_580','set_1')
fPm580_2, dataXsec_580_2, dataXsec_err_580_2 = load_xsec_file('pm_580','set_2')
fPm800_1, dataXsec_800_1, dataXsec_err_800_1 = load_xsec_file('pm_800','set_1')
fPm900_1, dataXsec_900_1, dataXsec_err_900_1 = load_xsec_file('pm_900','set_1')

dataXsec = {'pm_120': {'set_1':dataXsec_120_1,'set_2':dataXsec_120_2},
            'pm_580': {'set_1':dataXsec_580_1,'set_2':dataXsec_580_2},
            'pm_800': {'set_1':dataXsec_800_1},
            'pm_900': {'set_1':dataXsec_900_1},}

dataXsec_err = {'pm_120': {'set_1':dataXsec_err_120_1,'set_2':dataXsec_err_120_2},
                'pm_580': {'set_1':dataXsec_err_580_1,'set_2':dataXsec_err_580_2},
                'pm_800': {'set_1':dataXsec_err_800_1},
                'pm_900': {'set_1':dataXsec_err_900_1},}

# get bin information (for plotting)
ib = fPm120_1['i_b'] # bin number
xb = fPm120_1['xb'] # thnq bin center
yb = fPm120_1['yb'] # pm bin center

# Target Boiling
tgtBoil_syst = norm_syst(dataXsec,deep_runs,name='tgtBoil_corr')  
tgtBoil_syst.calc_syst(run_by_run=True)

tgtBoil_syst_avg = norm_syst(dataXsec,deep_runs,name='tgtBoil_corr') 
tgtBoil_syst_avg.calc_syst(avg=True)

tgtBoil_syst_wavg = norm_syst(dataXsec,deep_runs,name='tgtBoil_corr') 
tgtBoil_syst_wavg.calc_syst(wavg=True)

# Proton Absorption
pTr_syst = norm_syst(dataXsec,deep_runs,name='pTr')
pTr_syst.calc_syst(run_by_run=True)

# SHMS Tracking Efficiency
pTrk_syst = norm_syst(dataXsec,deep_runs,name='SHMS_TrkEff')
pTrk_syst.calc_syst(run_by_run=True)

pTrk_syst_avg = norm_syst(dataXsec,deep_runs,name='SHMS_TrkEff')
pTrk_syst_avg.calc_syst(avg=True)

pTrk_syst_wavg = norm_syst(dataXsec,deep_runs,name='SHMS_TrkEff')
pTrk_syst_wavg.calc_syst(wavg=True)

# HMS Tracking Efficiency
hTrk_syst = norm_syst(dataXsec,deep_runs,name='HMS_TrkEff')
hTrk_syst.calc_syst(run_by_run=True)

hTrk_syst_avg = norm_syst(dataXsec,deep_runs,name='HMS_TrkEff')
hTrk_syst_avg.calc_syst(avg=True)

hTrk_syst_wavg = norm_syst(dataXsec,deep_runs,name='HMS_TrkEff')
hTrk_syst_wavg.calc_syst(wavg=True)

# Total Livetime
tLT_syst = norm_syst(dataXsec,deep_runs,name='tLT')
tLT_syst.calc_syst(run_by_run=True)

tLT_syst_avg = norm_syst(dataXsec,deep_runs,name='tLT')
tLT_syst_avg.calc_syst(avg=True)

tLT_syst_wavg = norm_syst(dataXsec,deep_runs,name='tLT')
tLT_syst_wavg.calc_syst(wavg=True)

# Total Charge
Qtot_syst = norm_syst(dataXsec,deep_runs,name='BCM1_charge')
Qtot_syst.calc_syst(run_by_run=True)

Qtot_syst_avg = norm_syst(dataXsec,deep_runs,name='BCM1_charge')
Qtot_syst_avg.calc_syst(avg=True)

Qtot_syst_wavg = norm_syst(dataXsec,deep_runs,name='BCM1_charge')
Qtot_syst_wavg.calc_syst(wavg=True)

#%%
# write results to a file and plot
dir_name = 'SYSTEMATIC STUDIES/norm_systematics/'

header = """# NORMALIZATION SYSTEMATICS SUMMARY
#
# The normalization factors are: 
#    1. f_tgtBoil --> Target Boiling Correction
#    2. f_pTr --> Proton Transmission Correction
#    3. f_pTrk --> SHMS Tracking Efficiency
#    4. f_hTrk --> HMS Tracking Efficiency
#    5. f_tLT --> EDTM Total Live Time
#    6. f_Qtot --> Total accumulated charge
#
# header definitions
#   dsig2_syst_[x]: where x is any norm factor type, this is the square 
#                   of the syst error contribution, meaning 
#                   Xsec**2 * relative_error**2
#   dsig_norm_syst: total normalization systematic error, adding up all 
#                   contributions and taking the square root
#   rel_err_tot: total realtive error, adding up in quadrature all 
#               contributions relative errors.
#               relative errors are in %
# header line:
#! dsig2_syst_tgtBoil[f,0]/ dsig2_syst_pTr[f,1]/ dsig2_syst_pTrk[f,2]/ dsig2_syst_hTrk[f,3]/ dsig2_syst_tLT[f,4]/ dsig2_syst_Qtot[f,5]/ dsig_norm_syst[f,6]/ rel_err_tgtBoil[f,7]/ rel_err_pTr[f,8]/ rel_err_pTrk[f,9]/ rel_err_hTrk[f,10]/ rel_err_tLT[f,11]/ rel_err_Qtot[f,12]/ rel_err_total[f,13]/
"""

for s in deep_runs:
    for d in deep_runs[s]:
        xsec = dataXsec[s][d]
        # combine systematics to get total contribution
        a1 = tgtBoil_syst_wavg.syst_err[s][d]
        b1 = pTr_syst.syst_err[s][d]
        c1 = pTrk_syst_wavg.syst_err[s][d]
        d1 = hTrk_syst_wavg.syst_err[s][d]
        e1 = tLT_syst_wavg.syst_err[s][d]
        f1 = Qtot_syst_wavg.syst_err[s][d]
        
        norm_syst_total = np.sqrt(a1+b1+c1+d1+e1+f1)
        
        # calculate relative errors in %
        a2 = 100*np.sqrt(a1)/xsec
        b2 = 100*np.sqrt(b1)/xsec
        c2 = 100*np.sqrt(c1)/xsec
        d2 = 100*np.sqrt(d1)/xsec
        e2 = 100*np.sqrt(e1)/xsec
        f2 = 100*np.sqrt(f1)/xsec

        rel_error_total = 100*norm_syst_total/xsec        

        fname = f'{s}_norm_systematics_{d}.txt'
        f = open(dir_name+fname,'w')
        
        f.write(header)
        for i in range(len(ib)): # loop over bins
            # print(i)
            # print(f'{a1[i]:.4e} {b1[i]:.4e} {c1[i]:.4e} {d1[i]:.4e} {e1[i]:.4e} {f1[i]:.4e} {norm_syst_total[i]:.4e} {a2[i]:.4e} {b2[i]:.4e} {c2[i]:.4e} {d2[i]:.4e} {e2[i]:.4e} {f2[i]:.4e} {rel_error_total[i]:.4e}\n')
            f.write(f'{a1[i]:.4e} {b1[i]:.4e} {c1[i]:.4e} {d1[i]:.4e} {e1[i]:.4e} {f1[i]:.4e} {norm_syst_total[i]:.4e} {a2[i]:.4f} {b2[i]:.4f} {c2[i]:.4f} {d2[i]:.4f} {e2[i]:.4f} {f2[i]:.4f} {rel_error_total[i]:.4f}\n')
        f.close()    
        
#%%       
# for plotting
thnq = np.unique(xb)
pm = np.unique(yb)
y = np.zeros(pm.size)

# plot systematic errors
for s in deep_runs:
    for d in deep_runs[s]:
        xsec = dataXsec[s][d]
        dxsec = dataXsec_err[s][d]
        
        a1 = tgtBoil_syst_wavg.syst_err[s][d]
        b1 = pTr_syst.syst_err[s][d]
        c1 = pTrk_syst_wavg.syst_err[s][d]
        d1 = hTrk_syst_wavg.syst_err[s][d]
        e1 = tLT_syst_wavg.syst_err[s][d]
        f1 = Qtot_syst_wavg.syst_err[s][d]
        
        norm_syst_total = np.sqrt(a1+b1+c1+d1+e1+f1)
        
        # calculate relative errors in %
        a2 = 100*np.sqrt(a1)/xsec
        b2 = 100*np.sqrt(b1)/xsec
        c2 = 100*np.sqrt(c1)/xsec
        d2 = 100*np.sqrt(d1)/xsec
        e2 = 100*np.sqrt(e1)/xsec
        f2 = 100*np.sqrt(f1)/xsec

        rel_error_total = 100*norm_syst_total/xsec
        
        # calculate statistical rel error
        rel_error_stats = 100*dxsec/xsec
        
        for t in [35.]:
            B.pl.figure(figsize=(8,6))
            
            # B.plot_exp(pm,y,rel_error_stats[xb==t],label='Stat. Error',color='gray')
            
            B.plot_exp(pm,y,rel_error_total[xb==t],label='Total',color='black')
            B.plot_exp(pm,y,a2[xb==t],label='Target Boiling',color=cmp['tab20'](6))
            B.plot_exp(pm,y,b2[xb==t],label='Proton Absorption',color=cmp['tab20'](2))
            B.plot_exp(pm,y,c2[xb==t],label='SHMS Tracking Eff',color=cmp['tab20'](16))
            B.plot_exp(pm,y,d2[xb==t],label='HMS Tracking Eff',color=cmp['tab20'](4))
            B.plot_exp(pm,y,e2[xb==t],label='Livetime',color=cmp['tab20'](18))
            B.plot_exp(pm,y,f2[xb==t],label='Charge',color=cmp['tab20'](0))
            
            #'Normalization Systematic Uncertainties:'
            B.pl.title(f'\n{s} {d}'+' $\\theta_{nq}$'+f'={t}')
            B.pl.xlabel('$P_m$ [GeV]')
            B.pl.ylabel('Relative Error [%]')
            B.pl.legend(fontsize='small')
            # B.pl.yscale('log')
            B.pl.savefig(dir_name+f'plots_wavg_newq/{s}_{d}_rel_err_syst_thnq_{t}.pdf')
            
#%% calculate weighted averages of individual efficiencies/normalization values
header="""# Normalization Systematics Summary
#This file contains the weighted averages of each normalization factor 
#per mising momentum setting per run set.
#
#header line:
#! Pm[s,0]/ set[s,1]/ tgtBoil[f,2]/ pTr[f,3]/ pTrk[f,4]/ hTrk[f,5]/ tLT[f,6]/ Qtot[f,7]/
"""

fname = 'norm_systematics_normfac_summary.txt'
f = open(dir_name+fname,'w')
f.write(header)

for s in deep_runs:
    for d in deep_runs[s]: 
        #Target Boiling
        w = 1/tgtBoil_syst_wavg.syst_var_err[s][d]**2
        w[w==np.inf]=0.0
        v = tgtBoil_syst_wavg.syst_var[s][d]
        wavg_tgtBoil = np.average(v,weights=w)
        
        #Proton Transmission
        pTr = pTr_syst.syst_var
        
        #SHMS Tracking efficiency
        w = 1/pTrk_syst_wavg.syst_var_err[s][d]**2
        w[w==np.inf]=0.0
        v = pTrk_syst_wavg.syst_var[s][d]
        wavg_pTrk = np.average(v,weights=w)
        
        #HMS Tracking efficiency
        w = 1/hTrk_syst_wavg.syst_var_err[s][d]**2
        w[w==np.inf]=0.0
        v = hTrk_syst_wavg.syst_var[s][d]
        wavg_hTrk = np.average(v,weights=w)
        
        #Total EDTM Livetime
        w = 1/tLT_syst_wavg.syst_var_err[s][d]**2
        w[w==np.inf]=0.0
        v = tLT_syst_wavg.syst_var[s][d]
        wavg_tLT = np.average(v,weights=w)
        
        #Total Charge
        v = Qtot_syst_wavg.syst_var[s][d]
        Qtot = np.sum(v)
        
        f.write(f'{s} {d} {wavg_tgtBoil:.4f} {pTr:.4f} {wavg_pTrk:.4f} {wavg_hTrk:.4f} {wavg_tLT:.4f} {Qtot:.4f}\n')

f.close() 
        