#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:38:06 2026

@author: gvill
"""
import subprocess as sub

pm_set = input('Enter pm setting:\n')
model = input('Enter model name:\n')
data_set = int(input('Enter data sets to process:\n'))
if data_set == 1:
    data_set = ['set_1']
elif data_set == 2:
    data_set = ['set_1','set_2']

model_types = ['fsi_norad','pwia_norad']

# 1. calc_avg_kin.py
print('** Running average_kin/calc_avg_kin.py')
for mt in model_types:
    res1 = sub.run(['python3','yield_n_xsec/average_kin/calc_avg_kin.py'],
                   input=f'{pm_set} {model}{mt}',
                   text=True)
  
print('** Running average_xsec/calc_Xsec.py')
# 2. calc_Xsec.py
for s in data_set:
    res2 = sub.run(['python3','yield_n_xsec/average_xsec/calc_Xsec.py'],
                   input=f'{pm_set} {model} {s}',
                   text=True)

print('** Running jml_theory_xsec/calc_theory_Xsec.py')
# 3. calc_theory_Xsec.py
res3 = sub.run(['python3','yield_n_xsec/jml_theory_xsec/calc_theory_Xsec.py'],
               input=f'{pm_set}',
               text=True)

print('** Running bin_centering/calc_bc_corr.py')
# 4. calc_bc_corr.py
for s in data_set:
    res4 = sub.run(['python3','yield_n_xsec/bin_centering/calc_bc_corr.py'],
                   input=f'{pm_set} {model} {s}',
                   text=True)

print('** Running bin_centering/check_bc_corr.py')
# 5. check_bc_corr.py
for s in data_set:
    res5 = sub.run(['python3','yield_n_xsec/bin_centering/check_bc_corr.py'],
                   input=f'{pm_set} {model} {s}',
                   text=True)  
