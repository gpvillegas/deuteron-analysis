#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 10:18:29 2026

@author: gvill
"""

import LT.box as B

fname = 'SYSTEMATIC STUDIES/kin_systematics/residuals_model2_betterfit_5runs.txt'
f = B.get_file(fname)

x = f['x']
res = f['residual']
res_err = f['residual_err']
dOm = f['dObs_meas']
dOm_err = f['dObs_meas_err']
dOmod = f['dObs_model']
dOmod_err = f['dObs_model_err']

B.pl.figure()
B.plot_exp(x,res,res_err,label='Residuals')
B.plot_exp(x,dOm,dOm_err,label='dObs_meas')
B.plot_exp(x,dOmod,dOmod_err,label='dObs_model')
xlo,xhi=B.pl.xlim()
B.pl.hlines(0,xlo,xhi,color='r')
B.pl.legend()
B.pl.title('$\chi^2$ Minimization for H(e,e\'p) Coin Delta = [-8,-4,0,+4,+8]')
B.pl.ylabel('Residuals')
B.pl.xlabel('dObs #: [dW, dEm, dPmx, dPmz]')

# B.pl.ylim(-0.015,0.015)