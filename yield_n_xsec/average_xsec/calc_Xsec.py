#Code to extract the average DATA and SIMC cross sections from the 2D Pm vs. theta_nq Histos

import os
import sys
from sys import argv
import getopt

from LT import datafile
import bin_info2 as BI

#Set proper paths to import ROOT
# sys.path.append('../../../../pyroot/')
# sys.path.append('/apps/root/PRO/lib/')
# sys.path.append('/apps/root/PRO/')

# do the root operations directly here
# import ROOT as R
# from ROOT import *
# from ROOT import *

import pickle as P
import numpy as np

# from math import *
# import everything for handling 2d histograms
# import bin_info as BI
# use a corrected version (this ignores the overflow bins)
# import bin_info1 as BI

# some constants
dtr = np.pi/180.

#MeV
MP = 938.272
MN = 939.566
#MD = 1875.6127
MD = 1875.61
me = 0.51099; 


#------------------------------------------------------------
# header information for the output file
header = \
"""
# averaged kinematics results
# averaged kinematic varibles used as input to calculate the averaged 
# kinematics: Ei, omega, th_e, pf
#
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/ pwiaXsec[f,5]/  pwiaXsec_err[f,6]/  fsiXsec[f,7]/   fsiXsec_err[f,8]/  pwiaRC_dataXsec[f,9]/  pwiaRC_dataXsec_err[f,10]/  fsiRC_dataXsec[f,11]/  fsiRC_dataXsec_err[f,12]/ 
"""
#------------------------------------------------------------


#CODE USAGE: /apps/python/2.7.12/bin/python.py calc_Xsec.py 580  1 syst_ext

#User Input
pm_set, model, data_set = input('Enter pm setting, model, and data set: \n').split()
# pm_set = int(sys.argv[1])
# data_set = int(sys.argv[2])
# sys_ext = sys.argv[3]

# dir_name = "./%s"%(sys_ext)
# #check if output directory exists, else creates it.
# if not os.path.exists(dir_name):
#    os.makedirs(dir_name)

# print argv

#create output file to write Xsec avg over kin. bin
# if pm_set == 80:
#    output_file = './%s/pm%i_laget.txt'%(sys_ext, pm_set)
# else:
output_file = './yield_n_xsec/average_xsec/%s_%s_average_xsec_%s.txt'%(pm_set, model, data_set)


o = open(output_file,'w')

#----------------------------------Open root file that sets histo binning---------------------------------------

# bin_file = '../../root_files/average_kinematics/%s/deep_simc_histos_pm%i_lagetpwia_norad_set%i.root'%(sys_ext, pm_set, data_set) 

# bf = R.TFile(bin_file)
# bf.cd()
filename = f'yield_n_xsec/average_kin/{pm_set}_{model}_pwia_norad_2Dhistos_avgKin.pcl'
with open(filename, 'rb') as f:
    h = P.load(f)
all = BI.get_histo_data_arrays(h['Pm_vs_thnq_v'],LT_histo=True)
# bf.Close()
#----------------------------------------------------------------------------------------------------------------

#Open root file to read cross sections from 2D Pm vs. thnq
# root_file_pwia = '../../root_files/pm%i_pwiaXsec_set%i_%s/Xsec_pm%i_lagetpwia_dataset%i.root'%(pm_set, data_set, sys_ext, pm_set, data_set)  
# root_file_fsi = '../../root_files/pm%i_fsiXsec_set%i_%s/Xsec_pm%i_lagetfsi_dataset%i.root'%(pm_set, data_set, sys_ext, pm_set, data_set)   

# # open PWIA file
# rf_pwia = R.TFile(root_file_pwia)
# # open FSI file
# rf_fsi = R.TFile(root_file_fsi)

# #Change to FSI file
# rf_fsi.cd()

# open pickle file with xsec (all data, simc)
filename = f'yield_n_xsec/{pm_set}_data_{model}_histos.pcl'
with open(filename, 'rb') as f:
    h = P.load(f)

# write the necessary header parameters
o.write('# histogram parameters \n')
o.write('#\ dx = {0:}\n'.format(repr(all.dx)))
o.write('#\ dy = {0:}\n'.format(repr(all.dy)))
o.write('#\ nx = {0:}\n'.format(repr(all.nx)))
o.write('#\ ny = {0:}\n'.format(repr(all.ny)))
o.write('#\ xmin = {0:}\n'.format(repr(all.xmin)))
o.write('#\ ymin = {0:}\n'.format(repr(all.ymin)))
# write header
o.write(header)


#Get 2D Histogram Bin Info (Pm, thnq) Cross Section 
#bin_info_dataXsec        = BI.get_histo_data_arrays(rf_fsi.H_data2DXsec) 
# FSI histos        
bin_info_fsiRC_dataXsec  = BI.get_histo_data_arrays(h[f'{pm_set}_{model}_FSI_data_xsec'][data_set],LT_histo=True)         #FSI radiative corrected data Xsec
bin_info_fsiXsec         = BI.get_histo_data_arrays(h[f'{pm_set}_{model}_fsi_xsec'],LT_histo=True)         

# rf_fsi.Close()

#change to PWIA file and get simc histo
# rf_pwia.cd()
# bin_info_pwiaRC_dataXsec  = BI.get_histo_data_arrays(rf_pwia.H_data2DXsec)         #PWIA radiative corrected data Xsec
# bin_info_pwiaXsec        = BI.get_histo_data_arrays(rf_pwia.H_simc2DXsec)
bin_info_pwiaRC_dataXsec  = BI.get_histo_data_arrays(h[f'{pm_set}_{model}_PWIA_data_xsec'][data_set],LT_histo=True)         #PWIA radiative corrected data Xsec
bin_info_pwiaXsec        = BI.get_histo_data_arrays(h[f'{pm_set}_{model}_pwia_xsec'],LT_histo=True)  

# rf_pwia.Close()


#Retrieve the 2D bin content for each of the files
cont1        = bin_info_fsiRC_dataXsec.cont
cont2        = bin_info_fsiXsec.cont
cont3        = bin_info_pwiaRC_dataXsec.cont
cont4        = bin_info_pwiaXsec.cont

ibin                   = bin_info_fsiRC_dataXsec.i  #Get 2D bin number (should be the same regardless of ROOTfile, given that the binning is the same)
i_xbin                 = bin_info_fsiRC_dataXsec.ix
i_ybin                 = bin_info_fsiRC_dataXsec.iy
thnq_b                 = bin_info_fsiRC_dataXsec.xb
pm_b                   = bin_info_fsiRC_dataXsec.yb

#Loop over bin number (xbin, ybin)->(th_nq_bin, Pm_bin)
for i,ib in enumerate(ibin):

   
   #Check Bin Content
   if(cont1[i] == 0.):
      fsiRC_dataXsec        = -1.        #data Xsec radiatively corrected using the Laget FSI model    
      fsiRC_dataXsec_err    = -1.        #statistical error on data Xsec (radiatively corrected using the Laget FSI model)
   else:
      fsiRC_dataXsec        = bin_info_fsiRC_dataXsec.cont[i]
      fsiRC_dataXsec_err    = bin_info_fsiRC_dataXsec.dcont[i]

   if(cont2[i] == 0.):
      fsiXsec               = -1.
      fsiXsec_err           = -1.
   else:
      fsiXsec               = bin_info_fsiXsec.cont[i]
      fsiXsec_err           = bin_info_fsiXsec.dcont[i]

   if(cont3[i] == 0.):
      pwiaRC_dataXsec       = -1.        #data Xsec radiatively corrected using the Laget PWIA model    
      pwiaRC_dataXsec_err   = -1.        #statistical error on data Xsec (radiatively corrected using the Laget PWIA model)
   else:
      pwiaRC_dataXsec       = bin_info_pwiaRC_dataXsec.cont[i]
      pwiaRC_dataXsec_err   = bin_info_pwiaRC_dataXsec.dcont[i]

   if(cont4[i] == 0.):
      pwiaXsec               = -1.
      pwiaXsec_err           = -1.
   else:
      pwiaXsec               = bin_info_pwiaXsec.cont[i]
      pwiaXsec_err           = bin_info_pwiaXsec.dcont[i]
      

   l = "%i %i %i %f %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n"%(ib, i_xbin[i], i_ybin[i], thnq_b[i], pm_b[i],  pwiaXsec, pwiaXsec_err, fsiXsec, fsiXsec_err, pwiaRC_dataXsec, pwiaRC_dataXsec_err, fsiRC_dataXsec, fsiRC_dataXsec_err)
                                                                          
   o.write(l)
o.close()


