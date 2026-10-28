import numpy as np
import run_command as rc
# import bin_info2 as BI
import analyze_differ as ad
from LT.datafile import dfile
# import LT.box as B
# import argparse as AG
import sys
# from sys import argv
# import os
import shutil as SU

#load module to calculate norm. systematics
# import norm_systematics as ns

dtr = np.pi/180.

#User Input
pm_set = int(sys.argv[1])
model = sys.argv[2]
# data_set = int(sys.argv[3])
# sys_ext = sys.argv[4]

#Code usage example: ipython run_differ_general_v2.py 580 fsi 1 sys_ext


#--------------CALCULATE NORM. SYSTEMATICS ERROR ON Xsec (Used External Module 'norm_systematics')-----------------

#This code will write the systematic contribution and the quadrature sum of normalization factors as relative error, df/f
#To get the absolute error to plot, one would have to multiply by the cross section
# ns.get_average_systematics(sys_ext)

#------------------------------------------------------------------------------------------------------------------


#----------------------------CALCULATE KINEMATIC SYSTEMATICS ERRORS ON Xsec-----------------------

#average kinematics directory for avg. kin with the nominal cut settings
avg_kin_dir = './yield_n_xsec/average_kin/'   #'./avg_kin_dir/'
results_dir = './SYSTEMATIC STUDIES/kin_systematics/kin_syst_results/'

#create input/output filenames for average kinematics to be read/ to write out systematics results
input_fname = 'pm_%i_%s_norad_avgkin.txt'%(pm_set, model)
input_fname_syst = 'pm_%i_%s_norad_avgkin_systematics.txt'%(pm_set, model)   #average kinematics file copy (to write systematic results)
input_file = avg_kin_dir+input_fname
input_file_syst =  avg_kin_dir+input_fname_syst
output_file = results_dir+'pm_%i_%s_results.txt'%(pm_set, model)

#Clean summary files (since they open in append mode inside the loop, we do not want to append past data)
fsumin = open('./SYSTEMATIC STUDIES/kin_systematics/summary_files/SummaryInput_pm_%i_%s.txt'%(pm_set, model), 'w').close()
fsumout = open('./SYSTEMATIC STUDIES/kin_systematics/summary_files/SummaryOutput_pm_%i_%s.txt'%(pm_set, model), 'w').close()

#Create Summary table to write derivatives for each (Pm, thnq) bin
#This file will be used in analyze_differ.py (by write_table() function), but it must be passed to get_sig_tot() first
ftable_name= './SYSTEMATIC STUDIES/kin_systematics/summary_files/DervTable_pm_%i_%s.txt'%(pm_set, model)
ftable_out = open(ftable_name, 'w')
ftable_out.write('#This table contains the kinematics derivatives, where ds_dkin is the derivative of the cross section\n')
ftable_out.write('#with respect to the kinematic variable in units of %/mrad or %/MeV\n')
ftable_out.write('#There are also the relative systematic errors of individual kinematic variables, and the total error added in quadrature\n')
ftable_out.write('#to see which is the dominant source of the systematics errors\n')

ftable_out.write('#!ib[i,0]/  pm_bin[f,1]/   thnq_bin[f,2]/   ds_dthe[f,3]/   ds_dphe[f,4]/   ds_def[f,5]/   ds_dthp[f,6]/   ds_dphp[f,7]/  ds_dthb[f,8]/   ds_dphb[f,9]/    ds_dE[f,10]/    sig_the[f,11]/   sig_phe[f,12]/   sig_ef[f,13]/   sig_thp[f,14]/   sig_php[f,15]/  sig_thb[f,16]/   sig_phb[f,17]/    sig_E[f,18]/    sig_the_thp[f,19]/    sig_the_Ef[f,20]/     sig_the_Eb[f,21]/      sig_thp_Ef[f,22]/      sig_thp_Eb[f,23]/      sig_Ef_Eb[f,24]/      sig_tot[f,25]/ \n')
ftable_out.close()

#Make a copy of the original avg. kin file (to write additional systematics headers)
SU.copyfile(input_file, input_file_syst)

#Open original avg. kin file
#kin_file = dfile(input_file)

#open newly created kin. file to add more keys
syst_file = dfile(input_file_syst)
syst_file.add_key('sig_the', 'f')
syst_file.add_key('sig_phe', 'f')
syst_file.add_key('sig_thp', 'f')
syst_file.add_key('sig_php', 'f')
syst_file.add_key('sig_thb', 'f')
syst_file.add_key('sig_phb', 'f')
syst_file.add_key('sig_Ef', 'f')
syst_file.add_key('sig_Eb', 'f')
syst_file.add_key('sig_the_thp', 'f')
syst_file.add_key('sig_the_Ef', 'f')
syst_file.add_key('sig_the_Eb', 'f')
syst_file.add_key('sig_thp_Ef', 'f')
syst_file.add_key('sig_thp_Eb', 'f')
syst_file.add_key('sig_Ef_Eb', 'f')
syst_file.add_key('sig_kin_tot', 'f')  #add new header to store total absolute kinematic systematic error per bin:   dsig = (dsig/sig) * sig_exp

#syst_file.save(input_file_syst)

#file to write input specectrometer kinematics for each (Pm, thnq) bin
#Only used temporarily for that bin, as it is overwritten when the
#next bin is read
differ_file = 'differ.in'

#get average of cos(phi_pq)
cont = np.array( syst_file.get_data('cont') )  #2d bin content 
cphi = np.array( syst_file.get_data('cos_phi') )
mask = (cont==cont) & (cont>0.)
cphi_av = np.average(cphi[mask])  #average ONLY over non-zero bin content
print('cphi_av = ',cphi_av)
# phi_off = 0.
# if cphi_av >= 0 :
#    phi_off = np.pi
   
theory_fname = f'yield_n_xsec/jml_theory_xsec/pm_{pm_set}_jml_theory.txt'
theory_file = dfile(theory_fname)

#Loop over 2D (Pm, thnq) bins to get avg. kin for each bin
for ik, ka in enumerate(syst_file.data):
   
   cont = ka['cont']       #2d bin content
   ib = ka['i_b']          #2d bin
   Ei = ka['Ei']           #Avg. incident energy (MeV) / kinematic bin
   Ef = Ei - ka['omega']   #Avg. final electron energy (MeV) / kinematic bin
   pm_bin = ka['yb']     #pm_bin center value (bin center +/- 20 MeV)
   thnq_bin = ka['xb']   #theta_nq bin center value (bin center +/- 5 deg) 
       

   ftable_bin = np.array([ib, pm_bin, thnq_bin])   #array with bin info to pass to write_table() function in analyze_differ.py

   #create differ input file name
   outf = open(differ_file,"w")
   outf.write(f'{ib}\n')
   outf.write('theta_e, phi_e, final electron energy, second line:step sizes in mr and MeV \n')
   l = '%f , 0., %f \n'%(ka['th_e'], Ef)
   outf.write(l)
   outf.write('1.,1.,0.1\n')
   outf.write('theta_p, phi_p (dont change sign for other side) theta=hor, phi = vert \n')
   
   # calculate php for the proton. Use the average value of phi
   # an offst of pi (or 180 degrees) needs to be subtracted for phi centered around 0 degrees
   # coordinate system used: beam direction: z
   #                         left side     : x
   #                         vertical up   : y
   # this currently wotks only for in-plane settings
   # calculations in calc_kin.C
   #
   # Double_t p_par = ( pow(pf,2) + pow(q_lab,2) - pow(pmc, 2))/ (2.*q_lab) ;
   # Double_t p_perp2 = abs(pow(pf, 2) - pow(p_par, 2));
   # Double_t p_perp = sqrt(p_perp2);
   pf = ka['pf']
   q_lab = ka['q_lab']
   pmc = ka['pm_mc']
   #angle of q with resp. to beam
   cthq = (Ei - Ef*np.cos( ka['th_e']*dtr ) )/q_lab
   theta_q = np.arccos(cthq)
   #direction of the proton
   #phi_v = ka['phi']*dtr - phi_off   #'phi' is actually phi_pq (I dont have 'phi' directly in kin file)
   cosphi = ka['cos_phi']    #cos(phi_pq)
   sinphi = ka['sin_phi']    #sin(phi_pq)
   phi_off = 0.
   if cosphi >= 0.:
       phi_off = np.pi    
   #fix this: get correct avg phi
   # phi_v_3 = np.arcsin(sinphi)
   phi_v = np.atan2(sinphi,cosphi) # avg phi angle
   # phi_v = theory_file['phi_avg'][ik]*dtr
   phi_v = np.arccos(ka['cos_phi'])  #arccos(x) is returned in radians (no need to convert to radians)
   # if ik == 90 or ik == 91: print('phi_avg = ',phi_avg,'\nphi_v = ',phi_v,'\nphi_v2 = ',phi_v_2,'\nphi_v3 = ',phi_v_3)
   # phi_v = phi_v - phi_off   #'phi' is actually phi_pq
   # phi_v = (phi_v + 2*np.pi) % 2*np.pi  # change phi_v from (-pi,pi) to (0,2pi)
   
   p_par = (pf**2 + q_lab**2 - pmc**2)/(2.*q_lab)
   p_per = np.sqrt(pf**2 - p_par**2)
   #find phi
   py = p_per*np.sin(phi_v)
   sphp = py/pf
   php = np.arcsin(sphp)/dtr   #this is the actual out-of-plane proton angle
   
   #find theta
   p_xz = np.sqrt(pf**2 - py**2)
   if (p_par == 0.):
      cthpp = 0.
   else:
      cthpp = p_par/p_xz
   thpp = np.arccos(cthpp)
   if phi_v == np.pi/2.:
      thp = theta_q/dtr
   elif phi_v < np.pi/2.:
      thp = (theta_q - thpp)/dtr
   else:
      thp = (theta_q + thpp)/dtr
   l = '%f, %f \n'%(thp, php)
   outf.write(l)
   outf.write('1.,1.\n')
   outf.write('theta_b, phi_b, Incident Energy (theta = vertical, phi = hor !)\n')
   l = '0., 0., %f \n'%(ka['Ei'])
   outf.write(l)
   outf.write('1., 1., 0.1 \n')
   outf.write('omega independent of E i=1 (energy loss), i=0 otherwise \n')
   outf.write('0\n')
   #output files
   outf.close()   

   outf = open('./differ.in', 'r')
   data_summary= outf.read()
   outf.close()
   fsumin = open('./SYSTEMATIC STUDIES/kin_systematics/summary_files/SummaryInput_pm_%i_%s.txt'%(pm_set, model), 'a')
   fsumin.write(data_summary)
   fsumin.write('----------------------------------------------------------------------------------\n')
   fsumin.close()

   #ready to run differ
   run_differ = './SYSTEMATIC STUDIES/kin_systematics/differ3'
   rc.run_command(run_differ, 'differ.out','differ_err')
   # now analyze differ, output tot_err is in %
   sigma_the, sigma_phe, sigma_thp, sigma_php, sigma_thb, sigma_phb, sigma_ef, sigma_dE, sigma_the_thp, sigma_the_Ef, sigma_the_Eb, sigma_thp_Ef, sigma_thp_Eb,  sigma_Ef_Eb, tot_err = ad.get_sig_tot('differ.out',  ftable_name, ftable_bin, print_all = False)
   if np.isnan(tot_err):
      print ('total_err could not be calculated, it is set to 0 for kinematics : ', syst_file, ' bin ib = ', ib)
      tot_err = 0.0      


   
   #print('norm_syst=',ka['dsig_norm_tot'])
   #print('IB>>>>>=',ka['i_b'])

   #print('kin_syst=',relative_kin_err)
   
   #Add information to kinematic systematics data file (all errors are in %)
   ka['sig_the'] = sigma_the
   ka['sig_phe'] = sigma_phe
   ka['sig_thp'] = sigma_thp
   ka['sig_php'] = sigma_php
   ka['sig_thb'] = sigma_thb
   ka['sig_phb'] = sigma_phb
   ka['sig_Ef'] = sigma_ef
   ka['sig_Eb'] = sigma_dE
   ka['sig_the_thp'] = sigma_the_thp
   ka['sig_the_Ef'] =  sigma_the_Ef
   ka['sig_the_Eb'] = sigma_the_Eb
   ka['sig_thp_Ef'] = sigma_thp_Ef
   ka['sig_thp_Eb'] = sigma_thp_Eb
   ka['sig_Ef_Eb'] = sigma_Ef_Eb
   ka['sig_kin_tot'] = float(tot_err)
   
   #Add normalization and systematics errors in quadrature
   #tot_syst_err =np.sqrt(relative_kin_err**2 + ka['dsig_norm_tot']*ka['dsig_norm_tot']) 

   #ka['dsig_tot_syst'] = tot_syst_err



   fout_differ = open('./differ.out')
   data_summary= fout_differ.read()
   fsumout = open('./SYSTEMATIC STUDIES/kin_systematics/summary_files/SummaryOutput_pm_%i_%s.txt'%(pm_set, model), 'a')
   fsumout.write('2DBIN #: %i'%(ib))
   fsumout.write(data_summary)
   fsumout.write('----------------------------------------------------------------------------------\n')
   fsumout.close()

syst_file.save(input_file_syst)
syst_file.add_header_comment(' ------------------------------------------------------ ')
syst_file.add_header_comment(' All systematic errors are written in percent.')
syst_file.add_header_comment(' They are the relative error in cross section dsig/sig')
syst_file.add_header_comment(' To get the absolute error on the cross section, multiply by the cross section value. dsig/sig * sig_exp')
syst_file.add_header_comment(' ------------------------------------------------------ ')
syst_file.save(input_file_syst)


#write the new data into a file
#results = open(syst_file, 'w')
#syst_file.add_header_comment(' ------------------------- ')
#syst_file.add_header_comment(' systematic errors are in % !')
#syst_file.add_header_comment(' ------------------------- ')
#syst_file.write_all(results, complete_header = True)

