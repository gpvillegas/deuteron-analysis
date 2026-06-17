#Calculate the Laget Theory Cross Section at the Averaged Kinematics

# select of reading binary grid data
use_binary = 1

# import laget module
if use_binary:
    import Laget_Xsec_fp_1 as LX
else:
    import Laget_Xsec_fp as LX
        
import LT.box as B
import numpy as np
import pandas as pd
from deForest_Xsec import GMp,GEp,sigMott,deForest  #import deForest sig_cc1, mott, GMp, GEp functs.
#%%
# import sys
# import os
# from sys import argv

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr

# some constants
dtr = np.pi/180.
#MeV
MP = 938.272
MN = 939.566
#MD = 1875.6127
MD = 1875.61
me = 0.51099 

fm2ub = 1.e4   #fm^2 to ub
# initalize arrays: there is a link in the current directory
# assume it is in the current directory
deut_data_dir = './'

# select what kind of calculation and set flags
#do_fsi = 0
#do_pwia = 1

# select linear interpolation
lin_interp = 2
# do not save the grid as binary file
save_grid = 0
#%%

#------------------------------------------------------------
# header information for the output file
header = \
"""
# Laget Cross Section Results at the Averaged Kinematics
# averaged kinematic varibles used as input to calculate the cross section 
# kinematics: Ei, Q2, omega, th_pq_cm, phi_pq
#
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/  pwiaXsec[f,5]/  pwiaGEp[f,6]/   pwiaGMp[f,7]/   pwia_sigMott[f,8]/   pwia_Ksig_cc1[f,9]/    
"""
#------------------------------------------------------------
#%%
#create output file to write cross section @ avg kinematics

#usage: ipython calc_cross.py 580 1 sys_ext

#User Input
pm_set = input('Enter pm setting: \n')
# pm_set = int(sys.argv[1])
# data_set = int(sys.argv[2])
# sys_ext = sys.argv[3]

#check if directory exists, else creates it.
# if not os.path.exists(sys_ext):
#     os.makedirs(sys_ext)

# # print argv
# if pm_set == 80:
#    output_file = './%s/pm%i_laget_theory.txt'%(sys_ext, pm_set)
# else:
output_file = f'./yield_n_xsec/jml_theory_xsec/{pm_set}_jml_theory.txt'

o = open(output_file,'w')
# write header
o.write(header)

# intitialize the code
#if use_binary:
#    LX.init_laget('./', do_fsi, lin_interp, save_grid, use_binary)
#else:
#    LX.init_laget('./', do_fsi, lin_interp)


#======================================DO PWIA=============================================================

# read the averaged kinematics file
# if pm_set == 80:
#     f = B.get_file('../average_kinematics/%s/pm%i_pwia_norad_avgkin.txt'%(sys_ext,pm_set))
# else: 
pwia_fileloc = f'./yield_n_xsec/average_kin/{pm_set}_jmlpwia_norad_avgkin.txt'
print(f'**Getting PWIA data from\n{pwia_fileloc}...')
f = B.get_file(pwia_fileloc)


#Read Headers to be written to the output Xsec file
i_b = B.get_data(f, 'i_b')      #ith (Pm, th_nq) bin number
i_x = B.get_data(f, 'i_x')      #ith th_nq bin number
i_y = B.get_data(f, 'i_y')      #ith Pm bin number
xb  = B.get_data(f, 'xb')       #theta_nq central bin value [deg]
yb  = B.get_data(f, 'yb')       #Pmiss cental bin value [GeV]
cont = B.get_data(f, 'cont')    #2D bin content

#===============================================================================
#=========  Laget_Xsec_fp_1.f,  get_sigma_laget() Input Parameters =============
#===============================================================================
#units are in MeV and deg [convert to radians, as Laget code takes MeV and radians]

th_cm = convert2NaN(pd.Series(B.get_data(f, 'th_pq_cm'))*dtr)     #center of mass in-plane angle between proton and q-vector [rad]
omega = pd.Series(B.get_data(f,  'omega'))           #averaged calculated energy transfer [MeV]
q = pd.Series(B.get_data(f, 'q_lab'))                #averaged magnitude of 3-momentum transfer [MeV]
Q2 = pd.Series(B.get_data(f, 'Q2_calc'))             #averaged calculated 4-Momentum Transfer [MeV^2]
cphi = pd.Series(B.get_data(f, 'cos_phi'))           #averaged Cos(phi_pq) : out-of-plane angle between proton and q-vector, determined from SIMC  vertex quantities 
sphi = pd.Series(B.get_data(f, 'sin_phi'))           #averaged Sin(phi_pq) : out-of-plane angle between proton and q-vector, determined from SIMC  vertex quantities 
phi = np.atan2(sphi,cphi)                           # average phi calculated from the average sin and cos phi in radians
phi_wrong = np.arccos(cphi)                    #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
Ei = pd.Series(B.get_data(f, 'Ei'))                  #averaged Incident Beam Energy [MeV]
GE_p = B.get_data(f, 'GEp')                #Proton Electric Form Factor at the avg. kinematics
GM_p = B.get_data(f, 'GMp')                #Proton Magnetic Form Factor at the avg. kinematics
sig_Mott = B.get_data(f, 'sigMott')        #Mott cross section at the avg. kinematics
Ksig_cc1 = B.get_data(f, 'Ksig_cc1')      #Kinematic Factor K * deForest cc1 cross section at the avg. kinematics
kf = pd.Series(B.get_data(f, 'kf'))
the = pd.Series(B.get_data(f, 'th_e'))
Pf = pd.Series(B.get_data(f, 'pf'))
Pm = pd.Series(B.get_data(f, 'pm'))
# thp = pd.Series(B.get_data(f, 'th_p'))
thpq_calc = pd.Series(B.get_data(f, 'theta_pq_calc'))

# do a bit of interpolation to try to fill in gaps the data has
th_cm = th_cm.interpolate(method='cubicspline', limit=3, limit_area='inside')
omega = omega.interpolate(method='cubicspline', limit=3, limit_area='inside')
Q2 = Q2.interpolate(method='cubicspline', limit=3, limit_area='inside')
phi = phi.interpolate(method='cubicspline', limit=3, limit_area='inside')
Ei = Ei.interpolate(method='cubicspline', limit=3, limit_area='inside')
kf = kf.interpolate(method='cubicspline', limit=3, limit_area='inside')
the = the.interpolate(method='cubicspline', limit=3, limit_area='inside')
Pf = Pf.interpolate(method='cubicspline', limit=3, limit_area='inside')
Pm = Pm.interpolate(method='cubicspline', limit=3, limit_area='inside')
# thp = thp.interpolate(method='cubicspline', limit=3, limit_area='inside')
thpq_calc = thpq_calc.interpolate(method='cubicspline', limit=3, limit_area='inside')

# DO PWIA
print('**Initializing Laget...')
if use_binary:
    LX.init_laget('./', 0, lin_interp, save_grid, use_binary)
else:
    LX.init_laget('./', 0, lin_interp)
    
print('**Starting PWIA loop...')
sigma = []
for i, e0_i in enumerate(Ei):
    # check kinematics
    # if e0_i != 0.0:
    #     print("i = ",i," e0_i= ",e0_i, " Q2= ",Q2[i]," nu= ", omega[i]," th_cm=",th_cm[i], " phi=",phi[i])
    #     print ("p_rec = ", LX.p_recoil(Q2[i], omega[i], th_cm[i]))
    if i%100 == 0:
        print(f'**Entry: {i}')
    if cont[i]==0.0:
    # if False:
        sig         = -1.
        GE_p[i]      = -1.
        GM_p[i]      = -1.
        sig_Mott[i]  = -1.
        Ksig_cc1[i] = -1.
        #sigma.append(sig)
        # print ('GEp, GMp, sigMott, Ksig_cc1 = ', GEp[i], GMp[i], sigMott[i], Ksig_cc1[i])
        # Ef = np.sqrt(kf[i]*kf[i] + me*me)  
        
        # sig = LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi[i]) * fm2ub    #convert to (ub sr^-2 MeV^-1)
        
        #Calculate the deForest Cross Section Factor  K * sig_cc1,  where K = Pf * Ep , units: ub * MeV^2 / sr^2
        #from interpolated values
        # sig_Mott =  sigMott(kf[i], the[i], Q2[i])   #ub / sr
        # GE_p = GEp(Q2[i], 'JRA')
        # GM_p = GMp(Q2[i], 'JRA')
        # Kfact, f_rec, sig_eN, de_Forest = deForest(Ef, Q2[i], q[i], Pf[i], Pm[i], the[i], thp[i], cphi[i], thpq_calc[i], sig_Mott, GE_p, GM_p)
        
        # gep[i] = GE_p ; gmp[i] = GM_p ; sig_mott[i] = sig_Mott ; Ksig_cc1[i] = de_Forest
        
        if i == 100:
            print('xb = ', xb[i])
            print('yb = ', yb[i])
            print('e0_i, Q2, om, thpq, phi = ',e0_i,Q2[i],omega[i],th_cm[i],phi[i])
            print('sig = ', sig)
    else:
        # calculate cross sections (fm^2 sr^-2 mev^-1)
        sig = LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi[i]) * fm2ub    #convert to (ub sr^-2 MeV^-1)
        #sigma.append(sig)
        # if e0_i != 0.0:
        #     print('sig = ',sig)
        if yb[i] == 45. :
            print('sig = ',sig)
        # Write to File
    l = '%i %i %i %.6f %.6f %.12e  %.6f  %.6f  %.6f  %.6f\n'%(i_b[i], i_x[i], i_y[i], xb[i], yb[i], sig, GE_p[i], GM_p[i], sig_Mott[i], Ksig_cc1[i])
    o.write(l)
o.close()

#======================================DO FSI================================================

# read the averaged kinematics file, and output file to write cross sections
# if pm_set == 80:
#     f = B.get_file('../average_kinematics/%s/pm%i_fsi_norad_avgkin.txt'%(sys_ext, pm_set))
#     fname = './%s/pm%i_laget_theory.txt'%(sys_ext, pm_set)
#     output_file = B.get_file(fname)
# else: 
f = B.get_file(f'./yield_n_xsec/average_kin/{pm_set}_jmlfsi_norad_avgkin.txt')
fname = f'./yield_n_xsec/jml_theory_xsec/{pm_set}_jml_theory.txt'
output_file = B.get_file(fname)

print('**Preparing output file for FSI data...')
#Add Key to the output_file header
output_file.add_key('fsiXsec', 'f')
output_file.add_key('fsiGEp', 'f')
output_file.add_key('fsiGMp', 'f')
output_file.add_key('fsi_sigMott', 'f')
output_file.add_key('fsi_Ksig_cc1', 'f')
#Add additional avg. kinematics to output_file header
output_file.add_key('Ei_avg')     
output_file.add_key('th_pq_cm')   
output_file.add_key('omega_avg')  
output_file.add_key('q_avg')      
output_file.add_key('Q2_avg')     
output_file.add_key('cphi_pq_avg')
output_file.add_key('sphi_pq_avg')
output_file.add_key('phi_avg')
output_file.add_key('phi_wrong')
output_file.add_key('pf_avg')     
output_file.add_key('pm_avg')     
output_file.add_key('kf_avg')     
output_file.add_key('the_avg')    
output_file.add_key('xbj_avg')    








#NOTE: The PWIA and FSI avg. kin. files bin-numbering must be EXACTLY the same. Therefore, there is no need to read them a second time.
cont = B.get_data(f, 'cont')   #the bin content, however, might be different so it is read

#===============================================================================
#=========  Laget_Xsec_fp_1.f,  get_sigma_laget() Input Parameters =============
#===============================================================================
#units are in MeV and deg [convert to radians, as Laget code takes MeV and radians]
#All these quantities are averaged from SIMC or they have been calculated using the assumed masses
th_cm = B.get_data(f, 'th_pq_cm')*dtr     #center of mass in-plane angle between proton and q-vector [rad]
omega = B.get_data(f,  'omega')           #energy transfer [MeV]
q = B.get_data(f, 'q_lab')                #magnitude of 3-momentum transfer [MeV]
Q2 = B.get_data(f, 'Q2_calc')             #4-Momentum Transfer [MeV^2]
cphi = B.get_data(f, 'cos_phi')           #Cos(phi_pq)  out of plane angle between the proton and q-vector (lab frame)
sphi = pd.Series(B.get_data(f, 'sin_phi'))           #averaged Sin(phi_pq) : out-of-plane angle between proton and q-vector, determined from SIMC  vertex quantities 
phi = np.atan2(sphi,cphi)                           # average phi calculated from the average sin and cos phi in radians
phi_wrong = np.arccos(cphi)                    #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
Ei = B.get_data(f, 'Ei')                  #Incident Beam Energy [MeV]
GE_p = B.get_data(f, 'GEp')                #Proton Electric Form Factor at the avg. kinematics
GM_p = B.get_data(f, 'GMp')                #Proton Magnetic Form Factor at the avg. kinematics
sig_Mott = B.get_data(f, 'sigMott')        #Mott cross section at the avg. kinematics
Ksig_cc1 = B.get_data(f, 'Ksig_cc1')      #deForest cc1 cross section at the avg. kinematics
#Additional averaged kinematics to add
Pf = B.get_data(f, 'pf')                  #averaged proton final momentum
pm_avg = B.get_data(f, 'pm')              #averaged calculated missing momentum assuming deuteron proton and neutron mass
kf = B.get_data(f, 'kf')                  #avergaed final e- momentum
th_e = B.get_data(f, 'th_e')               #averaged in-plane electron angle
xbj = Q2 / (2.*MP*omega)                  #averaged x-Bkorker determined from calculated averaged Q2 and omega

# DO FSI
print('**Initializing Laget...')
if use_binary:
    LX.init_laget('./', 1, lin_interp, save_grid, use_binary)
else:
    LX.init_laget('./', 1, lin_interp)

print('**Starting FSI loop...')
#sigma = []
for i, e0_i in enumerate(Ei):
    # check kinematics
    #print("i = ",i," e0_i= ",e0_i, " Q2= ",Q2[i]," nu= ", omega[i]," th_cm=",th_cm[i], " phi=",phi[i])
    #print "p_rec = ", LX.p_recoil(Q2[i], omega[i], th_cm[i])
    #create new headers and set them to -1.
    if i%100 == 0:
        print(f'**Entry: {i}')
    
    if cont[i]==0.0:
    # if False:
        sig = -1.   
        output_file.data[i]['fsiXsec']      = -1.  #theory FSI Xsec (Paris)
        output_file.data[i]['fsiGEp']       = -1.  #proton Electric Form Factor
        output_file.data[i]['fsiGMp']       = -1.  #proton Magnetic Form Factor
        output_file.data[i]['fsi_sigMott']  = -1.  #mott cross section
        output_file.data[i]['fsi_Ksig_cc1'] = -1.  #K * f_rec *sig_cc1 
        #Write averaged kinematics
        output_file.data[i]['Ei_avg']      = -1.  #averaged incident beam energy  [MeV] 
        output_file.data[i]['th_pq_cm']      = -1.  #averaged th_pq in cm frame   [deg]
        output_file.data[i]['omega_avg']        = -1. #average energy transfer    [MeV]
        output_file.data[i]['q_avg']         = -1.  #average 3-momentum transfer, |q|  [MeV]
        output_file.data[i]['Q2_avg']        = -1.   #avg. calc. 4-momentum transfer [MeV^2]
        output_file.data[i]['cphi_pq_avg']    = -1.  #avg. cos(ph_pq) in the lab frame 
        output_file.data[i]['sphi_pq_avg']  = -1
        output_file.data[i]['phi_avg']      = -1
        output_file.data[i]['phi_wrong']    = -1
        output_file.data[i]['pf_avg']    = -1.      #avg. final proton momentum [MeV]
        output_file.data[i]['pm_avg']    = -1.     #avg. missing momentum [MeV]
        output_file.data[i]['kf_avg']    = -1.     #avg. final e- momentum [MeV]
        output_file.data[i]['the_avg']    = -1.     #avg. final e- angle [deg]
        output_file.data[i]['xbj_avg']    = -1.     #avg. X-Bjorken 

        

    else:
        # calculate cross sections (returns Xsec in fm^2 sr^-2 MeV^-1 (SIMC is in microbarn.  1 ub = 1e-4 fm^2)
        sig = LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi[i]) * fm2ub    #convert to (ub sr^-2 MeV^-1)
        output_file.data[i]['fsiXsec'] = float("%.12e"%(sig))                    #ub sr^-2 MeV^-1
        output_file.data[i]['fsiGEp'] = float("%.6f"%(GE_p[i]))
        output_file.data[i]['fsiGMp'] = float("%.6f"%(GM_p[i]))
        output_file.data[i]['fsi_sigMott'] = float("%.6f"%(sig_Mott[i]))           #ub / sr
        output_file.data[i]['fsi_Ksig_cc1'] = float("%.6f"%(Ksig_cc1[i]))         # ub MeV^2 / sr^2
        #Write averaged kinematics quantities 
        output_file.data[i]['Ei_avg']     = float("%.6f"%(Ei[i]))                    
        output_file.data[i]['th_pq_cm']   = float("%.6f"%(th_cm[i]))/dtr     
        output_file.data[i]['omega_avg']  = float("%.6f"%(omega[i]))     
        output_file.data[i]['q_avg']      = float("%.6f"%(q[i]))
        output_file.data[i]['Q2_avg']     = float("%.6f"%(Q2[i]))
        output_file.data[i]['cphi_pq_avg']= float("%.6f"%(cphi[i]))
        output_file.data[i]['sphi_pq_avg']  = float("%.6f"%(sphi[i]))
        output_file.data[i]['phi_avg']      = float("%.6f"%(phi[i]))/dtr
        output_file.data[i]['phi_wrong']    = float("%.6f"%(phi_wrong[i]))/dtr
        output_file.data[i]['pf_avg']     = float("%.6f"%(Pf[i]))
        output_file.data[i]['pm_avg']     = float("%.6f"%(pm_avg[i]))
        output_file.data[i]['kf_avg']     = float("%.6f"%(kf[i]))
        output_file.data[i]['the_avg']    = float("%.6f"%(th_e[i]))
        output_file.data[i]['xbj_avg']    = float("%.6f"%(xbj[i]))
    
    # if i == 99:
    #     print(i_b[i])
    #     print('xb = ', xb[i])
    #     print('yb = ', yb[i])
    #     print('Ei, Q2, omega, th_cm, phi, pm_avg = ',Ei[i],Q2[i],omega[i],
    #           th_cm[i]/dtr,phi[i],pm_avg[i])
        
# Write to File
output_file.save(fname)

#%%
# Calculate Laget Xsec using interpolated values
pm_set = input('Enter pm setting: \n')
output_file = f'./yield_n_xsec/jml_theory_xsec/{pm_set}_jml_theory_interp.txt'
o = open(output_file,'w')
o.write(header)

pwia_fileloc = f'./yield_n_xsec/average_kin/{pm_set}_jmlpwia_norad_avgkin.txt'
print(f'**Getting PWIA data from\n{pwia_fileloc}...')
f = B.get_file(pwia_fileloc)


#Read Headers to be written to the output Xsec file
i_b = B.get_data(f, 'i_b')      #ith (Pm, th_nq) bin number
i_x = B.get_data(f, 'i_x')      #ith th_nq bin number
i_y = B.get_data(f, 'i_y')      #ith Pm bin number
xb  = B.get_data(f, 'xb')       #theta_nq central bin value [deg]
yb  = B.get_data(f, 'yb')       #Pmiss cental bin value [GeV]
cont = B.get_data(f, 'cont')    #2D bin content

#===============================================================================
#=========  Laget_Xsec_fp_1.f,  get_sigma_laget() Input Parameters =============
#===============================================================================
#units are in MeV and deg [convert to radians, as Laget code takes MeV and radians]

th_cm = convert2NaN(pd.Series(B.get_data(f, 'th_pq_cm'))*dtr)     #center of mass in-plane angle between proton and q-vector [rad]
omega = pd.Series(B.get_data(f,  'omega'))           #averaged calculated energy transfer [MeV]
q = pd.Series(B.get_data(f, 'q_lab'))                #averaged magnitude of 3-momentum transfer [MeV]
Q2 = pd.Series(B.get_data(f, 'Q2_calc'))             #averaged calculated 4-Momentum Transfer [MeV^2]
cphi = pd.Series(B.get_data(f, 'cos_phi'))           #averaged Cos(phi_pq) : out-of-plane angle between proton and q-vector, determined from SIMC  vertex quantities 
phi = np.arccos(cphi)                    #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
Ei = pd.Series(B.get_data(f, 'Ei'))                  #averaged Incident Beam Energy [MeV]
GE_p = B.get_data(f, 'GEp')                #Proton Electric Form Factor at the avg. kinematics
GM_p = B.get_data(f, 'GMp')                #Proton Magnetic Form Factor at the avg. kinematics
sig_Mott = B.get_data(f, 'sigMott')        #Mott cross section at the avg. kinematics
Ksig_cc1 = B.get_data(f, 'Ksig_cc1')      #Kinematic Factor K * deForest cc1 cross section at the avg. kinematics
kf = pd.Series(B.get_data(f, 'kf'))
the = pd.Series(B.get_data(f, 'th_e'))
Pf = pd.Series(B.get_data(f, 'pf'))
Pm = pd.Series(B.get_data(f, 'pm'))
thp = pd.Series(B.get_data(f, 'th_p'))
thpq_calc = pd.Series(B.get_data(f, 'theta_pq_calc'))

# DO PWIA
print('**Initializing Laget...')
if use_binary:
    LX.init_laget('./', 0, lin_interp, save_grid, use_binary)
else:
    LX.init_laget('./', 0, lin_interp)
    
print('**Starting PWIA loop...')
sig = []
for i, e0_i in enumerate(Ei):
    # check kinematics
    # if e0_i != 0.0:
    #     print("i = ",i," e0_i= ",e0_i, " Q2= ",Q2[i]," nu= ", omega[i]," th_cm=",th_cm[i], " phi=",phi[i])
    #     print ("p_rec = ", LX.p_recoil(Q2[i], omega[i], th_cm[i]))
    if i%100 == 0:
        print(f'**Entry: {i}')
    if cont[i]==0.0:
    # if False:
        sig.append(np.nan)
        GE_p[i]      = np.nan
        GM_p[i]      = np.nan
        sig_Mott[i]  = np.nan
        Ksig_cc1[i] = np.nan
        #sigma.append(sig)
        # print ('GEp, GMp, sigMott, Ksig_cc1 = ', GEp[i], GMp[i], sigMott[i], Ksig_cc1[i])
        
    else:
        # calculate cross sections (fm^2 sr^-2 mev^-1)
        sig.append(LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi[i]) * fm2ub)    #convert to (ub sr^-2 MeV^-1)
        #sigma.append(sig)
        # if e0_i != 0.0:
        #     print('sig = ',sig)
        # Write to File
    
sig = pd.Series(sig)
sig_Mott = pd.Series(sig_Mott)
Ksig_cc1 = pd.Series(Ksig_cc1)

sig = sig.interpolate()
sig_Mott = sig_Mott.interpolate()
Ksig_cc1 = Ksig_cc1.interpolate()

for i, e0_i in enumerate(Ei):    
    l = '%i %i %i %.6f %.6f %.12e  %.6f  %.6f  %.6f  %.6f\n'%(i_b[i], i_x[i], i_y[i], xb[i], yb[i], sig[i], GE_p[i], GM_p[i], sig_Mott[i], Ksig_cc1[i])
    o.write(l)
o.close()

#======================================DO FSI================================================

f = B.get_file(f'./yield_n_xsec/average_kin/{pm_set}_jmlfsi_norad_avgkin.txt')
fname = f'./yield_n_xsec/jml_theory_xsec/{pm_set}_jml_theory_interp.txt'
output_file = B.get_file(fname)

print('**Preparing output file for FSI data...')
#Add Key to the output_file header
output_file.add_key('fsiXsec', 'f')
output_file.add_key('fsiGEp', 'f')
output_file.add_key('fsiGMp', 'f')
output_file.add_key('fsi_sigMott', 'f')
output_file.add_key('fsi_Ksig_cc1', 'f')
#Add additional avg. kinematics to output_file header
output_file.add_key('Ei_avg')     
output_file.add_key('th_pq_cm')   
output_file.add_key('omega_avg')  
output_file.add_key('q_avg')      
output_file.add_key('Q2_avg')     
output_file.add_key('cphi_pq_avg')
output_file.add_key('pf_avg')     
output_file.add_key('pm_avg')     
output_file.add_key('kf_avg')     
output_file.add_key('the_avg')    
output_file.add_key('xbj_avg')    

#NOTE: The PWIA and FSI avg. kin. files bin-numbering must be EXACTLY the same. Therefore, there is no need to read them a second time.
cont = B.get_data(f, 'cont')   #the bin content, however, might be different so it is read

#===============================================================================
#=========  Laget_Xsec_fp_1.f,  get_sigma_laget() Input Parameters =============
#===============================================================================
#units are in MeV and deg [convert to radians, as Laget code takes MeV and radians]
#All these quantities are averaged from SIMC or they have been calculated using the assumed masses
th_cm = B.get_data(f, 'th_pq_cm')*dtr     #center of mass in-plane angle between proton and q-vector [rad]
omega = B.get_data(f,  'omega')           #energy transfer [MeV]
q = B.get_data(f, 'q_lab')                #magnitude of 3-momentum transfer [MeV]
Q2 = B.get_data(f, 'Q2_calc')             #4-Momentum Transfer [MeV^2]
cphi = B.get_data(f, 'cos_phi')           #Cos(phi_pq)  out of plane angle between the proton and q-vector (lab frame)
phi = np.arccos(cphi)                     #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
Ei = B.get_data(f, 'Ei')                  #Incident Beam Energy [MeV]
GE_p = B.get_data(f, 'GEp')                #Proton Electric Form Factor at the avg. kinematics
GM_p = B.get_data(f, 'GMp')                #Proton Magnetic Form Factor at the avg. kinematics
sig_Mott = B.get_data(f, 'sigMott')        #Mott cross section at the avg. kinematics
Ksig_cc1 = B.get_data(f, 'Ksig_cc1')      #deForest cc1 cross section at the avg. kinematics
#Additional averaged kinematics to add
Pf = B.get_data(f, 'pf')                  #averaged proton final momentum
pm_avg = B.get_data(f, 'pm')              #averaged calculated missing momentum assuming deuteron proton and neutron mass
kf = B.get_data(f, 'kf')                  #avergaed final e- momentum
th_e = B.get_data(f, 'th_e')               #averaged in-plane electron angle
xbj = Q2 / (2.*MP*omega)                  #averaged x-Bkorker determined from calculated averaged Q2 and omega

# DO FSI
print('**Initializing Laget...')
if use_binary:
    LX.init_laget('./', 1, lin_interp, save_grid, use_binary)
else:
    LX.init_laget('./', 1, lin_interp)

print('**Starting FSI loop...')
sig = []
for i, e0_i in enumerate(Ei):
    # check kinematics
    #print("i = ",i," e0_i= ",e0_i, " Q2= ",Q2[i]," nu= ", omega[i]," th_cm=",th_cm[i], " phi=",phi[i])
    #print "p_rec = ", LX.p_recoil(Q2[i], omega[i], th_cm[i])
    #create new headers and set them to -1.
    if i%100 == 0:
        print(f'**Entry: {i}')
    
    if cont[i]==0.0:
    # if False:
        sig.append(np.nan)   
        # output_file.data[i]['fsiXsec']      = -1.  #theory FSI Xsec (Paris)
        GE_p[i]      = np.nan
        GM_p[i]      = np.nan
        sig_Mott[i]  = np.nan
        Ksig_cc1[i] = np.nan
        #Write averaged kinematics
        output_file.data[i]['Ei_avg']      = -1.  #averaged incident beam energy  [MeV] 
        output_file.data[i]['th_pq_cm']      = -1.  #averaged th_pq in cm frame   [deg]
        output_file.data[i]['omega_avg']        = -1. #average energy transfer    [MeV]
        output_file.data[i]['q_avg']         = -1.  #average 3-momentum transfer, |q|  [MeV]
        output_file.data[i]['Q2_avg']        = -1.   #avg. calc. 4-momentum transfer [MeV^2]
        output_file.data[i]['cphi_pq_avg']    = -1.  #avg. cos(ph_pq) in the lab frame 
        output_file.data[i]['pf_avg']    = -1.      #avg. final proton momentum [MeV]
        output_file.data[i]['pm_avg']    = -1.     #avg. missing momentum [MeV]
        output_file.data[i]['kf_avg']    = -1.     #avg. final e- momentum [MeV]
        output_file.data[i]['the_avg']    = -1.     #avg. final e- angle [deg]
        output_file.data[i]['xbj_avg']    = -1.     #avg. X-Bjorken 

        

    else:
        # calculate cross sections (returns Xsec in fm^2 sr^-2 MeV^-1 (SIMC is in microbarn.  1 ub = 1e-4 fm^2)
        sig.append(LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi[i]) * fm2ub)    #convert to (ub sr^-2 MeV^-1)
        output_file.data[i]['fsiXsec'] = float("%.12e"%(sig[i]))                    #ub sr^-2 MeV^-1
        output_file.data[i]['fsiGEp'] = float("%.6f"%(GE_p[i]))
        output_file.data[i]['fsiGMp'] = float("%.6f"%(GM_p[i]))
        output_file.data[i]['fsi_sigMott'] = float("%.6f"%(sig_Mott[i]))           #ub / sr
        output_file.data[i]['fsi_Ksig_cc1'] = float("%.6f"%(Ksig_cc1[i]))         # ub MeV^2 / sr^2
        #Write averaged kinematics quantities 
        output_file.data[i]['Ei_avg']     = float("%.6f"%(Ei[i]))                    
        output_file.data[i]['th_pq_cm']   = float("%.6f"%(th_cm[i]))/dtr     
        output_file.data[i]['omega_avg']  = float("%.6f"%(omega[i]))     
        output_file.data[i]['q_avg']      = float("%.6f"%(q[i]))
        output_file.data[i]['Q2_avg']     = float("%.6f"%(Q2[i]))
        output_file.data[i]['cphi_pq_avg']= float("%.6f"%(cphi[i]))
        output_file.data[i]['pf_avg']     = float("%.6f"%(Pf[i]))
        output_file.data[i]['pm_avg']     = float("%.6f"%(pm_avg[i]))
        output_file.data[i]['kf_avg']     = float("%.6f"%(kf[i]))
        output_file.data[i]['the_avg']    = float("%.6f"%(th_e[i]))
        output_file.data[i]['xbj_avg']    = float("%.6f"%(xbj[i]))
    
    # if i == 99:
    #     print(i_b[i])
    #     print('xb = ', xb[i])
    #     print('yb = ', yb[i])
    #     print('Ei, Q2, omega, th_cm, phi, pm_avg = ',Ei[i],Q2[i],omega[i],
    #           th_cm[i]/dtr,phi[i],pm_avg[i])
sig = pd.Series(sig)
sig_Mott = pd.Series(sig_Mott)
Ksig_cc1 = pd.Series(Ksig_cc1)

sig = sig.interpolate(method='cubicspline',limit_area='inside')
sig_Mott = sig_Mott.interpolate(method='cubicspline',limit_area='inside')
Ksig_cc1 = Ksig_cc1.interpolate(method='cubicspline',limit_area='inside')

for i, e0_i in enumerate(Ei):
    output_file.data[i]['fsiXsec'] = float("%.12e"%(sig[i]))                    #ub sr^-2 MeV^-1
    output_file.data[i]['fsi_sigMott'] = float("%.6f"%(sig_Mott[i]))           #ub / sr
    output_file.data[i]['fsi_Ksig_cc1'] = float("%.6f"%(Ksig_cc1[i]))         # ub MeV^2 / sr^2
       
# Write to File
output_file.save(fname)

#%%

sig_neg = LX.get_sigma_laget(10598.800659,4137222.290039,1673.210129,10.368263,-177*dtr)
sig_pos = LX.get_sigma_laget(10598.800659,4137222.290039,1673.210129,10.368263,183*dtr)

print(sig_neg)
print(sig_pos)
