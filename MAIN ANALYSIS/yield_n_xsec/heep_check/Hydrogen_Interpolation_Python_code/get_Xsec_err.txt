import numpy as np
import matplotlib.pyplot as plt
from interpolate_ff import *

#This code defines the H(e,e'p) elastic cross section in terms of
#the electric and magnetic form factors. Then apply error propagation
#in terms of the form factors to get uncertainty in cross section at thw
#gven beam energy, elecrtron angle and momentum.
#Define Physical Constnats

hbarc = 197.327053 / 1000.       #GeV * fm
alpha = 1./137.0359895    #structure constant
dtr = np.pi / 180.

#Masses in GeV
MP = 938.272 / 1000.
MN = 939.566 / 1000.
MD = 1875.61 / 1000.
me = 0.51099 / 1000.

#Template function that returns the H(e,e'p) cross section and its absolute uncertainty in units of ub / sr
def sig(kf, th_e, Q2, GEp_GD, GMp_muGD, dGEp_GD, dGMp_muGD):

    #Units: kf ~ Ef [GeV]   electron momentum  
    #th_e[deg] electron angle  
    #Q2[GeV2]   4-momentum transfer

    mu_N = 3.1524512326 *1e-14 / 1000.   #nuclear magneton [GeV T^-1]
    mu_p = 2.7928473508 * mu_N   #proton magnetic moment [GeV T^-1]

    #Mott cross section
    th_e = th_e * dtr    #convert degree to radians
    sigM = (2.*alpha*hbarc*kf*np.cos(th_e/2.)/Q2 )**2   #GeV^2 * fm^2 *GeV^2/GeV^4--> fm2 / sr
    sigMott = sigM * 1e4    #convert fm^2 -> microbarn

    tau =  Q2 /(4.*MP)**2
    L2 = 0.71   #Lambda parameter [GeV^2]
    GD = 1. / (1 + Q2/L2)**2     #dipole form factor
    muGD = mu_p * GD   #[GeV T^-1]
    
    #Calculate the form factors
    GEp = GEp_GD * GD
    GMp = GMp_muGD * muGD

    #Elastic H(e,e'p) Cross Section
    sig = sigMott * ( (GEp**2 + tau*GMp**2)/(1. + tau) + 2.*tau*GMp**2*np.tan(th_e/2.)**2 )

    #Error Propagation
    #derivatives of cross sectio w.r.to form factors
    dsig_dGEp = sigMott * ( (2.*GEp + tau*GMp**2) / (1. + tau))
    dsig_dGMp = sigMott * ( (GEp**2 + 2.*tau*GMp) / (1. + tau) + 4.*tau*GMp*np.tan(th_e/2.)**2 )

    #uncertainty in form factors
    dGEp = dGEp_GD * GD
    dGMp = dGMp_muGD * muGD

    dsig2 =  (dsig_dGEp* dGEp)**2 + (dsig_dGMp*dGMp)**2
    dsig = np.sqrt(dsig2)
    
    return sig, dsig



#------User Input:
#Q2 [GeV2]
#final e- momentum (kf) [GeV]
#final e- angle (th_e) [deg]

Eb = 2.070
th_e = np.linspace(20,16,100)
kf = np.linspace(1.8,1.9,100)
Q2 = 4*Eb*kf* pow(np.sin( th_e/2. * dtr), 2)
#th_e = 2.*asin( np.sqrt(Q2/(4.*Eb*kf))) / dtr # deg
#kf = Q2 / ( 4. * Eb * pow(np.sin( th_e/2. * dtr), 2) )

#Q2 = np.linspace(0.36, 0.42, 100) # generate equally spaced Q2 points

#th_e = np.linspace(17.4, 19.2, 100)                          # select ~ central e- angle based on distribution (this is an approximation)
#kf = np.linspace(1.82, 1.9, 100)  # calculate e- scattering angle [GeV]
#kf = Q2 / ( 4. * Eb * pow(np.sin( th_e/2. * dtr), 2) )

nu = Eb-kf  # energy transfer [GeV]

#Get the numerical values of the fofa uncertainties from JRA parametrization interpolation (see interpolate_ff.py)
GEp_GD = get_GEp_GD(Q2)
dGEp_GD = get_dGEp_GD(Q2)
GMp_muGD = get_GMp_muGD(Q2)
dGMp_muGD = get_dGMp_muGD(Q2)

Xsec, Xsec_err = sig(kf, th_e, Q2, GEp_GD, GMp_muGD, dGEp_GD, dGMp_muGD)
Xsec_rel_err = (Xsec_err / Xsec)*100.  #relatice cross section error in percent

print('H(e,e)p Elastic Cross Section')
print('Q2 = ', Q2, ' [GeV^2]')
print('e- Beam Energy (Eb) = ', Eb, ' [GeV]')
print('e- momentum (kf) = ', kf, ' [GeV/c]')
print('e- angle (th_e) = ', th_e, ' [deg]')
print('e- Energy Transfer (nu) = ', nu, ' [GeV]')
print('Xsec (ub/sr) = ',Xsec,' +/- ',Xsec_err, ' [microbarn/sr]')
print('Xsec_rel_err=', Xsec_rel_err, ' [%]' )

plt.errorbar(nu, Xsec, Xsec_err, marker='', linestyle='--', color='r')
plt.show()
