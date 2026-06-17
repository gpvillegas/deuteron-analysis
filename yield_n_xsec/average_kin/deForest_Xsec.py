#Code to calculate the deForest Cross Sections from the proton GMp, GEp form
#factors

import LT.box as B
import numpy as np
    
#Define Physical Constnats
hbarc = 197.327053       #MeV * fm
alpha = 1./137.0359895    #structure constant
dtr = np.pi / 180.

#Masses in MeV
MP = 938.272
MN = 939.566
MD = 1875.61
me = 0.51099; 

#define the proton magnetic form factor (See physics_proton.f)
#Assumes Q2 is in MeV2
def GMp(Q2, param=''):

    Q2 = Q2*1e-6    #convert from MeV2 to GeV2
    #Four momentum transfer to various powers (units GeV^{power})
    Q = np.sqrt(Q2)
    mu_p = 2.793      #proton magnetic moment

    #Peter Bosted Parametrization
    if param=='bosted':
        Q3 = Q**3
        Q4 = Q**4
        Q5 = Q**5    
        denom =  1. + 0.35*Q + 2.44*Q2 + 0.5*Q3 + 1.04*Q4 + 0.34*Q5
        GMp = mu_p / denom

    #John R. Arrington Parametrization
    if param=='JRA':
        Q22 = Q2**2
        Q23 = Q2**3 
        Q24 = Q2**4 	
        Q25 = Q2**5  
        Q26 = Q2**6
        GMp = mu_p/(1. + Q2 * (3.19) + Q22 * (1.355) + Q23 * (0.151) + Q24 * (-0.114E-01) + Q25 * (0.533E-03) + Q26 * (-0.900E-05) )

    return GMp

#define the proton electric form factor
#Assumes Q2 input is in MeV2
def GEp(Q2, param=''):
    #print ('Calculating Electric Form factor')

    Q2 = Q2*1e-6    #convert from MeV2 to GeV2

    #Four momentum transfer to various powers (units GeV^{power})
    Q = np.sqrt(Q2)
    
    #Peter Bosted Parametrization
    if param=='bosted':
        Q3 = Q**3
        Q4 = Q**4
        Q5 = Q**5
        
        denom =  1. + 0.62*Q + 0.68*Q2 + 2.8*Q3 + 0.83*Q4
        GEp = 1. / denom

    #John R. Arrington Parametrization
    if param=='JRA':
        Q22 = Q2**2  
        Q23 = Q2**3  
        Q24 = Q2**4 
        Q25 = Q2**5  
        Q26 = Q2**6 
        GEp=1./(1. + Q2 * (3.226) + Q22 * (1.508) + Q23 * (-0.3773) + Q24 * (0.611) + Q25 * (-0.1853) + Q26 * (0.1596E-01)  )

    return GEp


#Mott Cross Section
def sigMott(kf, th_e, Q2):

    #Units: kf ~ Ef [MeV]   electron momentum  
    #th_e[deg] electron angle  
    #Q2[MeV2]   4-momentum transfer
    if Q2<=0.:
        the_e = -1.
        sig = -1.
        sigMott = -1.
    else:
        th_e = th_e * dtr    # convert degree to radians
        sig = (2.*alpha*hbarc*kf*np.cos(th_e/2.)/Q2 )**2   #MeV^2 * fm^2 *MeV^2/MeV^4--> fm2
        sigMott = sig * 1e4    #convert fm^2 -> microbarn

    return sigMott


#deForest eN off-shell cross section (Refer to deForest 1983 paper)
#Input units: energy: MeV  angle: deg
#Use the average kinematics as inpu
def deForest(Ef, Q2, q, Pf, Pm, th_e, th_p, c_phi, thpq, sig_Mott, GE_p, GM_p):

    #Input Parameters: Q2, |q_lab|(3-vector q magnitude), Pf, Ef, theta_e, cos_phi, 

    #Ef: final e- energy, Pf: final proton energy, 
    #Er: recoil energy (neutron)
    #d5sig = k * sig_eN * S(Pm)
    # sig_eN = sigMott * { Vc * Wc + Vt * Wt + Vi * Wi  + Vs * Ws }, where the V's are the coefficients and W's are the respose functions

    #NOTE: The formulas are copied DIRECTLY from deForest (1983). For Mott scattering, deForest
    #did NOT put an (hbarc)^2 units, so I added it explicitly in sigMott(), in units of ub / sr
 
    #The 4-momentum transfer squared q_mu*q^mu, we define as -Q2
    q2mu = -Q2    #MeV^2 
    q4mu = q2mu**2

    #Convert Degrees to Radians
    th_e = th_e * dtr   #electron angle 
    gamma = thpq * dtr   #in-plane angle between struck proton in q-vector

    Epf = np.sqrt(MP*MP + Pf*Pf)  #final proton energy
    Er = np.sqrt(MN*MN + Pm*Pm)  #final neutron (recoil) energy

    #For a 6-Fold Differential XSec (d6sig/(dE'*dom *dOme* dOmp))   --See Hari's THesis
    Kfact = Pf*Epf     #kinematic factor, final proton momentum, energy  [MeV^2] (THis should also have a 1/(2*pi)^3)
 
    #For a 5-Fold Differential XSec (d5sig/(dom *dOme* dOmp))  : Missing Energy Integrated  (sig)deutpwia.f --Werner's code)
    #Kfact = Pf * Epf *MN / (8.*np.pi**3 * MD)  #from Hari's thesis
    
    pipf = 0.5 * (q**2 - (Pf**2 + Pm**2))
    rec = 1 -  (Epf/Er) * (pipf/Pf**2)  
    f_rec = 1./rec

    #Using definitions from William P. Ford, Sabine and J.W. Van Orden  (factor too small. check later on.)
    #om = np.sqrt(q*q - Q2)  #energy transfer
    #Kfact = MP*MN*Pf/(8.* np.pi**3 * MD)
    #rec = np.abs(1 + (om*Pf - Epf*q*np.cos(th_p)/(MP*Pf)))
    #f_rec=1. / rec

    # print('Kfact=',Kfact)
    # print('f_rec=',f_rec)
    #print('K*f_rec=', Kfact*f_rec)

    #Define the Response Functions Coefficients
    Vc = q4mu / q**4
    Vt = -q2mu / (2.*q**2) + np.tan(th_e/2.)**2
    Vi = -q2mu/q**2 * np.sqrt(-q2mu/q**2 + np.tan(th_e/2.)**2) * c_phi
    Vs = -q2mu/q**2 * c_phi**2 + np.tan(th_e/2.)**2

    #Define bar quantities
    Ebar_f = np.sqrt(Pm*Pm + MP*MP)   #struck proton energy
    om_bar = Epf - Ebar_f              #bar energy transfer (deForest calls omega bar)
    q2mu_bar = q**2 - om_bar**2
    EbarE = Ebar_f*Epf

    #Defing the Sachs Form Factors 
    #From: GEp = F1 - tau * F2,  GMp = F1 + F2,  where Tau = Q2 / 4Mp**2
    tau = Q2 / (4.*MP**2)
    F1 = (GE_p + tau*GM_p) / (1. + tau)
    kF2 = (GM_p - GE_p) /  (1. + tau)    

    sumFF1 = (F1 + kF2)**2
    sumFF2 = F1**2 + (q2mu_bar/(4.*MP**2))*kF2**2   

    Wc = 1./(4.*EbarE) * ( (Ebar_f + Epf)**2 * (F1**2 + (q2mu_bar/(4.*MP**2))*kF2**2) - q**2*(F1 + kF2)**2 )
    Wt = q2mu_bar/(2.*EbarE) * (F1 + kF2)**2
    Ws = Pf**2 * np.sin(gamma)**2 / (EbarE) * (F1**2 + (q2mu_bar/(4.*MP**2))*kF2**2)
    Wi = -Pf * np.sin(gamma) * (Ebar_f + Epf) / (EbarE) * (F1**2 + q2mu_bar/(4.*MP**2)*kF2**2)

    #Calculate eN offshell cross section (sig_cc1)
    sig_eN = sig_Mott * ( Vc*Wc + Vt*Wt + Vs*Ws + Vi*Wi ) 

    #DEBUG
    '''
    print('sig_Mott=',sig_Mott) 
    print(':Vc=',Vc, ':Wc=',Wc)
    print(':Vt=',Vt, ':Wt=',Wt)
    print(':Vs=',Vs, ':Ws=',Ws)
    print(':Vi=',Vi, ':Wi=',Wi)
    print('sumFF1=',sumFF1, ': sumFF2=', sumFF2)
    '''

    #Calculate the deForest Cross Section
    #d^5sigma/(dOm_e*dOm_p*dE') = K * sig_eN * S(Pm):  [ub/(sr^2*MeV)] = [MeV^2 * ub/sr^2 * 1/MeV^3]
    deForest = Kfact * f_rec * sig_eN      #this is actually: d^5sigma/(dOm_e*dOm_p*dE') / S(Pm) = K * sig_eN, so the units are: (ub * MeV^2) / sr^2   
    
    return Kfact, f_rec, sig_eN, deForest


def main():
    print('Calling main() . . .')

    kf =  8486.2257205347632
    Pf =  2866.5450933307575
    Pm =  218.97375008520876
    q = 3006.0732577503345
    th_e = 0.22681158350324249 / dtr
    Q2 =  4604071.9136472372
    GE_p = GEp(Q2)
    GM_p = GMp(Q2)
    sig_Mott = sigMott(kf, th_e, Q2)
    c_phi = 0.51778477668111345 
    sin_gamma = 5.7467453222693969E-002
    thpq = np.arcsin(sin_gamma) / dtr
    de_Forest = deForest(Q2, q, Pf, Pm, th_e, c_phi, thpq, sig_Mott, GE_p, GM_p)

    #print('GEp=',GE_p,':GMp=',GM_p,':sigMott=',sig_Mott,':deForest=',de_Forest)

if __name__ == "__main__":
    main()
