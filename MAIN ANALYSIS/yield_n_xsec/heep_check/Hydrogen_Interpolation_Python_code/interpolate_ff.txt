import LT.box as B
import numpy as np
from LT.datafile import dfile
from scipy.interpolate import interp1d


#open proton form factor data with errors (see supplemental material for definitions of headers)
f = dfile('proton_lookup_copy.data')

Q2 = f['Q2']   #GeV2
GEp_GD = f['GEp_GD']
dGEp_GD = f['dGEp_GD']
GMp_muGD = f['GMp_muGD']
dGMp_muGD = f['dGMp_muGD']

#Interpolate Form factors
f_GEp_GD = interp1d(Q2, GEp_GD, fill_value='extrapolate') 
f_dGEp_GD = interp1d(Q2, dGEp_GD, fill_value='extrapolate') 

f_GMp_muGD = interp1d(Q2, GMp_muGD, fill_value='extrapolate') 
f_dGMp_muGD = interp1d(Q2, dGMp_muGD, fill_value='extrapolate') 


#define function to get form factors
def get_GEp_GD(Q2):
    result = f_GEp_GD(Q2)
    return result

def get_dGEp_GD(Q2):
    result = f_dGEp_GD(Q2)
    return result


def get_GMp_muGD(Q2):
    result = f_GMp_muGD(Q2)
    return result

def get_dGMp_muGD(Q2):
    result = f_dGMp_muGD(Q2)
    return result
