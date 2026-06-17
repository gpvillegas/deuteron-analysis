#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:36:59 2026

@author: gvill

Module for load_data_init:
    function: loads data or SIMC files based on selection.
    DATA can be selected by run, kin_study, or setting,
    SIMC files can be selected by theory model and kin_study or 
    setting.
"""
import data_init as D
import root_util as R

import numpy as np
import database_operations as db
import vector as vec
from scipy.spatial.transform import Rotation as Rot
from typing import Any, Union

#%%
"""
CLASS load_data_init:
    Load selected setting data_init objects with selected branches.
    
    setting: (str) options: 'pm_[#]', 'delta_scan_[#]', 'lumi_[tar]_[#]'
        where [#] -> specific setting number, [tar] -> target
    data_dir: (str) directory where data is, leave empty for default location
    simc_dir: (str) directory where simc files are
    simc_models: (list[(str),...]) list of models to be loaded. 
        options: 'JML', 'Paris', 'V18', 'CD-Bonn'
        if empty list no SIMC files will be loaded
    load_comm_data: (bool) set True to load commissioning data
    load_data: (bool) set True to load data
"""

class load_data_init:
    def __init__(self,setting='pm_120',data_dir='',simc_dir='',
                 simc_models=['JML'],load_comm_data=True,load_data=True):
        self.setting = setting
        self.data_dir = data_dir
        self.simc_dir = simc_dir
        self.comm_data_dir = '/home/gvill/deuteron/deut_offline_replay/'+\
                                'code_development/commissioning_data/'
        self.simc_models = simc_models
        self.load_data = load_data
        self.load_cd = load_comm_data
        self.load_JML = False
        self.load_PAR = False
        self.load_V18 = False
        self.load_CDB = False
        # select appropriate HCANA variables 
        kin_var = ['H.kin.secondary.th_bq','H.kin.secondary.pmiss',
                   'H.kin.secondary.emiss_nuc','P.kin.primary.Q2']

        sieve_var = ['H.extcor.xsieve', 'H.extcor.ysieve',
                          'P.extcor.xsieve', 'P.extcor.ysieve']

        coinTime_var = ['CTime.epCoinTime_ROC2']

        acc_var = ['H.gtr.dp','P.gtr.dp','H.gtr.th','H.gtr.ph',
                   'P.gtr.th','P.gtr.ph']

        calPID_var = ['P.cal.etottracknorm']

        z_react_var = ['H.react.z','P.react.z']
        
        EDTM_var = ['T.coin.pEDTM_tdcTimeRaw']

        T_sel = kin_var + coinTime_var + acc_var + calPID_var + sieve_var +\
            z_react_var + EDTM_var

        TSP_sel = ['evNumber','P.BCM1.scalerCurrent']

        self.br_sel = {'T':T_sel,'TSP':TSP_sel}
        self.t_sel = ['T','TSP']
        
        self.br_sel_sim = ['theta_rq','Pm','Em','Q2','e_delta','h_delta',
                           'Weight','Normfac','e_xptar','e_yptar','h_xptar',
                           'h_yptar','tar_x','h_zv','e_zv','h_ytar','e_ytar',
                           'sig','sigcc','xcoll','ycoll','Jacobian_corr',
                           'probabs','Ein_v','pf_v','q_lab_v','Q2_v','nu_v',
                           'pm_v','e_xptar_v','e_yptar_v','h_xptar_v',
                           'h_yptar_v','Ef_v','sig']
        if self.load_data:
            self.__load_data()
        
        if self.simc_models:
            self.__load_simc()
            
        if self.load_cd:
            self.__load_comm_data()
        
        print('Finished loading all data.')
    
    def __load_data(self):
        if self.setting:
            print('**Loading data...')
            self.load_data = True
            data = D.DATA_INIT(data_type='deut23_data',
                                     kin_study='deep',
                                     setting=self.setting,
                                     select_branches=self.br_sel,
                                     select_trees=self.t_sel,
                                     ROOTfiles_path=self.data_dir)
            self.data = data
        else:
            print('No data setting selected. Not loading data.')                
                
    def __load_simc(self):
        print('** Loading SIMC files...')
        for sm in self.simc_models:
            match sm:
                case 'JML':
                    self.load_JML = True
                    JML_FSIrad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='jmlfsi_rad')
                    JML_FSInorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='jmlfsi_norad')
                    JML_PWIAnorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='jmlpwia_norad')
                    JML_PWIArad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='jmlpwia_rad')
        
                    self.JML_FSIrad = JML_FSIrad
                    self.JML_FSInorad = JML_FSInorad
                    self.JML_PWIAnorad = JML_PWIAnorad
                    self.JML_PWIArad = JML_PWIArad
                    
                case 'Paris':
                    self.load_PAR = True
                    PAR_FSIrad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='Paris_fsi_rad')
                    PAR_FSInorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='Paris_fsi_norad')
                    PAR_PWIAnorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='Paris_pwia_norad')
                    PAR_PWIArad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='Paris_pwia_rad')
                    self.PAR_FSIrad = PAR_FSIrad
                    self.PAR_FSInorad = PAR_FSInorad
                    self.PAR_PWIAnorad = PAR_PWIAnorad
                    self.PAR_PWIArad = PAR_PWIArad
                case 'V18':
                    self.load_V18 = True
                    V18_FSIrad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='V18_fsi_rad')
                    V18_FSInorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='V18_fsi_norad')
                    V18_PWIAnorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='V18_pwia_norad')
                    V18_PWIArad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='V18_pwia_rad')
                    self.V18_FSIrad = V18_FSIrad
                    self.V18_FSInorad = V18_FSInorad
                    self.V18_PWIAnorad = V18_PWIAnorad
                    self.V18_PWIArad = V18_PWIArad
                case 'CD-Bonn':
                    self.load_CDB = True
                    CDB_FSIrad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='CD-Bonn_fsi_rad')
                    CDB_FSInorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='CD-Bonn_fsi_norad')
                    CDB_PWIAnorad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='CD-Bonn_pwia_norad')
                    CDB_PWIArad = D.DATA_INIT(data_type='SIMC',setting = self.setting,
                                            select_branches={'SNT':self.br_sel_sim},
                                            SIMC_ROOTfiles_path=self.simc_dir+self.setting+'/',
                                                  simc_type='CD-Bonn_pwia_rad')
                    self.CDB_FSIrad = CDB_FSIrad
                    self.CDB_FSInorad = CDB_FSInorad
                    self.CDB_PWIAnorad = CDB_PWIAnorad
                    self.CDB_PWIArad = CDB_PWIArad
                case _:
                    print('Model not recognized or not selected. No SIMC files loaded.')
        
    def __load_comm_data(self):
        match self.setting:
            case 'pm_120':
                print('**Loading commissioning data...')
                pm80_filename = 'pm80/Xsec_pm80_lagetfsi_dataset1.root'
                pm80_data = R.root_histos(self.comm_data_dir+pm80_filename)
                pm80_data.load_histos()
                
                self.pm80_data = pm80_data
            case 'pm_580':
                print('**Loading commissioning data...')
                pm580_filename1 = 'pm580/Xsec_pm580_lagetfsi_dataset1.root'
                
                pm580_comm_data1 = R.root_histos(self.comm_data_dir+pm580_filename1)
                pm580_comm_data1.load_histos()
                
                pm580_filename2 = 'pm580/Xsec_pm580_lagetfsi_dataset2.root'
                
                pm580_comm_data2 = R.root_histos(self.comm_data_dir+pm580_filename2)
                pm580_comm_data2.load_histos()
                
                self.pm580_comm_data1 = pm580_comm_data1
                self.pm580_comm_data2 = pm580_comm_data2
            case 'pm_800':
                print('**Loading commissioning data...')
                pm580_filename1 = 'pm580/Xsec_pm580_lagetfsi_dataset1.root'
                
                pm580_comm_data1 = R.root_histos(self.comm_data_dir+pm580_filename1)
                pm580_comm_data1.load_histos()
                
                pm580_filename2 = 'pm580/Xsec_pm580_lagetfsi_dataset2.root'
                
                pm580_comm_data2 = R.root_histos(self.comm_data_dir+pm580_filename2)
                pm580_comm_data2.load_histos()
                
                pm750_filename1 = 'pm750/Xsec_pm750_lagetfsi_dataset1.root'
                
                pm750_data1 = R.root_histos(self.comm_data_dir+pm750_filename1)
                pm750_data1.load_histos()
                
                pm750_filename2 = 'pm750/Xsec_pm750_lagetfsi_dataset2.root'
                
                pm750_data2 = R.root_histos(self.comm_data_dir+pm750_filename2)
                pm750_data2.load_histos()
                
                pm750_filename3 = 'pm750/Xsec_pm750_lagetfsi_dataset3.root'
                
                pm750_data3 = R.root_histos(self.comm_data_dir+pm750_filename3)
                pm750_data3.load_histos()
                
                self.pm580_comm_data1 = pm580_comm_data1
                self.pm580_comm_data2 = pm580_comm_data2
                self.pm750_data1 = pm750_data1
                self.pm750_data2 = pm750_data2
                self.pm750_data3 = pm750_data3
            case 'pm_900':
                print('**Loading commissioning data...')
                pm580_filename1 = 'pm580/Xsec_pm580_lagetfsi_dataset1.root'
                
                pm580_comm_data1 = R.root_histos(self.comm_data_dir+pm580_filename1)
                pm580_comm_data1.load_histos()
                
                pm580_filename2 = 'pm580/Xsec_pm580_lagetfsi_dataset2.root'
                
                pm580_comm_data2 = R.root_histos(self.comm_data_dir+pm580_filename2)
                pm580_comm_data2.load_histos()
                
                pm750_filename1 = 'pm750/Xsec_pm750_lagetfsi_dataset1.root'
                
                pm750_data1 = R.root_histos(self.comm_data_dir+pm750_filename1)
                pm750_data1.load_histos()
                
                pm750_filename2 = 'pm750/Xsec_pm750_lagetfsi_dataset2.root'
                
                pm750_data2 = R.root_histos(self.comm_data_dir+pm750_filename2)
                pm750_data2.load_histos()
                
                pm750_filename3 = 'pm750/Xsec_pm750_lagetfsi_dataset3.root'
                
                pm750_data3 = R.root_histos(self.comm_data_dir+pm750_filename3)
                pm750_data3.load_histos()
                
                self.pm580_comm_data1 = pm580_comm_data1
                self.pm580_comm_data2 = pm580_comm_data2
                self.pm750_data1 = pm750_data1
                self.pm750_data2 = pm750_data2
                self.pm750_data3 = pm750_data3
            case _:
                print('Commisioning data not loaded.')

#%%
MP = 0.938272     #proton mass
MN = 0.939566     #neutron mass
MD = 1.87561      #deuteron mass
me = 0.00051099   #electron mass
dtr = np.pi/180   #degree-to-radian

# function to convert geographical angles to spherical, in radian
def GeoToSph(thetaGeo, phiGeo):
    ct = np.cos(thetaGeo)
    cp = np.cos(phiGeo)
    tmp = ct*cp
    thetaSph = np.acos(tmp)
    tmp = np.sqrt(1. - tmp*tmp)
    
    # change angle to 0 if really small, otherwise calculate phi
    phiSph = np.where(np.abs(tmp) < 1e-6, 0.0, np.acos(np.sqrt(1.0-ct*ct)*cp/tmp))
    
    # check if thetaGeo is from 0 to pi
    condition = (thetaGeo/(2*np.pi) - np.floor(thetaGeo/(2*np.pi))) > 0.5
    phiSph = np.where(condition, np.pi-phiSph, phiSph)
    
    # check if phiGeo is from -pi to pi
    condition = (phiGeo/(2*np.pi) - np.floor(phiGeo/(2*np.pi))) > 0.5
    phiSph = np.where(condition, -phiSph, phiSph)
    return (thetaSph, phiSph)

# function to set up rotation matrix to lab coordinates based 
# on spectrometer coords
def SetCentralAngles(th, ph):
    th_geo = th*dtr
    ph_geo = ph*dtr
    th_sph, ph_sph = GeoToSph(th_geo,ph_geo)
    
    st = np.sin(th_sph)
    ct = np.cos(th_sph)
    sp = np.sin(ph_sph)
    cp = np.cos(ph_sph)
    norm = np.sqrt(ct*ct + st*st*cp*cp)
    nx = np.array([st*st*sp*cp/norm, -norm, st*ct*sp/norm])
    ny = np.array([ct/norm, 0.0, -st*cp/norm])
    nz = np.array([st*cp, st*sp, ct])
    
    rot_matrix = {'xx':nx[0],'xy':ny[0],'xz':nz[0],
                  'yx':nx[1],'yy':ny[1],'yz':nz[1],
                  'zx':nx[2],'zy':ny[2],'zz':nz[2]}
    return rot_matrix

# function to transform a momentum vector into lab coords  
def TransportToLab(p, xptar, yptar, rot):
    v = vec.array({'px':xptar, 'py':yptar, 'pz':np.ones(xptar.size)})
    v *= p/np.sqrt(1.0+xptar*xptar+yptar*yptar)
    pvec = v.transform3D(rot)
    return pvec

# function to transform from lab to q
def _to_numpy_vectors(obj: Any) -> np.ndarray:
    """
    Convert a 3D vector or vector array to a NumPy array of shape (..., 3).
    Supports:
      - numpy arrays
      - lists/tuples
      - vector library objects with .to_numpy()
      - vector library objects with .x, .y, .z attributes
    """
    if hasattr(obj, "to_numpy"):
        arr = np.asarray(obj.to_numpy(), dtype=float)
    elif hasattr(obj, "x") and hasattr(obj, "y") and hasattr(obj, "z"):
        x = np.asarray(obj.x, dtype=float)
        y = np.asarray(obj.y, dtype=float)
        z = np.asarray(obj.z, dtype=float)
        arr = np.stack([x, y, z], axis=-1)
    else:
        arr = np.asarray(obj, dtype=float)

    if arr.ndim == 1 and arr.shape[0] == 3:
        return arr.reshape(1, 3)
    if arr.ndim == 2 and arr.shape[1] == 3:
        return arr
    raise ValueError("Input must be a 3-vector or an array of 3-vectors")

def align_vectors(v1: Any, v2: Any) -> np.ndarray:
    v1 = _to_numpy_vectors(v1).reshape(3)
    v2 = _to_numpy_vectors(v2).reshape(3)

    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)

    if np.allclose(v1_norm, v2_norm):
        return np.eye(3)

    if np.allclose(v1_norm, -v2_norm):
        perp = np.array([1.0, 0.0, 0.0]) if abs(v1_norm[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(v1_norm, perp)
        axis /= np.linalg.norm(axis)
        return _rotation_matrix_from_axis_angle(axis, np.pi)

    axis = np.cross(v1_norm, v2_norm)
    axis /= np.linalg.norm(axis)
    angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0))
    return _rotation_matrix_from_axis_angle(axis, angle)

def align_vector_pairs(v1_array: Any, v2_array: Any) -> np.ndarray:
    """
    This function uses Rodriguez Formula 
    (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)
    to align two vectors. I use this to transform momentum vectors from
    the lab to q frame.
    """
    v1 = _to_numpy_vectors(v1_array)
    v2 = _to_numpy_vectors(v2_array)

    if v1.shape != v2.shape:
        raise ValueError("Both inputs must have the same shape")
    if v1.ndim != 2 or v1.shape[1] != 3:
        raise ValueError("Inputs must be arrays of 3-vectors")

    n = v1.shape[0]
    v1_norm = v1 / np.linalg.norm(v1, axis=1, keepdims=True)
    v2_norm = v2 / np.linalg.norm(v2, axis=1, keepdims=True)

    axes = np.cross(v1_norm, v2_norm)
    dots = np.sum(v1_norm * v2_norm, axis=1)
    angles = np.arccos(np.clip(dots, -1.0, 1.0))

    axis_norms = np.linalg.norm(axes, axis=1, keepdims=True)
    axis_norms[axis_norms == 0] = 1.0
    axes_normed = axes / axis_norms

    K = np.zeros((n, 3, 3), dtype=float)
    K[:, 0, 1] = -axes_normed[:, 2]
    K[:, 0, 2] = axes_normed[:, 1]
    K[:, 1, 0] = axes_normed[:, 2]
    K[:, 1, 2] = -axes_normed[:, 0]
    K[:, 2, 0] = -axes_normed[:, 1]
    K[:, 2, 1] = axes_normed[:, 0]

    K2 = np.matmul(K, K)
    I = np.eye(3)
    rotation_matrix = (
        I
        + np.sin(angles)[:, None, None] * K
        + (1 - np.cos(angles))[:, None, None] * K2
    )

    parallel = axis_norms.flatten() < 1e-12
    for i in np.where(parallel)[0]:
        if np.allclose(dots[i], 1.0):
            rotation_matrix[i] = np.eye(3)
        else:
            perp = np.array([1.0, 0.0, 0.0]) if abs(v1_norm[i, 0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            axis = np.cross(v1_norm[i], perp)
            axis /= np.linalg.norm(axis)
            rotation_matrix[i] = _rotation_matrix_from_axis_angle(axis, np.pi)
    
    xx = rotation_matrix[:,0,0]
    xy = rotation_matrix[:,0,1]
    xz = rotation_matrix[:,0,2]
    yx = rotation_matrix[:,1,0]
    yy = rotation_matrix[:,1,1]
    yz = rotation_matrix[:,1,2]
    zx = rotation_matrix[:,2,0]
    zy = rotation_matrix[:,2,1]
    zz = rotation_matrix[:,2,2]
    
    rot_matrix = {'xx':xx,'xy':xy,'xz':xz,
                  'yx':yx,'yy':yy,'yz':yz,
                  'zx':zx,'zy':zy,'zz':zz}
    return rot_matrix

def _rotation_matrix_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0]
    ])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

def LabtoQ(qvec,kvec):   
    xx = []
    xy = []
    xz = []
    yx = []
    yy = []
    yz = []
    zx = []
    zy = []
    zz = []
    
    for i in range(qvec.size):
        q = np.array([qvec[i].x,qvec[i].y,qvec[i].z])
        k = np.array([kvec[i].x,kvec[i].y,kvec[i].z])
        
        rot, rssd = Rot.align_vectors(k, q)
        r = rot.as_matrix()
        
        xx.append(r[0][0])
        xy.append(r[0][1])
        xz.append(r[0][2])
        yx.append(r[1][0])
        yy.append(r[1][1])
        yz.append(r[1][2])
        zx.append(r[2][0])
        zy.append(r[2][1])
        zz.append(r[2][2])
        
    rot_matrix = {'xx':xx,'xy':xy,'xz':xz,
                  'yx':yx,'yy':yy,'yz':yz,
                  'zx':zx,'zy':zy,'zz':zz}    
    return rot_matrix
    

def calc_vertex_kin(data_obj):
    if type(data_obj.many) is list:
        first = True
        for m in data_obj.many:
            # get SIMC variables from data_obj and convert MeV -> GeV
            Ein_v = data_obj.Branches[m]['Ein_v']/1e3
            Ef_v = data_obj.Branches[m]['Ef_v']/1e3
            Pf_v = data_obj.Branches[m]['pf_v']/1e3
            Q2_v = data_obj.Branches[m]['Q2_v']/1e6
            omega_v = data_obj.Branches[m]['nu_v']/1e3
            Pm_v = data_obj.Branches[m]['pm_v']/1e3
            q_v = data_obj.Branches[m]['q_lab_v']/1e3
            e_xptar_v = data_obj.Branches[m]['e_xptar_v']
            e_yptar_v = data_obj.Branches[m]['e_yptar_v']
            h_xptar_v = data_obj.Branches[m]['h_xptar_v']
            h_yptar_v = data_obj.Branches[m]['h_yptar_v']
            
            n = Ein_v.size  # number of events
            
            ki_v = np.sqrt(Ein_v*Ein_v - me*me)     # initial e momentum
            kf_v = np.sqrt(Ef_v*Ef_v - me*me)       # final e momentum
            xbj_v = Q2_v/(2*MP*omega_v)             # Bjorken-x
            
            # get e and p central angles from db
            e_th, h_th = db.retrieve('deuteron_db.db', 'SHMS_Angle, HMS_Angle', 'RUN_LIST_UPDATED',
                                     where= f"kin_study=\'{data_obj.kin_study[m]}\'"+\
                                         f" AND setting=\'{data_obj.setting[m]}\'",distinct=True)[0]
            e_rot_matrix = SetCentralAngles(e_th, 0.) # get rotation matrix from e-spec coords to lab
            kf_lab = TransportToLab(kf_v, e_xptar_v, e_yptar_v, e_rot_matrix) # e 3-momentum vector in lab frame
            h_rot_matrix = SetCentralAngles(-h_th, 0.) # get rotation matrix from h-spec coords to lab
            pf_lab = TransportToLab(Pf_v, h_xptar_v, h_yptar_v, h_rot_matrix) # h 3-momentum vector in lab frame
                
            #create 4-vectors
            # initial e 4-momentum vector
            P0_v = vec.array({'px':np.zeros(n),'py':np.zeros(n),'pz':ki_v,'m':np.ones(n)*me})
            # final e 4-momentum vector
            P1_v = vec.array({'px':kf_lab.px,'py':kf_lab.py,'pz':kf_lab.pz,'m':np.ones(n)*me})
            # initial target (assumed at rest)
            A_v = vec.array({'px':np.zeros(n),'py':np.zeros(n),'pz':np.zeros(n),'m':np.ones(n)*MD})
            Q_v = P0_v - P1_v   # q 4-vector
            A1_v = A_v - Q_v    # final target 4-momentum based on conservation
            # final detected system 4-momentum vector
            X_v = vec.array({'px':pf_lab.px,'py':pf_lab.py,'pz':pf_lab.pz,'m':np.ones(n)*MP})
            # undetected system 4-momentum
            B_v = A1_v - X_v
            
            # Missing momentum components (lab)
            Pmx_lab_v = B_v.px
            Pmy_lab_v = B_v.py
            Pmz_lab_v = B_v.pz
            
            # Check if calcualted Pm_v is the same as SIMC pm_v
            Pm_v_calc = np.sqrt(Pmx_lab_v*Pmx_lab_v + Pmy_lab_v*Pmy_lab_v + Pmz_lab_v*Pmz_lab_v)
            
            # e and p in-plane angles
            the_v = kf_lab.theta
            thp_v = pf_lab.theta
            
            # Now rotate recoil system from lab to q
            thq_v = Q_v.theta  # in-plane angle of q in lab
            # phq_v = Q_v.phi    # out-of-plane angle of q in lab
            
            xq_v = X_v.rotateZ(thq_v)   # detected system rotated to q frame
            bq_v = B_v.rotateZ(thq_v)   # undetected system rotated to q frame
            
            # angles of rotated recoil system
            thpq_v = xq_v.theta
            phpq_v = xq_v.phi
            thnq_v = bq_v.theta
            phnq_v = bq_v.phi
            
            #check th_nq
            
            # missing momentum in q frame
            Pmq_v = -bq_v
            Pmx_q_v = Pmq_v.px
            Pmy_q_v = Pmq_v.py
            Pmz_q_v = Pmq_v.pz
            
            # Compare q_v_calc with SIMC q_lab_v
            # Compare Pm_v_calc with SIMC pm_v
            kin = {'xbj_v':xbj_v,
                   'q_v_calc':Q_v.p,
                   'pm_v_calc':Pm_v_calc,
                   'the_v':the_v/dtr,
                   'thp_v':thp_v/dtr,
                   'th_pq_v':thpq_v/dtr,
                   'ph_pq_v':phpq_v/dtr,
                   'th_nq_v':thnq_v/dtr,
                   'ph_nq_v':phnq_v/dtr,
                   'Pm_q_v':Pmq_v,
                   'Pmx_q_v':Pmx_q_v,
                   'Pmy_q_v':Pmy_q_v,
                   'Pmz_q_v':Pmz_q_v,
                   'kf_v':kf_v,
                   'thq_v':thq_v/dtr,
                   'pm_v_GeV':Pm_v,
                   'Q2_v_GeV2':Q2_v,
                   'nu_v_GeV':omega_v,
                   'q_lab_v_GeV':q_v,
                   'pf_v_GeV':Pf_v,
                   'Ein_v_GeV':Ein_v
                   }
            # update data_init object with new branches
            data_obj.Branches[m].update(kin)
            if first:
                print(f'data_init object at {m} updated with the following branches:',
                      kin.keys())
                first = False
            else:
                print(f'data_init object at {m} updated')
    else:
        # get SIMC variables from data_obj and convert MeV -> GeV
        Ein_v = data_obj.Branches['Ein_v']/1e3
        Ef_v = data_obj.Branches['Ef_v']/1e3
        Pf_v = data_obj.Branches['pf_v']/1e3
        Q2_v = data_obj.Branches['Q2_v']/1e6
        omega_v = data_obj.Branches['nu_v']/1e3
        Pm_v = data_obj.Branches['pm_v']/1e3
        q_v = data_obj.Branches['q_lab_v']/1e3
        e_xptar_v = data_obj.Branches['e_xptar_v']
        e_yptar_v = data_obj.Branches['e_yptar_v']
        h_xptar_v = data_obj.Branches['h_xptar_v']
        h_yptar_v = data_obj.Branches['h_yptar_v']
        
        n = Ein_v.size  # number of events
        
        ki_v = np.sqrt(Ein_v*Ein_v - me*me)     # initial e momentum
        kf_v = np.sqrt(Ef_v*Ef_v - me*me)       # final e momentum
        xbj_v = Q2_v/(2*MP*omega_v)             # Bjorken-x
        
        # get e and p central angles from db
        e_th, h_th = db.retrieve('deuteron_db.db', 'SHMS_Angle, HMS_Angle', 'RUN_LIST_UPDATED',
                                 where= f"kin_study=\'{data_obj.kin_study}\'"+\
                                     f" AND setting=\'{data_obj.setting}\'",distinct=True)[0]
        e_rot_matrix = SetCentralAngles(e_th, 0.) # get rotation matrix from e-spec coords to lab
        kf_lab = TransportToLab(kf_v, e_xptar_v, e_yptar_v, e_rot_matrix) # e 3-momentum vector in lab frame
        h_rot_matrix = SetCentralAngles(-h_th, 0.) # get rotation matrix from h-spec coords to lab
                                                    # set HMS angle to negative
        pf_lab = TransportToLab(Pf_v, h_xptar_v, h_yptar_v, h_rot_matrix) # h 3-momentum vector in lab frame  
        #create 4-vectors
        # initial e 4-momentum vector
        P0_v = vec.array({'px':np.zeros(n),'py':np.zeros(n),'pz':ki_v,'m':np.ones(n)*me})
        # final e 4-momentum vector
        P1_v = vec.array({'px':kf_lab.px,'py':kf_lab.py,'pz':kf_lab.pz,'m':np.ones(n)*me})
        # initial target (assumed at rest)
        A_v = vec.array({'px':np.zeros(n),'py':np.zeros(n),'pz':np.zeros(n),'m':np.ones(n)*MD})
        # final detected system 4-momentum vector
        X_v = vec.array({'px':pf_lab.px,'py':pf_lab.py,'pz':pf_lab.pz,'m':np.ones(n)*MP})
        Q_v = P0_v - P1_v   # q 4-vector
        A1_v = A_v + Q_v    # final target 4-momentum based on conservation
        B_v = A1_v - X_v    # undetected system 4-momentum
        
        # check Q2 and omega
        # Q2_calc_v = Q_v.p*Q_v.p - Q_v.e*Q_v.e
        Q2_calc_v = -Q_v.tau2
        om_calc_v = Q_v.energy

        # Missing momentum components (lab)
        Pmx_lab_v = B_v.px
        Pmy_lab_v = B_v.py
        Pmz_lab_v = B_v.pz
        
        # Check if calculated Pm_v is the same as SIMC pm_v
        Pm_v_calc = np.sqrt(Pmx_lab_v*Pmx_lab_v + Pmy_lab_v*Pmy_lab_v + Pmz_lab_v*Pmz_lab_v)
        
        # e and p in-plane angles
        the_v = kf_lab.theta
        thp_v = pf_lab.theta
        
        #another check for Q2
        Q2_calc2 = 4*P0_v.e*P1_v.e*(np.sin(the_v/2.)**2)
        
        # Now rotate recoil system from lab to q
        thq_v = Q_v.theta  # in-plane angle of q in lab
        # phq_v = Q_v.phi    # out-of-plane angle of q in lab
        
        #get x and b vectors respect to q by rotating the z axis in the lab
        # to q, this method uses a loop over entries 
        Rot_to_q = align_vector_pairs(Q_v.to_Vector3D(),P0_v.to_Vector3D())
        
        xq_v = X_v.transform3D(Rot_to_q)
        bq_v = B_v.transform3D(Rot_to_q)
        # xq_v = X_v.rotateZ(thq_v)   # detected system rotated to q frame
        # bq_v = B_v.rotateZ(thq_v)   # undetected system rotated to q frame
        
        # angles of rotated recoil system
        thpq_v = xq_v.theta
        phpq_v = xq_v.phi
        thnq_v = bq_v.theta
        phnq_v = bq_v.phi
        
        # check th_nq
        costhnq = (Q_v.p**2 + Pm_v_calc**2 - Pf_v**2)/(2.*Q_v.p*Pm_v_calc)
        thnq_calc = np.arccos(costhnq)
        
        # check th_pq
        Pf_par = ( Pf_v**2 + Q_v.p**2 - Pm_v_calc**2)/ (2.*Q_v.p)
        costhpq = Pf_par/Pf_v
        thpq_calc = np.arccos(costhpq)
        
        # check phi_pq
        # y_sct = P0_v.to_Vector3D().unit().cross(P1_v.to_Vector3D().unit())
        # y_react = Q_v.to_Vector3D().unit().cross(pf_lab.unit())
        # cosphipq = y_sct.dot(y_react)
        # phpq_calc = np.arccos(cosphipq)
        
        # wrong thnq
        # xq_v = X_v.rotateZ(thq_v)   # detected system rotated to q frame
        # bq_v = B_v.rotateZ(thq_v)
        
        # wrong_thnq = bq_v.theta
        
        # missing momentum in q frame
        Pmq_v = -bq_v
        Pmx_q_v = Pmq_v.px
        Pmy_q_v = Pmq_v.py
        Pmz_q_v = Pmq_v.pz
        
        # Compare q_v_calc with SIMC q_lab_v
        # Compare Pm_v_calc with SIMC pm_v
        # convert all angles to degrees here
        kin = {'xbj_v':xbj_v,
               'q_v_calc':Q_v.p,
               'pm_v_calc':Pm_v_calc,
               'the_v':the_v/dtr,
               'thp_v':thp_v/dtr,
               'th_pq_v':thpq_v/dtr,
               'ph_pq_v':phpq_v/dtr,
               'th_nq_v':thnq_v/dtr,
               'ph_nq_v':phnq_v/dtr,
               'Pm_q_v':Pmq_v,
               'Pmx_q_v':Pmx_q_v,
               'Pmy_q_v':Pmy_q_v,
               'Pmz_q_v':Pmz_q_v,
               'kf_v':kf_v,
               'thq_v':thq_v/dtr,
               'pm_v_GeV':Pm_v,
               'Q2_v_GeV2':Q2_v,
               'nu_v_GeV':omega_v,
               'q_lab_v_GeV':q_v,
               'pf_v_GeV':Pf_v,
               'Ein_v_GeV':Ein_v,
               'Q2_calc':Q2_calc_v,
               'nu_calc':om_calc_v,
               'Q2_calc2':Q2_calc2,
               'thnq_calc':thnq_calc/dtr,
               'thpq_calc':thpq_calc/dtr
               }
        # update data_init object with new branches
        data_obj.Branches.update(kin)
        print('data_init object updated with the following branches:',
              kin.keys())
        