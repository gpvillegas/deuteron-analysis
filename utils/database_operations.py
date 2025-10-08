"""
WB SQLITe data base initialization and operations to hande parameters needed for data analysis

- create('/Users/boeglinw/my_db.db'):  creates an initial data base

"""
import sqlite3 as lite

import os
import sys


import re
#import os

#DATA_BASE_DIR = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/DataAnalysis_Diamond/'
#DATA_BASE_DIR= ''
DATA_BASE_DIR= '/home/gvill/deuteron/deuteron_db/'

NOSHOT = -99999  # Default number for new shot
DB_ERROR = None

#%% connect to db
def connect_sqlite(dbfile, new = False):
    if os.path.isfile(dbfile) or new :
        return lite.connect(dbfile)
    else:
        raise FileNotFoundError
        return None
        
#%%
# perpare data to be stored in data base
class db_table_data:
    def __init__(self, name, fields, types, values = None, special = None):
        self.name = name
        if values == None:
            self.values = [dict.fromkeys(fields)]
        else:
            self.values = [dict(zip(fields, v)) for v in values]  # convert array of values to array of dict
        self.field_names = self.values[0].keys()
        self.special = special
        self.types = types


    def create_table(self, db_connect):
        global DB_ERROR
        cur = db_connect.cursor()
        sql = f'CREATE TABLE IF NOT EXISTS {self.name} ('
        for i, s in enumerate(self.field_names):
            if i>0:
                sql += ', '
            sql += f'{s} {self.types[i]}'
        # add special commands, e.f.define primary key
        if self.special is not None:
            sql += ', '
            sql += self.special
        sql += ')'
        try:
            cur.execute(sql)
        except Exception as err:
                DB_ERROR = err
                print(f'---> create_table problem with: {sql}, {err}')
        DB_ERROR = None

    def insert_into(self, db_connect):
        global DB_ERROR
        cur = db_connect.cursor()
        for d in self.values:
            sql = f'INSERT INTO {self.name} VALUES ('
            for i,s in enumerate(self.field_names):
                if i>0:
                    sql += ', '
                sql += f'{d[s]}'
            sql += ')'
            try:
                cur.execute(sql)
            except Exception as err:
                DB_ERROR = err
                print(f'---> insert_into problem with: {sql}, {err}')
            DB_ERROR = None
        # all done

#%%
#creates a starting database with necessary tables and fields
def cast_to_list(obj):
    if type(obj) == list:
        return obj
    else:
        return [obj]

def hcana_dict(b=None):   
    # assess calibration status:
    # Reference Time branches
    hREF_branches = ['T.coin.hDCREF1_tdcMultiplicity','T.coin.hDCREF1_tdcTime',
            'T.coin.hDCREF1_tdcTimeRaw','T.coin.hDCREF2_tdcMultiplicity',
            'T.coin.hDCREF2_tdcTime','T.coin.hDCREF2_tdcTimeRaw',
            'T.coin.hDCREF3_tdcMultiplicity','T.coin.hDCREF3_tdcTime',
            'T.coin.hDCREF3_tdcTimeRaw','T.coin.hDCREF4_tdcMultiplicity',
            'T.coin.hDCREF4_tdcTime','T.coin.hDCREF4_tdcTimeRaw',
            'T.coin.hDCREF5_tdcMultiplicity','T.coin.hDCREF5_tdcTime',
            'T.coin.hDCREF5_tdcTimeRaw', 'T.coin.hT1_tdcMultiplicity',
            'T.coin.hT1_tdcTime','T.coin.hT1_tdcTimeRaw',
            'T.coin.hT2_tdcMultiplicity','T.coin.hT2_tdcTime',
            'T.coin.hT2_tdcTimeRaw','T.coin.hFADC_TREF_ROC1_adcMultiplicity',
            'T.coin.hFADC_TREF_ROC1_adcPulseTime',
            'T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw']
    
    pREF_branches = ['T.coin.pDCREF10_tdcMultiplicity','T.coin.pDCREF10_tdcTime',
              'T.coin.pDCREF10_tdcTimeRaw','T.coin.pDCREF1_tdcMultiplicity',
              'T.coin.pDCREF1_tdcTime','T.coin.pDCREF1_tdcTimeRaw',
              'T.coin.pDCREF2_tdcMultiplicity','T.coin.pDCREF2_tdcTime',
              'T.coin.pDCREF2_tdcTimeRaw','T.coin.pDCREF3_tdcMultiplicity',
              'T.coin.pDCREF3_tdcTime','T.coin.pDCREF3_tdcTimeRaw',
              'T.coin.pDCREF4_tdcMultiplicity','T.coin.pDCREF4_tdcTime',
              'T.coin.pDCREF4_tdcTimeRaw','T.coin.pDCREF5_tdcMultiplicity',
              'T.coin.pDCREF5_tdcTime','T.coin.pDCREF5_tdcTimeRaw',
              'T.coin.pDCREF6_tdcMultiplicity','T.coin.pDCREF6_tdcTime',
              'T.coin.pDCREF6_tdcTimeRaw','T.coin.pDCREF7_tdcMultiplicity',
              'T.coin.pDCREF7_tdcTime','T.coin.pDCREF7_tdcTimeRaw',
              'T.coin.pDCREF8_tdcMultiplicity','T.coin.pDCREF8_tdcTime',
              'T.coin.pDCREF8_tdcTimeRaw','T.coin.pDCREF9_tdcMultiplicity',
              'T.coin.pDCREF9_tdcTime','T.coin.pDCREF9_tdcTimeRaw',
              'T.coin.pT1_tdcMultiplicity','T.coin.pT1_tdcTime',
              'T.coin.pT1_tdcTimeRaw','T.coin.pT2_tdcMultiplicity',
              'T.coin.pT2_tdcTime','T.coin.pT2_tdcTimeRaw',
              'T.coin.pFADC_TREF_ROC2_adcMultiplicity',
              'T.coin.pFADC_TREF_ROC2_adcPulseTime',
              'T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw']
    
    
        
    # HODOSCOPE branches
    Phodo_branches = ['P.hod.beta','P.hod.betanotrack','P.hod.betachisqnotrack',
                     'P.hod.adctdc_offset','P.hod.2y.totNumTdcHits',
                     'P.hod.2y.totNumAdcHits','P.hod.2y.TrackXPos',
                     'P.hod.2y.TrackYPos','P.hod.2x.totNumTdcHits',
                     'P.hod.2x.totNumAdcHits','P.hod.2x.TrackXPos',
                     'P.hod.2x.TrackYPos','P.hod.1y.totNumTdcHits',
                     'P.hod.1y.totNumAdcHits','P.hod.1y.TrackXPos',
                     'P.hod.1y.TrackYPos','P.hod.1x.totNumTdcHits',
                     'P.hod.1x.totNumAdcHits','P.hod.1x.TrackXPos',
                     'P.hod.1x.TrackYPos']
    
    Hhodo_branches = ['H.hod.beta','H.hod.betanotrack','H.hod.betachisqnotrack',
                     'H.hod.adctdc_offset','H.hod.2y.totNumTdcHits',
                     'H.hod.2y.totNumAdcHits','H.hod.2y.TrackXPos',
                     'H.hod.2y.TrackYPos','H.hod.2x.totNumTdcHits',
                     'H.hod.2x.totNumAdcHits','H.hod.2x.TrackXPos',
                     'H.hod.2x.TrackYPos','H.hod.1y.totNumTdcHits',
                     'H.hod.1y.totNumAdcHits','H.hod.1y.TrackXPos',
                     'H.hod.1y.TrackYPos','H.hod.1x.totNumTdcHits',
                     'H.hod.1x.totNumAdcHits','H.hod.1x.TrackXPos',
                     'H.hod.1x.TrackYPos']
    
    # DRIFT CHAMBER branches
    Pdc_branches = ['P.dc.x_fp','P.dc.y_fp','P.dc.residualExclPlane',
                    'P.dc.2x2.dist', 'P.dc.2x1.dist', 'P.dc.2v2.dist', 
                    'P.dc.2v1.dist', 'P.dc.2u2.dist', 'P.dc.2u1.dist',
                    'P.dc.1x2.dist', 'P.dc.1x1.dist', 'P.dc.1v2.dist',
                    'P.dc.1v1.dist', 'P.dc.1u2.dist', 'P.dc.1u1.dist'] 
    
    Hdc_branches = ['H.dc.x_fp','H.dc.y_fp','H.dc.residualExclPlane',
                    'H.dc.2x2.dist', 'H.dc.2x1.dist', 'H.dc.2v2.dist', 
                    'H.dc.2v1.dist', 'H.dc.2u2.dist', 'H.dc.2u1.dist',
                    'H.dc.1x2.dist', 'H.dc.1x1.dist', 'H.dc.1v2.dist',
                    'H.dc.1v1.dist', 'H.dc.1u2.dist', 'H.dc.1u1.dist']
    
    # CALORIMETER branches
    Pcal_branches = ['P.cal.etottracknorm', 'P.cal.etotnorm', 'P.cal.ytrack',
                     'P.cal.xtrack']
    
    Hcal_branches = ['H.cal.etottracknorm', 'H.cal.etotnorm', 'H.cal.ytrack',
                     'H.cal.xtrack']
    
    # CERENKOV branches
    Phgcer_branches = ['P.hgcer.npeSum', 'P.hgcer.npe']
    
    Pngcer_branches = ['P.ngcer.npeSum', 'P.ngcer.npe']
    
    Hcer_branches = ['H.cer.npeSum', 'H.cer.npe']
    
    # KINEMATIC branches
    
    kin_branches= ['H.kin.secondary.Erecoil','H.kin.secondary.MMK',
                   'H.kin.secondary.MMp','H.kin.secondary.MMpi',
                   'H.kin.secondary.MandelS','H.kin.secondary.MandelT',
                   'H.kin.secondary.MandelU','H.kin.secondary.Mrecoil',
                   'H.kin.secondary.Prec_x','H.kin.secondary.Prec_y',
                   'H.kin.secondary.Prec_z','H.kin.secondary.emiss',
                   'H.kin.secondary.emiss_nuc','H.kin.secondary.ph_bq',
                   'H.kin.secondary.ph_xq','H.kin.secondary.phb_cm',
                   'H.kin.secondary.phx_cm','H.kin.secondary.pmiss',
                   'H.kin.secondary.pmiss_x','H.kin.secondary.pmiss_y',
                   'H.kin.secondary.pmiss_z','H.kin.secondary.px_cm',
                   'H.kin.secondary.t_tot_cm','H.kin.secondary.tb',
                   'H.kin.secondary.tb_cm','H.kin.secondary.th_bq',
                   'H.kin.secondary.th_xq','H.kin.secondary.thb_cm',
                   'H.kin.secondary.thx_cm','H.kin.secondary.tx',
                   'H.kin.secondary.tx_cm','H.kin.secondary.xangle',
                   'P.kin.primary.Q2','P.kin.primary.W','P.kin.primary.W2',
                   'P.kin.primary.epsilon','P.kin.primary.nu',
                   'P.kin.primary.omega','P.kin.primary.ph_q','P.kin.primary.q3m',
                   'P.kin.primary.q_x','P.kin.primary.q_y','P.kin.primary.q_z',
                   'P.kin.primary.scat_ang_deg','P.kin.primary.scat_ang_rad',
                   'P.kin.primary.th_q','P.kin.primary.x_bj']
    
    # GOLDEN TRACK branches
    
    gtr_branches= ['H.dc.gtrack_nsp','H.gtr.beta','H.gtr.dp','H.gtr.index',
                   'H.gtr.ok','H.gtr.p','H.gtr.ph','H.gtr.px','H.gtr.py',
                   'H.gtr.pz','H.gtr.th','H.gtr.x','H.gtr.y','P.dc.gtrack_nsp',
                   'P.gtr.beta','P.gtr.dp','P.gtr.index','P.gtr.ok','P.gtr.p',
                   'P.gtr.ph','P.gtr.px','P.gtr.py','P.gtr.pz','P.gtr.th',
                   'P.gtr.x','P.gtr.y']
    
    # REACTION VERTEX branches
    
    react_branches= ['H.react.ok','H.react.x','H.react.y','H.react.z',
                     'P.react.ok','P.react.x','P.react.y','P.react.z']
    
    # CUT VARIABLES (not included in the groups above)
    cut_branches= ['P.gtr.dp','H.gtr.dp','P.react.z','H.react.z',
                   'P.gtr.th','P.gtr.ph','H.gtr.th','H.gtr.ph']
    
    # optics variables: x/yfp, x/yptar, dp, x/ysieve
    optics_branches= ['P.dc.x_fp','P.dc.y_fp','P.dc.xp_fp','P.dc.yp_fp',
                      'H.dc.x_fp','H.dc.y_fp','H.dc.xp_fp','H.dc.yp_fp',
                      'P.gtr.th','P.gtr.ph','H.gtr.th','H.gtr.ph',
                      'P.extcor.xsieve','P.extcor.ysieve',
                      'H.extcor.xsieve','H.extcor.ysieve',
                      'P.gtr.dp','H.gtr.dp']
    
    branch_types= {'phodo':Phodo_branches,
                     'hhodo':Hhodo_branches,
                     'pdc':Pdc_branches,
                     'hdc':Hdc_branches,
                     'pcal':Pcal_branches,
                     'hcal':Hcal_branches,
                     'phgcer':Phgcer_branches,
                     'pngcer':Pngcer_branches,
                     'hcer':Hcer_branches,
                     'kin':kin_branches,
                     'gtr': gtr_branches,
                     'react': react_branches,
                     'cut':cut_branches,
                     'hREF':hREF_branches,
                     'pREF':pREF_branches,
                     'optics':optics_branches}

    calib_group =  branch_types['phodo'] + branch_types['hhodo'] +\
        branch_types['pdc'] +  branch_types['hdc'] + branch_types['pcal'] +\
        branch_types['hcal'] + branch_types['phgcer'] +\
        branch_types['pngcer'] + branch_types['hcer']
    
    kins_group = branch_types['kin'] + branch_types['cut']
    
    reftime_group = branch_types['hREF'] + branch_types['pREF']
    
    branches_groups= {'calib': calib_group, 'kins': kins_group,
                      'reftime': reftime_group}
    
    all_branches = []
    for br in branch_types:
        all_branches = all_branches + branch_types[br]

    if b is list:
        list_b = cast_to_list(b)
        branches_sel = []
        for name in list_b:
            branches_sel = branches_sel + branches_groups[name]        
        return branches_sel
    elif b == 'types':
        return branch_types
    elif b == 'groups':
        return branches_groups
    else: 
        return all_branches

#%%
def create(db_file):
    conn = connect_sqlite(db_file, new = True)
    # default values for creating tables:
    common_parameters_fields = ['Root_Folder']
    common_parameters_types = ['TEXT']
    common_parameters_values = ['"./"']
    
    common_parameters = db_table_data('Common_Parameters', common_parameters_fields, common_parameters_types, [common_parameters_values])
    
    
    # Shot list table
    shot_list_fields = ['Shot',         'Date', 'File_Name',' Folder', 'RP_position', 'RP_setpoint','t_offset', 'N_chan', 'Comment']
    shot_list_types =  ['INT not NULL', 'TEXT', 'TEXT',      'TEXT',   'REAL',        'REAL',         'REAL',   'INT',    'TEXT']
    shot_list_values = []
    shot_list_values.append( [29975,      '"22-Aug-2013"', '"29975_DAQ_220813-141746.hws"','"Data/"', 1.65, 0., 0., 6,   '"No comment"' ] )  # strings need to be enclosed in ""
    shot_list_values.append ([29879,      '"19-Aug-2013"', '"DAQ_190813-112521.hws"','"Data/"', 1.83, 0., 0., 6,   '"No comment"' ])   # strings need to be enclosed in ""
    shot_list_values.append([29880,      '"19-Aug-2013"', '"DAQ_190813-114059.hws"','"Data/"', 1.65, 0., 0., 6,   '"No comment"' ])   # strings need to be enclosed in "
    shot_list_values.append( [99999,      '"22-Aug-2013"', '"29975_DAQ_220813-141746.hws"','"Data/"', 1.65, 0., 0., 6,   '"No comment"' ] )  # strings need to be enclosed in ""

    shot_list = db_table_data('Shot_List', shot_list_fields, shot_list_types, shot_list_values, special = 'PRIMARY KEY (Shot)') 

    # raw_fitting_table
    raw_fit_fields = ['Shot', 'Channel', 'Version',               'Comment', 'dtmin', 'dtmax', 'n_peaks_to_fit', 'poly_order', 'add_pulser', 'pulser_rate','P_amp', 'use_threshold', 'Vth', 'Vstep', 'n_sig_low', 'n_sig_high', 'n_sig_boundary', 'sig', 'min_delta_t', 'max_neg_V', 'Result_File_Name', 'Corrected_Data_File_Name', 'Input_File_Name']
    raw_fit_types = ['INT not NULL','INT not NULL','INT not NULL', 'TEXT',    'REAL',  'REAL',  'INT',            'INT',        'TEXT',       'REAL',       'REAL',  'TEXT',          'REAL','REAL',  'REAL',     'REAL',       'REAL'          ,  'REAL',  'REAL',    'REAL',           'TEXT',    'TEXT',    'TEXT']
    raw_fit_values = [99999,          0,             0,        '"No Comment"', 0.01,   0.1,     10,              10,           '"True"',        1000.,        1.0,    '"True"',          0.2,   0.2,       3.,           3.,         3.,              0.3,      1e-7,       -.3 ,        '"No File Saved"', '"No File Saved"', '" "']

    raw_fitting = db_table_data('Raw_Fitting', raw_fit_fields, raw_fit_types, [raw_fit_values], special = 'PRIMARY KEY (Shot, Channel, Version)')  # if only 1 row of values enter in backets

    # peak sampling
    peak_sampling_fields = ['Shot', 'Channel', 'Version',               'Comment', 'Vstep', 'Vth', 'Chi2', 'tmin', 'tmax', 'decay_time', 'rise_time', 'position', 'n_samp', 'n_max', 'n_below', 'n_above']
    peak_sampling_types = ['INT not NULL','INT not NULL','INT not NULL', 'TEXT',   'REAL',  'REAL', 'REAL', 'REAL', 'REAL', 'REAL',       'REAL',      'REAL',     'INT',    'INT',    'INT',     'INT']
    peak_sampling_values = [99999,          0,             0,        '"No Comment"', .1,     .3,     .1,     .1,    .15,    200e-9,      100e-9,      350e-9,    120,       20,      15,        50]


    peak_sampling_fields += ['b1',  'e1',  'b2',  'e2',  'b3',  'e3',  'b4',  'e4',  'b5',  'e5',  'b6',  'e6']
    peak_sampling_types +=  ['REAL','REAL','REAL','REAL','REAL','REAL','REAL','REAL','REAL','REAL','REAL','REAL']
    peak_sampling_values += [0.0264465750973, 0.0264511717629,0.0264640275167, 0.0264695346329, 0.0265033656474, 0.0265088727637, 0.0265219299587, 0.0265274370749, 0.0265635664222, 0.0265690735384,0.0266368021858, 0.0266423093020]

    peak_sampling_fields += ['b7',  'e7',  'b8',  'e8',  'b9',  'e9',  'b10',  'e10',  'b11',  'e11',  'b12',  'e12']
    peak_sampling_types +=  ['REAL','REAL','REAL','REAL','REAL','REAL','REAL', 'REAL', 'REAL', 'REAL', 'REAL', 'REAL']
    peak_sampling_values += [0.0266432752680, 0.0266487823842,0.0266662141836, 0.0266717212999, 0.0267147567486, 0.0267202638648, 0.0268340257862, 0.0268395329025, 0.0268832345346, 0.0268887416508,0.0268869429556, 0.0268924500719]

    peak_sampling =db_table_data('Peak_Sampling', peak_sampling_fields, peak_sampling_types, [peak_sampling_values], special = 'PRIMARY KEY (Shot, Channel, Version)')

    # rate plotting
    Rate_Analysis_fields = ['Shot', 'Channel', 'Version',               'Comment', 'time_slice_width', 'h_min', 'h_max' , 'h_bins', 'draw_p', 'draw_t', 'draw_pul' , 'draw_sum']
    Rate_Analysis_types = ['INT not NULL','INT not NULL','INT not NULL', 'TEXT',   'REAL',             'REAL',  'REAL',   'INT',   'TEXT',   'TEXT',   'TEXT',        'TEXT']
    Rate_Analysis_values = [99999,          0,             0,        '"No Comment"', 1.e-3,            0.,      1.4,      160,      '"True"',  '"False"',  '"False"',  '"False"']

    Rate_Analysis_fields+= ['p_min', 'p_max', 't_min', 't_max', 'pul_min', 'pul_max',  'sig_ratio', 'Result_File_Name', 'Input_File_Name']
    Rate_Analysis_types += ['REAL',  'REAL',  'REAL',  'REAL',  'REAL',    'REAL',       'REAL', 'TEXT', 'TEXT' ]
    Rate_Analysis_values += [0.6,     0.8,    0.3,      0.5,     0.9,        1.2,        100,  '"No File Saved"',  '" "']

    Rate_Analysis = db_table_data('Rate_Analysis', Rate_Analysis_fields, Rate_Analysis_types, [Rate_Analysis_values], special = 'PRIMARY KEY (Shot, Channel, Version)')


    # corrected shot list
    shot_list_corrected_fields = ['Shot',       'Channel',     'Version','      Iteration',   'File_Name', 'Folder', 'Comment']
    shot_list_corrected_types = ['INT not NULL','INT not NULL','INT not NULL', 'INT not NULL', '""',       '""',      '""']
    shot_list_corrected_values = [99999,          0,             0,              0,              '"No File"',  '"NO Folder"', '"No comment"']
    
    
    shot_list_corrected = db_table_data('Shot_List_Corrected', shot_list_corrected_fields, shot_list_corrected_types, values = [shot_list_corrected_values], special = 'PRIMARY KEY (Shot, Channel, Version, Iteration)' )

    # ocmbine rates 
    comb_rates_fields = ['Shot', 'Channels', 't_min', 't_max', 'A_min', 'A_max', 'd_time', 'view_dir', 'view_names', 'r_min', 'r_max', 'use_all_variables', 'calc_rate', 'model', 'Comment']
    comb_rates_types = [ 'INT',  'TEXT',     'REAL',  'REAL',  'REAL',  'REAL',  'REAL',   'TEXT',     'TEXT',       'REAL',  'REAL',  'TEXT',              'TEXT',      'TEXT',  'TEXT']
    comb_rates_values = [99999,  '"0,1,2,3"',  0.,       0.5,     0.,      150.e3,  5.e-3,    '"./orbit_public/NSTX_output"', '"nml_orb_NSTX-Case_3_0.3"', 0., 1.5, '"True"',  '"True"', '"Simple Gauss"', '"No Comment"']

    comb_rates =db_table_data('Combined_Rates', comb_rates_fields, comb_rates_types, [comb_rates_values])

    with conn:
        
        common_parameters.create_table(conn)
        common_parameters.insert_into(conn)
        
        shot_list.create_table(conn)
        shot_list.insert_into(conn)
        #
        raw_fitting.create_table(conn)
        raw_fitting.insert_into(conn)
        #
        peak_sampling.create_table(conn)
        peak_sampling.insert_into(conn)
        #
        Rate_Analysis.create_table(conn)
        Rate_Analysis.insert_into(conn)
        #
        comb_rates.create_table(conn)
        comb_rates.insert_into(conn)
        #
        shot_list_corrected.create_table(conn)
        shot_list_corrected.insert_into(conn)
        
    print("Created DB")
    conn.close()

def create_HcanaBranches(db_file):
    conn = connect_sqlite(db_file)
    
    # default values for creating tables:
    common_parameters_fields = ['Root_Folder']
    common_parameters_types = ['TEXT']
    common_parameters_values = ['"./"']
    
    common_parameters = db_table_data('Common_Parameters', 
                                      common_parameters_fields,
                                      common_parameters_types, 
                                      [common_parameters_values])
    
    # Table with all T Tree Branches
    T_branches_fields = ['Branches','Description','Type','Groups','Comment']
    T_branches_types = ['TEXT','TEXT','TEXT','TEXT','TEXT']
    T_branches_values = []
    with open('/home/gvill/deuteron/root_analysis/hcana_branches.txt','r')\
        as f:
        for line in f:
            T_branches_values.append([f'"{line.strip()}"','""','""','""','""'])
    
    TTree_branches = db_table_data('TTree_Branches', T_branches_fields, 
                                   T_branches_types, T_branches_values,
                                   special = 'PRIMARY KEY(Branches)')
    
    with conn:
        TTree_branches.create_table(conn)
        TTree_branches.insert_into(conn)
        
    print("Created DB")
    conn.close()


#%% create help text table 
#  this function extracts the table information and creates a new one with the same column names but with
# data type TEXT

def make_help_table(db_file, table):
    """
    create a help type table

    Parameters
    ----------
    db_file : str
        data base file
    table : str 
        table name to create a help table from.

    Returns
    -------
    None.

    """
    col_names = []
    col_types = []
    
    qline = f"PRAGMA table_info('{table}');"
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(qline)
        except Exception as e:
            print(f'probem with command: {qline}')
            print(f'Error: {e}')
            return 
        rows = cur.fetchall()
    conn.close()
    print(rows)
    if rows == []:
        print(f'No table information received from {qline}')
        return 
    else:
        for row in rows:
            col_names.append(row[1])
            col_types.append('TEXT')
    # setup column data
    help_table = db_table_data(f'{table}_help', col_names, col_types)
    
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    with conn:
        help_table.create_table(conn)
        print(f'{table}_help was created')
    conn.close()
    

#%% Get table information

def get_table_information(db_file, table):
    global DB_ERROR
    qline = f"PRAGMA table_info('{table}');"
    rows = []
    try:
        conn = connect_sqlite(DATA_BASE_DIR + db_file)
        with conn:
            cur = conn.cursor()
            cur.execute(qline)
            rows = cur.fetchall()
    except Exception as e:
        DB_ERROR = e
        print(f"Error in accessing table information {e}")
    return rows         


#%% get list of tables

def get_list_of_tables(db_file):
    global DB_ERROR
    qline = "SELECT name FROM sqlite_master WHERE type ='table' AND name NOT LIKE 'sqlite_%';"
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(qline)
            res = cur.fetchall()
        except Exception as e:
            print(f'probem with command: {qline}')
            print(f'Error: {e}')
            res = []
    conn.close()
    return [r[0] for r in res]
    

#%% check if data exist
def check_condition(db_file, table, where):
    global DB_ERROR
    """
    Check if a record satsifying the condition where exists

    Parameters
    ----------
    db_file : str
        database file name.
    table : str
        table name.
    where : str
        SQL string for conditional selection.

    Returns
    -------
    bool
        True: a record exists
        False: no reecord exists

    """
    qline = f'select 1 from {table} where exists( select 1 from {table} where {where})'
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(qline)
        except Exception as e:
            print(f'probem with command: {qline}')
            print(f'Error: {e}')
            return False
        res = cur.fetchall()
    conn.close()
    if res == []:
        return False
    else:
        return True

#%% raise sys.exit with error message if data do not exist

def check_data(db_file, table, where, message = None):
    """
    Check if a record satsifying the condition where exists.
    If not raises sys.exit with custom error message

    Parameters
    ----------
    db_file : str
        database file name.
    table : str
        table name.
    where : str
        SQL string for conditional selection.
    message: str, optional
        Error message to be displayed

    Returns
    -------

    """ 
    global DB_ERROR
    if check_condition(db_file, table, where):
        return
    else:
        if message is None:
            message = f'No data in {table} for condition {where}'
        sys.exit(message)
        

#%%retrieve data
def retrieve(db_file, params, table, where = None, distinct = False):
    global DB_ERROR
    """
    Retrieves paremeters from table in database.

    Parameters
    ----------
    db_file : str
        database file name.
    params : str
        parameter (field) name to retreive.
    table : str
        table name where the parameter is located.
    where : str
        SQL string for conditional selection.

    Returns
    -------
    
    list of values   

    Example: retrieve(db_file, 'Folder, File_Name', 'Shot_List', 'Shot = 29975')
    
    returns:
        
        [('Data/', '29975_DAQ_220813-141746.hws')]

    """
    if distinct:
        qdist = 'DISTINCT '
    else:
        qdist = ''
        
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    #create query line
    if where is None:
        qline='SELECT ' + qdist + params +' FROM ' + table
    else:
        qline='SELECT ' + qdist + params +' FROM ' + table + ' WHERE ' + where
    print(f'Using database file : {DATA_BASE_DIR + db_file}')
    print(f'query = {qline}')
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(qline)
        except Exception as e:
            print(70*'-')
            print(f'probem with command: {qline}')
            print(f'Error: {e}')
            print(70*'-')
            return []
        result = cur.fetchall()
        return result
    conn.close()


#%% retreive a whole row and return a col. name dictionary value

def get_row(db_file, table, where = None, return_dict = True):
    # get table information
    t_info = get_table_information(db_file, table)
    field_names = [s[1] for s in t_info]
    field_types = [s[2] for s in t_info]
    # get the raw data        
    rows = retrieve(db_file, '*', table, where = where, distinct = False)
    if rows == []:
        return None
    else:
        if len(rows) > 1:
            print(f'---> get_row returned {len(rows)} rows only the first will be used !')
        row_values = list(rows[0])
    for i, tt in enumerate(field_types):
        if tt == 'TEXT':
            row_values[i] = "'" + row_values[i] + "'"  # make sure string are handled correctly
    if return_dict:
        return dict(zip(field_names, row_values))
    else:
        return field_names, row_values
    
#%%write to database
def writetodb(db_file, params, table, where):
    global DB_ERROR
    """
    
    Write parameters to data base

    Parameters
    ----------
    db_file : str
        database file name.
    params : str
        parameter (field) name to set.
    table : str
        table containing hte parameter
    where : str
        selection criteria

    Returns
    -------
    None.

    Example:
        
        writetodb(db_file, 'Folder = "Data/"', 'Shot_list', 'Shot = 29975')
        
    """
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    #create query line
    qline='UPDATE '+ table +' SET ' + params + ' WHERE ' + where
    print(f'---->writetodb: {qline}')
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(qline)
        except Exception as e:
            DB_ERROR = e
            print(f'probem with command: {qline}')
            print(f'Error: {e}')
    conn.close()
    DB_ERROR = None

#%% insert a new row, this produces an error if there is a uniqueness error

def insert_row_into(db_file, table, names, values):
    global DB_ERROR
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    cur = conn.cursor()
    p_1 = f'INSERT INTO  {table} ('
    p_2 = ' VALUES ('
    for i,d in enumerate(names):
        p_1 += d + ','
        p_2 += f'{values[i]}' + ','
    p1 = p_1[:-1] + ')'
    p2 = p_2[:-1] + ')'
    sql = p1 + p2
    print(sql)
    with conn:
        try:
            cur.execute(sql)
        except Exception as err:
            print(f'---> insert_into problem with: {sql}, {err}')
    # all done
    last_row = cur.lastrowid
    conn.close()
    return last_row

#%%  update_row
def update_row(db_file, table, names, values, where):
    global DB_ERROR
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    cur = conn.cursor()
    cmd0 = f'UPDATE  {table} SET ' 
    cmd1 = ','.join( [f'{nn} = {values[i]}' for i,nn in enumerate(names)] )
    cmd2 = f' WHERE {where}'
    q_line  = cmd0 + cmd1 + cmd2
    print(f'update_row: {q_line}')
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(q_line)
        except Exception as e:
            print(f'probem with command: {q_line}')
            print(f'Error: {e}')
    conn.close()    
    
    
#%% some utility functions
def get_shot_list(db_file):
    sl = retrieve(db_file, 'Shot, N_chan', 'Shot_List')
    return sl
    

def get_folder(db_file, shot):
    (dd,) = retrieve(db_file, 'Folder', 'Shot_List', f'Shot = {shot}')[0]
    return dd

def set_folder(db_file, shot, dir_name):
    q_what = f'Folder = "{dir_name}"'
    q_table = 'Shot_List'
    q_where = f'Shot = {shot}'
    writetodb(db_file, q_what, q_table, q_where)    


#%% get previous shot
def prevshot(db_file, shot):
    """
    
    Retrieves previous existing shot number, takes one input - shot number, will return
    shot number which exist in database before entered input shot number.

    Parameters
    ----------
    db_file : str
        database file name.
    shot : str
        shot number

    Returns
    -------
    TYPE
        previous shot number
        
    Example: prevshot(29976)

    """
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    try:
        (rowid,) = retrieve(db_file, 'ROWID', 'Shot_List', f'Shot = {shot}')[0]
    except:
        print(f'--> could not find shot {shot} in {db_file}')
        return None
    print(rowid)
    #create query line
    qline='SELECT Shot FROM Shot_List WHERE ROWID = ' + str(rowid-1)
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(qline)
        except Exception as e:
            print(f'probem with command: {qline}')
            print(f'Error: {e}')
        return cur.fetchall()
    conn.close()


#%% Check for Version in Query

def find_version(q, is_where = True, called_from = ''):
    """
    Find a version statement, and return a new statemen without it.
    If is_where == True: q is a where statement from a query
    if is_where == False: q is a substitution statement for db_write

    Parameters
    ----------
    q : strt
        where or substitution statement
    is_where : bool, optional
        if true, q is a where statement. The default is True.

    Returns
    -------
    bool
        q contains a version statement.
    str
        statement without version.
    str
        version statement

    """
    has_version = False
    version_value = ''    
    if is_where:
        # for find an remove Version statement in where statement
        # if there is a version statement in the where statement remove it
        fields = [ss.strip() for ss in re.split('AND', q, flags=re.IGNORECASE)]
        has_v = ['Version' in s for s in fields]
        for i, v in enumerate(has_v):
            if v:
                version_value = fields.pop(i)
                has_version = True
        q_no_version = ' AND '.join(fields)
    else:
        # find version in subs statement
        # get fields 
        qf = [s.strip() for s in q.split(',')]
        # find version
        has_v = ['Version' in s for s in qf]
        for i, v in enumerate(has_v):
            if v:
                version_value = qf.pop(i)
                has_version = True
        q_no_version = ', '.join(qf)
    print(f'-> find_version : ({called_from}): {has_version}, {q_no_version}, {version_value}')
    print('-> find_version: return values')
    return has_version, q_no_version, version_value

#%% duplicate row
def duplicate_row(db_file, table, where_cp):
    print(f'--------------- > duplicate_row intable {table} where {where_cp}')
    global DB_ERROR
    """
    duplicates a row. If it is versioned it increments the version number

    Parameters
    ----------
    db_file : str
        database file name.
    table : str
        table name
    where_cp : str
        where statement to select row

    Returns
    -------
    has_version : Bool
        True if thie is a versioned row
    max_version : Int
        latest version number

    Be careful when using queries with versions

    """
    # get table information
    t_info = get_table_information(db_file, table)
    field_names = [s[1] for s in t_info]
    field_types = [s[2] for s in t_info]    

    # start connection for fyurhter analysis
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    with conn:
        cur = conn.cursor()
        # check if a version field exists'
        has_version = 'Version' in field_names
        # for finding the max row remove the Version statement from the query if it contains one
        found_version, where_cp_no_version, version_statement =  find_version(where_cp, is_where = True, called_from = 'duplcate_row')       
        # if there is a version statement in the where statement remove it
        if found_version:
            print('duplicate_row--->' +  70*'-')
            print(f'{version_statement} in WHERE statement will be ignored,\n'+\
                  'new version will bey assigned automatically')
            print('duplicate_row---> ' + 70*'-')
        # if it is a versioned row get all field names except version
        if has_version:
            # check if record exists
            if not (check_condition(db_file, table, where_cp) ):
                print(f'--> No data for condition: {where_cp} in {table},  cannot duplicate')
                return has_version, -1                    
            # get the row with the largest version number
            query = f'SELECT *, max(Version) FROM {table}' + ' WHERE ' + where_cp_no_version
            max_version_row = list(cur.execute(query).fetchall()[0])
            # increment the version number
            max_version_row[field_names.index('Version')] += 1
            q_names = '('+ ''.join([f'{x},' for x in field_names])[:-1] + ')'
            # convert values to string
            value_list = [f'"{x}",' if field_types[i] == 'TEXT' else f'{x},' for i,x in enumerate(max_version_row[:-1])]
            q_values = '('+ ''.join(value_list)[:-1] + ')'
            # insert new row with new version number
            q_insert = 'INSERT INTO '+ table + ' ' + q_names +' VALUES ' + q_values
            max_version = max_version_row[field_names.index('Version')]
            print(70*'-')      
            print(f'duplicate_row---> New Version = {max_version}')
            print(f'duplicate_row---> new q_insert {q_insert}')
            try:
                cur.execute(q_insert)
            except Exception as err:
                print(f'duplicate_row---> Cannot execute {q_insert} : {err}')
                max_version = -cur.lastrowid
                return True, max_version, where_cp_no_version
        else:
            # get the current selected values
            row_dict = get_row(db_file, table, where = where_cp)
            row_dict['Shot'] = NOSHOT   # set temporary shot number
            names = list(row_dict.keys())
            values = [row_dict[k] for k in row_dict]
            # get the latest row id
            lastrowid = insert_row_into(db_file, table, names, values)            
            max_version = -lastrowid  # make the last rowid negative to idntify the value as return value
    conn.close()
    return has_version, max_version, where_cp_no_version


#%%duplicate certain row in databse
def copy_row(db_file, table, where_cp, substitutions):
    global DB_ERROR
    """
    Creates a copy of the row in table and changes some parameter values in it.

    Parameters
    ----------
    db_file : str
        database file name.
    table : str
        table containing hte parameter.
    where_cp : str
        selection criteria.
    substitutions : str
        parameter setting command.

    Returns
    -------
    None.
    
    Example: copy_row(db_file, 'Raw_Fitting', 'Shot = 29975 AND Channel = 0', 'Channel = 1')

    """
    # check if row exists
    if not (check_condition(db_file, table, where_cp) ):
        print(f'copy_row---> No data for condition: {where_cp} in {table},  cannot copy')
        return
    # check if version is in the where_cp
    # first duplicate row
    has_version, max_version, where_cp_no_version = duplicate_row(db_file, table, where_cp)
    # insert new values
    # find new version number
    found_version, _ , version_statement =  find_version(substitutions, is_where = False, called_from = 'copy_row')
    if has_version:
        vv = int(version_statement.split('=')[-1])
        print(f'copy_row---> desired version = {vv}, {where_cp_no_version}')
        # find version numbers for this query
        versions = [int(r[0]) for r in retrieve(db_file, 'Version', table, where = where_cp_no_version)]
        print(f'copy_row---> versions = {versions}')
        # check if the chosen version value exists already, if so 
        while vv in versions:
            # the desired version exists already, increment it by one and try again
            vv += 1
        # valid version found
        print(f'copy_row---> New version selected = {vv}')
        writetodb(db_file, substitutions, table, where_cp )
    else:
        writetodb(db_file, substitutions, table,  f'ROWID = {-max_version}')
    print(f'----> copy_row: has_version = {has_version}, max_version = {max_version}')

#%% delete rows in a table

def delete_row(db_file, table, where_del):
    global DB_ERROR
    """
    Delete a row in a table 

    Parameters
    ----------
    db_file : str
        database fi
        le name.
    table : str
        table where the row should be deleted.
    where_del : str
        selection criteria for row to be deleted (careful !).

    Returns
    -------
    None.

    """
    conn = connect_sqlite(DATA_BASE_DIR + db_file)
    
    if not (check_condition(db_file, table, where_del) ):
        print(f'--> No data for condition: {where_del} in {table}, nothing done')
        return 
    q_line = f'delete from {table} where {where_del};'
    print(f'{q_line}')
    with conn:
        cur = conn.cursor()
        try:
            cur.execute(q_line)
        except Exception as e:
            print(f'probem with command: {q_line}')
            print(f'Error: {e}')
    conn.close()   
    

#%% Tests
"""
if __name__ == "__main__":
    db_file = 'New_MainDB1.db'
    wheredb_cp = ('Shot = 29975 AND Channel = 0')
    wheredb_cp_new = ('Shot = 29980 AND Channel = 3')
    print(wheredb_cp)
    # create(db_file)
    # try:
    copyrow(db_file, 'Peak_Sampling', wheredb_cp, 'Shot = 29980, Channel = 3')
    print(retrieve(db_file, '*','Peak_Sampling', wheredb_cp))
    print(retrieve(db_file, '*','Peak_Sampling', wheredb_cp_new))
    # except:
    #    print("Didnt work")
"""


