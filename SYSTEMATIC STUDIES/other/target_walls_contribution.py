import LT.box as B
import numpy as np

import data_init as D
import cut_handler as C

# utility functions
###
# function to combine histograms from different runs, will be normalized by
# efficiencies per run, then divided by total charge of combined runs.
# updated to deal with different data sets
###        
        
def comb_n_norm_histos(histos,many=[],rt='coin'):   
    char = []
    first = True
    for m in many:
        enorm = D.get_eff_norm(m,run_type=rt)
        char.append(D.get_charge_norm(m,curr='BCM1'))
        
        if first:
            enorm_h = histos[m]/enorm
            first = False
        else:
            enorm_h += histos[m]/enorm
        
    tot_char = np.sum(np.array(char))
    combined_histos = enorm_h/tot_char
    
    return combined_histos

# runs to analyze: pm580 set_1, Al_dummy run 20857
R = [20873,20874,20875,20876,20877,20878,20880,20881,20882,20883]

# branches to get
T_branches = ['H.kin.secondary.emiss_nuc','P.react.z','P.gtr.dp','H.gtr.dp',
              'H.extcor.xsieve','H.extcor.ysieve','CTime.epCoinTime_ROC2',
              'P.extcor.xsieve','P.extcor.ysieve']

# load branches in RUN from root_DIR directory
root_DIR = "/media/gvill/Gema's T7/ROOTfiles/pass_4/"

RUN = D.DATA_INIT(data_type='deut23_data', run=R, 
                    select_branches={'T':T_branches},ROOTfiles_path=root_DIR)

#%%
# select cuts to apply and combine them all into all cuts
cuts_list = [C.hms_delta, C.shms_delta] # std acceptance cuts  
cuts_to_apply = {}
for r in RUN.many:
    c_list = []
    print(f'Cuts for Run {r}')
    for cut in cuts_list:
        cut.init()
        br = RUN.Branches[r][C.HCANA_names[cut.name]]
        cut_array = cut(br)
        cut.stats()
        c_list.append(cut_array)
        
    cuts_to_apply[r] = c_list

# HMS collimator cut    
hcoll_cut = C.coll_cut(RUN,spec='HMS',many=True)
hcoll_cut_arrays = hcoll_cut()
for r in RUN.many:
    cuts_to_apply[r].append(hcoll_cut_arrays[r]) 

# # SHMS collimator cut
# pcoll_cut = C.coll_cut(RUN,spec='SHMS',many=True)
# pcoll_cut_arrays = pcoll_cut()
# for r in RUN.many:
#     cuts_to_apply[r].append(pcoll_cut_arrays[r]) 

# # coincidence time cut
ct_cut = C.CTime_cut(RUN,many=True,cut_bounds=C.CT_CUT_BOUND)
ct_cut_arrays = ct_cut() 
for r in RUN.many:
    cuts_to_apply[r].append(ct_cut_arrays[r]) 

# combine all cuts here
all_cuts = {}
for r in RUN.many:
    all_cuts_arr = cuts_to_apply[r][0]
    for arr in cuts_to_apply[r]:    
        all_cuts_arr = all_cuts_arr & arr    
    
    all_cuts[r] = all_cuts_arr 
#%%
# define aluminum endcap cut
Al_endcap_cut = C.WCUT(-0.2,-0.05,name='Al_endcaps_cut')

# get desired branches
Em_h = {}
Zt_h = {}
Zt_h_cut = {}
for m in RUN.many:
    Em_nuc = RUN.Branches[m]['H.kin.secondary.emiss_nuc'][all_cuts[m]]
    Z_tar = RUN.Branches[m]['P.react.z'][all_cuts[m]]
    
    Em_h[m] = B.histo(Em_nuc, range=(-0.2,0.1), bins=50)
    Zt_h[m] = B.histo(Z_tar, range=(-10,10), bins=100)
    
    Al_endcap_cut_arr = Al_endcap_cut(Em_nuc)
    Z_tar_cut = Z_tar[Al_endcap_cut_arr]
    
    Zt_h_cut[m] = B.histo(Z_tar_cut, range=(-10,10), bins=100)

# these are the histograms with std cuts only, combined and normalized 
all_Em_h = comb_n_norm_histos(Em_h,RUN.many)
all_Zt_h = comb_n_norm_histos(Zt_h,RUN.many)    

# this is the Ztar histogram with the Al endcap cut
all_Zt_h_cut = comb_n_norm_histos(Zt_h_cut,RUN.many)  

#%%
#plot histograms
B.pl.figure()
all_Em_h.plot()
B.pl.ylim(0.0015,0.5)
B.pl.yscale('log')

B.pl.figure()
all_Zt_h.plot()
all_Zt_h_cut.plot()

B.pl.ylim(0.0001,0.5)
B.pl.yscale('log')

all_Zt_sum = all_Zt_h.sum(-6.5,6.5)
Zt_cut_sum = all_Zt_h_cut.sum(-6.5,6.5)

tgt_wall_contr = Zt_cut_sum[0]/all_Zt_sum[0]

print('tgt wall contribution = ', tgt_wall_contr)
    
