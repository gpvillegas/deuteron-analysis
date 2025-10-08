import numpy as np
import matplotlib.pyplot as plt
import LT.box as B

#from matplotlib import style
#style.use('seaborn-v0_8-poster')

import deut_analysis_tools as dat

#%%
pm580 = dat.DEUT_DATA(branches_sel=dat.branches_types('kins'),
                      Runs=20873)
pm580_SIMC = dat.DEUT_SIMC(kin_study='deep',setting='pm_580')
pm580_SIMC = pm580_SIMC.Branches['pm_580']

b580 = pm580.Branches
r580 = pm580.Runs

pmS_pm580 = dat.make_plots_SIMC('Pm', pm580_SIMC)
pm_pm580 = dat.make_plots('H.kin.secondary.pmiss', b580, r580)

pm_pm580.make_1Dhistos([(-0.1,2.5),100],apply_norm=True,SIMC=pmS_pm580)

pmS_pm580 = dat.make_plots_SIMC('Em', bS580)
pm_pm580 = dat.make_plots('H.kin.secondary.emiss_nuc', b580, r580)

pm_pm580.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,SIMC=pmS_pm580)

empmS_pm580 = dat.make_plots_SIMC(['Pm','Em'],bS580)
empm_pm580 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.emiss_nuc'],b580,r580)

h = empm_pm580.make_2Dhistos([(-0.1,2.5),100,(-0.1,1.0),100],
                             fig_title='Emiss vs. Pmiss',
                             apply_norm=True,plot=True)

hS = empmS_pm580.make_2Dhistos([(-0.1,2.5),100,(-0.1,1.0),100],
                               fig_title='Emiss vs. Pmiss SIMC',
                               apply_weights=True,plot=True)
#%%
pm800 = dat.DEUT_DATA(branches_sel=['kins'],
                      Runs=20886)
pm800_SIMC = dat.DEUT_SIMC(kin_study='deep',setting='pm_800')
pm800_SIMC = pm800_SIMC.Branches['pm_800']

b800 = pm800.Branches
r800 = pm800.Runs
bS800 = pm800_SIMC

pmS_pm800 = dat.make_plots_SIMC('Pm', bS800)
pm_pm800 = dat.make_plots('H.kin.secondary.pmiss', b800, r800)

pm_pm800.make_1Dhistos([(-0.1,2.5),100],apply_norm=True,SIMC=pmS_pm800)

pmS_pm800 = dat.make_plots_SIMC('Em', bS800)
pm_pm800 = dat.make_plots('H.kin.secondary.emiss_nuc', b800, r800)

pm_pm800.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,SIMC=pmS_pm800)

empmS_pm800 = dat.make_plots_SIMC(['Pm','Em'],bS800)
empm_pm800 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.emiss_nuc'],b800,r800)

h = empm_pm800.make_2Dhistos([(-0.1,2.5),100,(-0.1,1.0),100],
                             fig_title='Emiss vs. Pmiss',
                             apply_norm=True,plot=True)

hS = empmS_pm800.make_2Dhistos([(-0.1,2.5),100,(-0.1,1.0),100],
                               fig_title='Emiss vs. Pmiss SIMC',
                               apply_weights=True,plot=True)
#%%
pm900 = dat.DEUT_DATA(branches_sel=['kins'],
                      Runs=20958)
pm900_SIMC = dat.DEUT_SIMC(kin_study='deep',setting='pm_900')
pm900_SIMC = pm900_SIMC.Branches['pm_900']

b900 = pm900.Branches
r900 = pm900.Runs
bS900 = pm900_SIMC

pmS_pm900 = dat.make_plots_SIMC('Pm', bS900)
pm_pm900 = dat.make_plots('H.kin.secondary.pmiss', b900, r900)

pm_pm900.make_1Dhistos([(-0.1,2.5),100],apply_norm=True,SIMC=pmS_pm900)

pmS_pm900 = dat.make_plots_SIMC('Em', bS900)
pm_pm900 = dat.make_plots('H.kin.secondary.emiss_nuc', b900, r900)

pm_pm900.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,SIMC=pmS_pm900)

empmS_pm900 = dat.make_plots_SIMC(['Pm','Em'],bS900)
empm_pm900 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.emiss_nuc'],b900,r900)

h = empm_pm900.make_2Dhistos([(-0.1,2.5),100,(-0.1,1.0),100],
                             fig_title='Emiss vs. Pmiss',
                             apply_norm=True,plot=True)

hS = empmS_pm900.make_2Dhistos([(-0.1,2.5),100,(-0.1,1.0),100],
                               fig_title='Emiss vs. Pmiss SIMC',
                               apply_weights=True,plot=True)


#%%
b120 = dat.pm120.Branches
r120 = dat.pm120.Runs
bS120 = dat.pm120_SIMC

rad2deg=180/np.pi
for run in r120:
    b120[run]['H.kin.secondary.th_bq'] =\
        dat.pm120.Branches[run]['H.kin.secondary.th_bq']*rad2deg
bS120['theta_rq'] = dat.pm120_SIMC['theta_rq']*rad2deg

#Pmiss
pmS_pm120 = dat.make_plots_SIMC('Pm', bS120)
pm_pm120 = dat.make_plots('H.kin.secondary.pmiss', b120, r120)

pm_pm120.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,SIMC=pmS_pm120)

#Emiss
pmS_pm120 = dat.make_plots_SIMC('Em', bS120)
pm_pm120 = dat.make_plots('H.kin.secondary.emiss_nuc', b120, 20871)

h = pm_pm120.make_1Dhistos([(-0.1,1.0),100],apply_norm=True,plot=False)
h = h['H.kin.secondary.emiss_nuc'][20871]

hS = pmS_pm120.make_1Dhistos([(-0.1,1.0),100],apply_weights=True,plot=False)
hS = hS['Em']

#W
WS_pm120 = dat.make_plots_SIMC('W', bS120)
W_pm120 = dat.make_plots('P.kin.primary.W', b120, 20871)

h = W_pm120.make_1Dhistos([(0.3,2.1),100],apply_norm=True,plot=False)
h = h['P.kin.primary.W'][20871]

hS = WS_pm120.make_1Dhistos([(0.3,2.1),100],apply_weights=True,plot=False)
hS = hS['W']

#Emiss vs Pmiss
empmS_pm120 = dat.make_plots_SIMC(['Pm','Em'],bS120)
empm_pm120 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.emiss_nuc'],b120,r120)

h = empm_pm120.make_2Dhistos([(-0.1,1.0),100,(-0.1,1.0),100],
                             fig_title='Emiss vs. Pmiss',
                             apply_norm=False,plot=True)

hS = empmS_pm120.make_2Dhistos([(-0.1,1.0),100,(-0.1,1.0),100],
                               fig_title='Emiss vs. Pmiss SIMC',
                               apply_weights=True,plot=False)

#%%
f = open('YieldSummary.txt','w')
Em_max = [0.25,0.15,0.1,0.05,0.04,0.035]

for MAX in Em_max:
    nevt = []
    f.write(f'Pm(-0.1,{MAX})\n')
    for run in [21075]:
        hx = h[run].project_x((-0.1,MAX))
        nevt.append(hx.sum()[0])
        f.write(f'Nevt [{run}]= {hx.sum()}\n')

    hSx = hS.project_x((-0.1,MAX))    
    f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

    ratio = nevt/hSx.sum()[0]
    f.write(f'{ratio}\n\n')
    
    B.plot_exp(h[21075].x_bin_center, y, kwargs)
print()
print('Pm(-100,150)')
f.write('Pm(-100,150)\n')
nevt = []
for run in [20871,20872]:
    hx = h[run].project_x((-0.1,0.15))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x((-0.1,0.15))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio150 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio150.append(rat)

print(ratio150)   
f.write(f'{ratio150}\n\n')

print()
print('Pm(-100,100)')
f.write('Pm(-100,100)\n')
nevt = []
for run in [20871,20872]:
    hx = h[run].project_x((-0.1,0.1))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x((-0.1,0.1))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio100 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio100.append(rat)

print(ratio100) 
f.write(f'{ratio100}\n\n')

print()
print('Pm(-100,50)')
f.write('Pm(-100,50)\n')
nevt = []
for run in [20871,20872]:
    hx = h[run].project_x((-0.1,0.05))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x((-0.1,0.05))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio50 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio50.append(rat)

print(ratio50)
f.write(f'{ratio50}\n\n')

print()
print('Pm(-100,40)')
f.write('Pm(-100,40)\n')
nevt = []
for run in [20871,20872]:
    hx = h[run].project_x((-0.1,0.04))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x((-0.1,0.04))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio40 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio40.append(rat)

print(ratio40) 
f.write(f'{ratio40}\n\n')
 
print()
print('Pm(-100,35)')
f.write('Pm(-100,35)\n')
nevt = []
for run in [20871,20872]:
    hx = h[run].project_x((-0.1,0.035))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x((-0.1,0.035))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio35 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio35.append(rat)

print(ratio35)
f.write(f'{ratio35}\n\n')

print()
print('Em(-100,100)')
f.write('Em(-100,100)\n')
nevt = []
for run in [20871,20872]:
    hy = h[run].project_y((-0.1,0.1))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')

hSy = hS.project_y((-0.1,0.1))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio100 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio100.append(rat)

print(ratio100) 
f.write(f'{ratio100}\n\n')\
    
print()
print('Em(-100,80)')
f.write('Em(-100,80)\n')
nevt = []
for run in [20871,20872]:
    hy = h[run].project_y((-0.1,0.08))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')
    
hSy = hS.project_y((-0.1,0.08))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio80 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio80.append(rat)

print(ratio80)
f.write(f'{ratio80}\n\n')

print()
print('Em(-100,60)')
f.write('Em(-100,60)\n')
nevt = []
for run in [20871,20872]:
    hy = h[run].project_y((-0.1,0.06))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')
    
hSy = hS.project_y((-0.1,0.06))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio60 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio60.append(rat)

print(ratio60)
f.write(f'{ratio60}\n\n')

print()
print('Em(-100,50)')
f.write('Em(-100,50)\n')
nevt = []
for run in [20871,20872]:
    hy = h[run].project_y((-0.1,0.05))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')
    
hSy = hS.project_y((-0.1,0.05))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio50 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio50.append(rat)

print(ratio50)
f.write(f'{ratio50}\n\n')

print()
print('Em(-100,40)')
f.write('Em(-100,40)\n')
nevt = []
for run in [20871,20872]:
    hy = h[run].project_y((-0.1,0.04))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')

hSy = hS.project_y((-0.1,0.04))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio40 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio40.append(rat)

print(ratio40)
f.write(f'{ratio40}\n\n')
  
print()
print('Em(-100,35)')
f.write('Em(-100,35)\n\n')
nevt = []
for run in [20871,20872]:
    hy = h[run].project_y((-0.1,0.035))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')

hSy = hS.project_y((-0.1,0.035))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio35 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio35.append(rat)

print(ratio35)
f.write(f'{ratio35}')

f.close()
#%%
#Theta_nq
tS_pm120 = dat.make_plots_SIMC('theta_rq', bS120)
t_pm120 = dat.make_plots('H.kin.secondary.th_bq', b120, r120)

t_pm120.make_1Dhistos([(-0.1,3.1),100],apply_norm=False)
tS_pm120.make_1Dhistos([(-0.1,3.1),100],apply_weights=True)
t_pm120.make_1Dhistos([(-0.1,3.1),100],apply_norm=True,SIMC=tS_pm120)

t_pm120.make_1Dhistos([(-0.1,180),100],apply_norm=True,SIMC=tS_pm120)
tS_pm120.make_1Dhistos([(-0.1,180),100],apply_weights=False)

run = 20871

#Theta_nq vs Pmiss
tpmS_pm120 = dat.make_plots_SIMC(['Pm','theta_rq'],bS120)
tpm_pm120 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.th_bq'],b120,run)

h = tpm_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,180.),100],
                             fig_title=r'$\theta_{nq}$ vs. Pmiss',
                             apply_norm=True)

hS = tpmS_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,180),100],
                               fig_title=r'$\theta_{nq}$ vs. Pmiss SIMC',
                               apply_weights=True)

#Theta_nq vs Pmiss, Emiss cut: (-0.1,0.25)
tpmS_pm120 = dat.make_plots_SIMC(['Pm','theta_rq'],bS120,
                                 custom_cuts={'Em':(-0.1,0.25)})
tpm_pm120 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.th_bq'],b120,r120,
                         custom_cuts={'H.kin.secondary.emiss_nuc':(-0.1,0.25)})

h = tpm_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,3.1),100],
                             fig_title=r'$\theta_{nq}$ vs. Pmiss',
                             apply_norm=True)

hS = tpmS_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,3.1),100],
                               fig_title=r'$\theta_{nq}$ vs. Pmiss SIMC',
                               apply_weights=True)

#Theta_nq vs Pmiss, Emiss cut: (-0.1,0.15)
tpmS_pm120 = dat.make_plots_SIMC(['Pm','theta_rq'],bS120,
                                 custom_cuts={'Em':(-0.1,0.15)})
tpm_pm120 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.th_bq'],b120,r120,
                         custom_cuts={'H.kin.secondary.emiss_nuc':(-0.1,0.15)})

h = tpm_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,3.1),100],
                             fig_title=r'$\theta_{nq}$ vs. Pmiss',
                             apply_norm=True)

hS = tpmS_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,3.1),100],
                               fig_title=r'$\theta_{nq}$ vs. Pmiss SIMC',
                               apply_weights=True)

#Theta_nq vs Pmiss, Emiss cut: (-0.1,0.10)
tpmS_pm120 = dat.make_plots_SIMC(['Pm','theta_rq'],bS120,
                                 custom_cuts={'Em':(-0.1,0.1)})
tpm_pm120 = dat.make_plots(['H.kin.secondary.pmiss',
                             'H.kin.secondary.th_bq'],b120,r120,
                         custom_cuts={'H.kin.secondary.emiss_nuc':(-0.1,0.1)})

h = tpm_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,180),100],
                             fig_title=r'$\theta_{nq}$ vs. Pmiss',
                             apply_norm=True)

hS = tpmS_pm120.make_2Dhistos([(-0.1,1.0),100,(0.,180),100],
                               fig_title=r'$\theta_{nq}$ vs. Pmiss SIMC',
                               apply_weights=True)

#%%
f = open('AngleSummary.txt','w')
nevt = []
print('Pm(0.,3.1)')
f.write('Pm(0.,3.1)\n')
for run in r120:
    hx = h[run].project_x(range=(0.,3.1))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,3.1))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio180 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio180.append(rat)

print(ratio180)
f.write(f'{ratio180}\n\n')

nevt = []
print('Pm(0.,2.6)')
f.write('Pm(0.,2.6)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,2.6))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,2.6))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,2.1)')
f.write('Pm(0.,2.1)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,2.1))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,2.1))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,1.6)')
f.write('Pm(0.,1.6)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,1.6))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,1.6))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,1.1)')
f.write('Pm(0.,1.1)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,1.1))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,1.1))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,0.6)')
f.write('Pm(0.,0.6)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,0.6))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,0.6))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,0.5)')
f.write('Pm(0.,0.5)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,0.5))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,0.5))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,0.4)')
f.write('Pm(0.,0.4)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,0.4))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,0.4))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

nevt = []
print('Pm(0.,0.3)')
f.write('Pm(0.,0.3)\n')
for run in [21074,21075]:
    hx = h[run].project_x(range=(0.,0.3))
    nevt.append(hx.sum()[0])
    print(f'Nevt [{run}]=', hx.sum())
    f.write(f'Nevt [{run}]= {hx.sum()}\n')

hSx = hS.project_x(range=(0.,0.3))    
print('Nevt [SIMC]=', hSx.sum())
f.write(f'Nevt [SIMC]= {hSx.sum()}\n')

ratio130 = []
for SUM in nevt:
    rat = SUM/hSx.sum()[0]
    ratio130.append(rat)

print(ratio130)
f.write(f'{ratio130}\n\n')

print()
print('theta_nq(-0.1,100)')
f.write('theta_nq(-0.1,100)\n')
nevt = []
for run in [21074,21075]:
    hy = h[run].project_y((-0.1,0.1))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')

hSy = hS.project_y((-0.1,0.1))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio100 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio100.append(rat)

print(ratio100) 
f.write(f'{ratio100}\n\n')
    
print()
print('theta_nq(-100,80)')
f.write('theta_nq(-100,80)\n')
nevt = []
for run in [21074,21075]:
    hy = h[run].project_y((-0.1,0.08))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')
    
hSy = hS.project_y((-0.1,0.08))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio80 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio80.append(rat)

print(ratio80)
f.write(f'{ratio80}\n\n')

print()
print('theta_nq(-100,60)')
f.write('theta_nq(-100,60)\n')
nevt = []
for run in [21074,21075]:
    hy = h[run].project_y((-0.1,0.06))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')
    
hSy = hS.project_y((-0.1,0.06))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio60 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio60.append(rat)

print(ratio60)
f.write(f'{ratio60}\n\n')

print()
print('theta_nq(-100,50)')
f.write('theta_nq(-100,50)\n')
nevt = []
for run in [21074,21075]:
    hy = h[run].project_y((-0.1,0.05))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')
    
hSy = hS.project_y((-0.1,0.05))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio50 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio50.append(rat)

print(ratio50)
f.write(f'{ratio50}\n\n')

print()
print('theta_nq(-100,35)')
f.write('theta_nq(-100,35)\n')
nevt = []
for run in [21074,21075]:
    hy = h[run].project_y((-0.1,0.035))
    nevt.append(hy.sum()[0])
    print(f'Nevt [{run}]=', hy.sum())
    f.write(f'Nevt [{run}]= {hy.sum()}\n')

hSy = hS.project_y((-0.1,0.035))    
print('Nevt [SIMC]=', hSy.sum())
f.write(f'Nevt [SIMC]= {hSy.sum()}\n')

ratio35 = []
for SUM in nevt:
    rat = SUM/hSy.sum()[0]
    ratio35.append(rat)

print(ratio35)
f.write(f'{ratio35}\n\n')

f.close()


#%% PLOTTING ==================================================================
#   analysis of 13 pm120 deep runs + pm80 from spring2017

# Em Pm
pm120_empm = make_plots(['H.kin.secondary.pmiss','H.kin.secondary.emiss_nuc'],
                        pm120_branches,pm120_runs)
pm120_empm.make_1Dhistos([(-0.1,1.0),100])

pm120_empm.make_2Dhistos([(-0.1,1.0),100,(-0.1,1.0),100],'Emiss vs. Pmiss')

pm80_empm = make_plots(['H.kin.secondary.pmiss','H.kin.secondary.emiss_nuc'],
                        pm80_branches,pm80_runs)
pm80_empm.make_1Dhistos([(-0.1,1.0),100])

pm80_empm.make_2Dhistos([(-0.1,1.0),100,(-0.1,1.0),100],'Emiss vs. Pmiss')

pm120_SIMC_pm = make_plots_SIMC('Pm',pm120_SIMC_branches)

pm120_SIMC_pm.make_1Dhistos([(-0.1,2.0),100])

pm120_SIMC_em = make_plots_SIMC('Em',pm120_SIMC_branches)

pm120_SIMC_em.make_1Dhistos([(-0.1,1.0),100])

pm120_SIMC_empm = make_plots_SIMC(['Pm','Em'],pm120_SIMC_branches,apply_cuts=False)

pm120_SIMC_empm.make_2Dhistos([(-0.1,2.0),100,(-0.1,1.0),100])

# Mrecoil

pm120_mrec_50 = make_plots('H.kin.secondary.Mrecoil',
                        pm120_branches,pm120_runs,
                        custom_cuts={'H.kin.secondary.pmiss':(0.,0.5)})
pm120_mrec_100 = make_plots('H.kin.secondary.Mrecoil',
                        pm120_branches,pm120_runs,
                        custom_cuts={'H.kin.secondary.pmiss':(0.5,1.0)})
pm120_mrec_150 = make_plots('H.kin.secondary.Mrecoil',
                        pm120_branches,pm120_runs,
                        custom_cuts={'H.kin.secondary.pmiss':(1.0,1.5)})
pm120_mrec_200 = make_plots('H.kin.secondary.Mrecoil',
                        pm120_branches,pm120_runs,
                        custom_cuts={'H.kin.secondary.pmiss':(1.5,2.0)})

pm120_mrec_50.make_1Dhistos([(-0.1,2.5),100])
pm120_mrec_100.make_1Dhistos([(-0.1,2.5),100])
pm120_mrec_150.make_1Dhistos([(-0.1,2.5),100])
pm120_mrec_200.make_1Dhistos([(-0.1,2.5),100])

pm120_W2 = make_plots('P.kin.primary.W',
                        pm120_branches,pm120_runs)
pm120_W2.make_1Dhistos([(0.1,2.5),100])