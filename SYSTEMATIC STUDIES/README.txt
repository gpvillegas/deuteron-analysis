#######################################
# SYSTEMATICS STUDIES FOR DEUTERON 23 #
#######################################
 
CUT SYSTEMATICS
code location: cut_systematics/cut_systematics.py
this code uses the class R_Barlow to perform the Barlow test on different
analysis cuts.

1. Nominal Cuts
 Window cut: hms_delta; xmin = -8.0, xmax = 8.0
 Window cut: shms_delta; xmin = -10.0, xmax = 22.0
 Window cut: shms_calPID; xmin = 0.7, xmax = 1.3
 Window cut: Em_cut; xmin = -0.05, xmax = 0.05
 Window cut: Q2_cut; xmin = 4.0, xmax = 5.0
 Value cut: noEDTM; value = 0.0
 HMS collimator cut ON
 BCM1_cut
 ztar_diff_cut
 epCoinTime_cut

2. Varying Cuts *DONE
    *1 - HCOLL: HMS Collimator Cut - vary scale
        A) 1.1, 1.0 (Nominal),  B) 0.9,  C) 0.5
    *2 - shms_calPID: SHMS Calorimeter Cut
        ON (Nominal) / OFF  
    *3 - Em_cut: Missing Energy Cut
        A) (-0.08,0.08), B) (-0.07,0.07), C) (-0.06,0.06), (-0.05,0.05) (Nominal),
        D) (-0.04,0.04), E) (-0.03,0.03)
    *4 - ZT: Ztar Difference (HMS Ztar - SHMS Ztar) cut - cut around mean value
        A) (-2.5,2.5), (-2,2) (Nominal), B) (-1.5,1.5), C) (-1,1)
    *5 - epCoinTime_cut: e,p Coincidence Time Cut
        ON (Nominal) / OFF
    *6 - noEDTM: EDTM TDC Signal Cut 
        ON (Nominal) / OFF
    *7 - BCM1_cut: BCM1 Current Cut 
        ON (Nominal) / OFF
        
Study on Cut systematics: Why are R Barlow positive for ON/OFF cut tests?
Trying using ranges instead:
    1 - shms_calPID:
        (0.7,1.3) (nominal), A) (0.8,1.2), B) (0.9,1.1), C) (0.95,1.05)
    2 - epCoinTime:
        +/- 2 (nominal), A) 1, B) 0.5, C) 0.25 

NORMALIZATON SYSTEMATICS: 

code location: norm_systematics/norm_systematics.py

This code calculates the systematic uncertainty due
to each normalization factor on the cross section. 3 errors are calculated: 
a run-by-run error added in quadrature at the end, an error on the regular 
average, and an error on the weighted average of a set per setting.

The normalization factors are: 
    1. f_tgtBoil --> Target Boiling Correction
    2. f_pTr --> Proton Transmission Correction
    3. f_pTrk --> SHMS Tracking Efficiency
    4. f_hTrk --> HMS Tracking Efficiency
    5. f_tLT --> EDTM Total Live Time
    6. f_Qtot --> Total accumulated charge
    
The total systematic uncertainty on the cross section due to these factors 
would be:
    dsig^2 = \sum[i=6 (dsig/df_i)^2 * df_i^2 ]
where sig is the experimental cross section (normalized) 
    and dsig/df_i = sig/f_i

So we end up with a simpler form:
    dsig^2 = sig^2 * \sum[i=6 (df_i/f_i)^2 ]

Explicitly this is:
    dsig^2 = sig^2 *( [df/f]_tgtBoil^2 + [df/f]_pTr^2 + [df/f]_pTrk^2 + \
         [df/f]_hTrk^2 + [df/f]_tLT^2 + [df/f]_Qtot^2)

So the systematic error squared is the sum of the squares of the relative errors 
of all efficiency and correction factors times the cross section squared.

The uncertainties of each factor are located:
    1. tgtBoil: written in DB per run
    2. pTr: one overall = 0.953 +/- 0.011
    3. pTrk: written in DB per run
    4. hTrk: written in DB per run
    5. tLT: written in DB per run
    6. Qtot: TBD, one overall, estimate dI/I ~ 2%   