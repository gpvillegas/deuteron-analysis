############################################################################
Calculating bin centering corrections and momentum distributions
WORKFLOW
############################################################################

Before staring this workflow, you must have a file containing the 2D average 
kinematic histograms, I use a pickle (*.pcl) file for this which contains a 
dictionary with the desired histograms in it.

Run the following programs in this order:

1. 'calc_avg_kin_py': this code gets the average kinematic histograms and 
extracts all relevant information from them. The result is a text file with the
relevant kinematics for each bin, the Mott xsec and deForest xsec get 
calculated here. This code has to be run separately for fsi and pwia files.
    - Usage: when prompted enter
        pm_setting: pm_120, pm_580, pm_800, or pm_900
        model: enter whole model name, for example: 'jmlpwia_norad' or 
        'Paris_fsi_norad'

2. 'average_xsec/calc_Xsec.py': this code reads the file with the average xsec 
histograms and extracts all relevant information in the same fashion as 
'calc_avg_kin.py'.
The result is also a text file with the pwia and fsi data xsec for a given
pm setting and set.
    - Usage: when prompted enter
    pm_setting: pm_120, pm_580, pm_800, or pm_900
    model: enter just the name of the model, for example: 'jml' or 'Paris'
    data_set: enter data set number as 'set_1' or 'set_2', etc...
    
3. 'jml_theory_xsec/calc_theory_Xsec.py': this code runs the Laget model xsec 
to calculate the xsec at the average kinematics for each bin.
    - Usage: when prompted enter
    pm_setting: pm_120, pm_580, pm_800, or pm_900
    
4. 'bin_centering/calc_bc_corr.py': this code calculates bin centering 
corrections and applies them to the fsi norad data. The reduced xsec (momentum
distribution) is also calculated here.
    - Usage: when prompted enter
    pm_setting: pm_120, pm_580, pm_800, or pm_900
    model: enter just the name of the model, for example: 'jml' or 'Paris'
    data_set: enter data set number as 'set_1' or 'set_2', etc...
    
5. 'bin_centering/check_bc_corr.py': this code plots the data xsec, theory xsec
and bin centering corr to check the correction.
 

