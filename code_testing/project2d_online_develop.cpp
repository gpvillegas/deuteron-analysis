//#ifndef PROJECT_2D_H
//#define PROJECT_2D_H

/*
  Author: C. Yero
  Date: Feb 15, 2023 
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../UTILS/parse_utils.h" //useful C++ string parsing utilities

using namespace std;



// helper function : check if a string only contain numebers or not
int is_digits(string& str)
{
    for (char ch : str) {
        int v = ch; // ASCII Val converted
        if (!(ch >= 48 && ch <= 57)) {   // The ASCII value of '0' is 48 , '9' is 57.
            
	  return 0;  // characters in string have letters
        }
	
	
    }
 
    return 1; // no charactes have letter (therefore true number)
}


// helper function : combined 2d histos for a given run range based on a .csv file,
// and looks for the histos in a pre-defined path
TH2F* combine_2dhistos(int run_min=0, int run_max=99999, TString hist2d_name="randSub_plots/H_Pm_vs_thrq_rand_sub", bool apply_corrections=false)
{

  // cout << Form("Will skip "" runs", skip_len) << endl;
  
  //experiment online file to read run numbers (assumes .csv format, and runs to be 1st column)
  ifstream file("UTILS_DEUT/runlist/deut-2023_runlist.csv"); 
  string line;
  string token;
  int run_num;

  // define generic root file name to be read 
  string root_file_path; // generic path to where root files are written

  // define tfile to be used to read .root
  TFile *data_file = NULL;

  // define generic histogram to be combined
  TH2F *myhist2d = 0;
  TH2F *myhist2d_total = 0;

  TH2F *myhist2d_corr = 0;
  TH2F *myhist2d_corr_total = 0;

  //define total charge and track inefficiencies for correction later on
  // (charge should be summed separately and scaled by total counts after combinin histos)
  // inefficiencies should be applied on a run by run basis
  double Q, hms_trk_eff, hms_trk_eff_err, shms_trk_eff, shms_trk_eff_err, total_live_time, total_live_time_err;
  double Q_tot = 0; // charge counter
  
  int cnt=0; // good run counter
  
    // read line by line
  while (getline(file, line)) {

    // get token (parsed line string up to the 1st comma)
    token = line.substr(0, line.find(','));

    // check if token is a number
    if( is_digits(token) ){


      // convert token to int
      run_num = stoi(token);

      
      // define generic rootfile name to path
      root_file_path =Form("./DEUT_OUTPUT/ROOT/deut_prod_LD2_deep_%d_-1_histos.root", run_num);

      // if file does not exist continue
      if(gSystem->AccessPathName(root_file_path.c_str())) continue;

      // skip runs outside range
      if( (run_num<run_min) or (run_num>run_max) ) continue;
      
      
      cout << "---> reading ROOTfile: " << root_file_path.c_str() << endl;

      TString report_output_path = Form("./DEUT_OUTPUT/REPORT/deut_prod_LD2_deep_report_%d_-1.txt", run_num);
      // read the relevant inefficiencies and charge (FOR ONLINE BCM4A column was reading valued from BCM4C, which are no much different ~2-3 mC)
      Q = stod(split(FindString("BCM4C_Charge [mC]", report_output_path.Data())[0])[1]);
      Q_tot = Q_tot + Q;
      
      hms_trk_eff = stod(split(split(FindString("hms_had_track_eff", report_output_path.Data())[0], ':')[1], '+')[0]);
      hms_trk_eff_err = stod(split(split(split(FindString("hms_had_track_eff", report_output_path.Data())[0], ':')[1], '+')[1], '-')[1]);

      shms_trk_eff = stod(split(split(FindString("shms_elec_track_eff", report_output_path.Data())[0], ':')[1], '+')[0]);
      shms_trk_eff_err = stod(split(split(split(FindString("shms_elec_track_eff", report_output_path.Data())[0], ':')[1], '+')[1], '-')[1]);
      
      total_live_time = stod(split(split(FindString("T6_tLT", report_output_path.Data())[0], ':')[1], '+')[0]);
      total_live_time_err = stod(split(split(split(FindString("T6_tLT", report_output_path.Data())[0], ':')[1], '+')[1], '-')[1]);

      cout << "------------------" << endl;
      cout << Form("run number: %d", run_num) << endl;
      cout << "------------------" << endl;
      cout << Form("Q [mC] = %.3f, Q_tot = %.3f", Q, Q_tot) << endl;
      cout << Form("htrack_eff = %.3f +/- %.3f ",  hms_trk_eff,  hms_trk_eff_err) << endl;
      cout << Form("ptrack_eff = %.3f +/- %.3f ",  shms_trk_eff,  shms_trk_eff_err) << endl;
      cout << Form("tLT = %.3f +/- %.3f ",  total_live_time,  total_live_time_err) << endl;



      // for each root file, get the desired histogram to be combined
      data_file =  new TFile(root_file_path.c_str(), "READ");

      data_file->GetObject(Form("%s", hist2d_name.Data()), myhist2d);
      
      // clone histo to be corrected (so we may have before/after corrections)
      myhist2d_corr = (TH2F*)myhist2d->Clone();


      // apply inefficiency corrections (run-by-run)
      double eff_tot = 1./(hms_trk_eff*shms_trk_eff*total_live_time);
      myhist2d_corr->Scale(eff_tot);
      
      // only for 1st run, clone histogram to the total
      if(cnt==0) { 
	myhist2d_total = (TH2F*)myhist2d->Clone(hist2d_name.Data()); 
	myhist2d_corr_total = (TH2F*)myhist2d_corr->Clone(hist2d_name.Data()); 
	
      }

      cout << "passed L1" << endl;
      // add subseqquent histos
      if(cnt>0) { 
	myhist2d_total->Add(myhist2d); 
	myhist2d_corr_total->Add(myhist2d_corr); 
      }
      
      // increment counter
      cnt++;
      
    } // end check if token is number
    
  } // end readlines

  cout << "passed L2" << endl;

  // apply charge normalization to total histogram
  myhist2d_corr_total->Scale(1./Q_tot);

  if(apply_corrections){
    return myhist2d_corr_total;  // efficiency-corrected, charge-normalized yield
  }
  else{
    return myhist2d_total; 
  }

}


TH2F* get_simc_2d_histos(TString setting="pm120", TString model="fsi", TString hist_type="rad_corr_ratio"){

  // Brief: this function gets specific histogram from SIMC for
  // radiative corrections and phase space, for a given deuteron kin
  // the arguments are:
  // setting="pm120", "pm580", "pm800", "pm900"
  // model="fsi" or "pwia"
  // hist_type="norad", "rad", "rad_corr_ratio" , "phase_space", 

  cout << "calling get_simc_2d_histos()" << endl;
  
  // simc no_rad / rad 2D histo (FSI)
  TH2F *dummy = 0;
  TH2F *H2_Pm_vs_thrq_simc_fsi_rad = 0;
  TH2F *H2_Pm_vs_thrq_simc_fsi_norad = 0;
  TH2F *H2_Pm_vs_thrq_simc_fsi_ps = 0; //phase space
  TH2F *H2_Pm_vs_thrq_simc_fsi_ratio = 0; // for nonrad/rad ratio 

  // simc PWIA 
  TH2F *H2_Pm_vs_thrq_simc_pwia_norad = 0;
  TH2F *H2_Pm_vs_thrq_simc_pwia_ps = 0; //phase space

  
  // define generic rootfile name to path (ln -sf ../../hallc_simulations/worksim/analyzed/pass1 SIMC)
  TString root_file_path_fsi_rad = Form("SIMC/d2_%s_jmlfsi_rad_analyzed.root", setting.Data());
  TString root_file_path_fsi_norad = Form("SIMC/d2_%s_jmlfsi_norad_analyzed.root", setting.Data());

  // SIMC PWIA
  TString root_file_path_pwia_norad = Form("SIMC/d2_%s_jmlpwia_norad_analyzed.root", setting.Data());

  
  // if file does not exist, EXIT
  if(gSystem->AccessPathName(root_file_path_fsi_rad.Data()) || gSystem->AccessPathName(root_file_path_fsi_norad.Data()) || gSystem->AccessPathName(root_file_path_pwia_norad.Data()) ) { gSystem->Exit(kTRUE);}


  // for each root file, get the desired histogram to be combined
  TFile *simc_file_fsi_rad =  new TFile(root_file_path_fsi_rad.Data(), "READ");
  TFile *simc_file_fsi_norad =  new TFile(root_file_path_fsi_norad.Data(), "READ");

  TFile *simc_file_pwia_norad =  new TFile(root_file_path_pwia_norad.Data(), "READ");

 
  cout << Form("Read data file . . . %s",root_file_path_fsi_rad.Data()) << endl;
  
    // Retrieve FSI histograms
  simc_file_fsi_rad->cd();
  simc_file_fsi_rad->GetObject("kin_plots/H_Pm_vs_thrq", H2_Pm_vs_thrq_simc_fsi_rad);  //radiative

  cout << Form("Read data file . . . %s",root_file_path_fsi_norad.Data()) << endl;
  simc_file_fsi_norad->cd();
  simc_file_fsi_norad->GetObject("kin_plots/H_Pm_vs_thrq", H2_Pm_vs_thrq_simc_fsi_norad); //non-radiative
  simc_file_fsi_norad->GetObject("kin_plots/H_Pm_vs_thrq_ps", H2_Pm_vs_thrq_simc_fsi_ps); //non-radiative phase-space


  // Retrieve PWIA histograms
  simc_file_pwia_norad->cd();

  simc_file_pwia_norad->GetObject("kin_plots/H_Pm_vs_thrq", H2_Pm_vs_thrq_simc_pwia_norad); //non-radiative
  simc_file_pwia_norad->GetObject("kin_plots/H_Pm_vs_thrq_ps", H2_Pm_vs_thrq_simc_pwia_ps); //non-radiative phase-space

  
  int xnbins = H2_Pm_vs_thrq_simc_fsi_norad->GetXaxis()->GetNbins();
  float xmin = H2_Pm_vs_thrq_simc_fsi_norad->GetXaxis()->GetXmin();
  float xmax = H2_Pm_vs_thrq_simc_fsi_norad->GetXaxis()->GetXmax();

  int ynbins = H2_Pm_vs_thrq_simc_fsi_norad->GetYaxis()->GetNbins();
  float ymin = H2_Pm_vs_thrq_simc_fsi_norad->GetYaxis()->GetXmin();
  float ymax = H2_Pm_vs_thrq_simc_fsi_norad->GetYaxis()->GetXmax();

  H2_Pm_vs_thrq_simc_fsi_ratio = new TH2F("H2_Pm_vs_thrq_simc_fsi_ratio", "SIMC Y_{norad}/Y_{rad}; Y_{norad}/Y_{rad}; #theta_{rq} [deg] ", xnbins, xmin, xmax, ynbins, ymin, ymax);
  H2_Pm_vs_thrq_simc_fsi_ratio->Divide(H2_Pm_vs_thrq_simc_fsi_norad, H2_Pm_vs_thrq_simc_fsi_rad);

  
  if(hist_type=="norad" && model=="fsi"){
    return H2_Pm_vs_thrq_simc_fsi_norad;
  }
    
  else if(hist_type=="norad" && model=="pwia"){
    return H2_Pm_vs_thrq_simc_pwia_norad;
  }

  else if(hist_type=="rad" && model=="fsi"){
    return H2_Pm_vs_thrq_simc_fsi_rad;
  }
    
  else if(hist_type=="rad_corr_ratio" && model=="fsi"){
    return H2_Pm_vs_thrq_simc_fsi_ratio;
  }

  else if(hist_type=="phase_space" && model=="fsi"){
    return  H2_Pm_vs_thrq_simc_fsi_ps;
  }

  else if(hist_type=="phase_space" && model=="pwia"){
    return  H2_Pm_vs_thrq_simc_pwia_ps;
  }
  
  else{
    return dummy;
  }
  
  

}


TH2F* get_comm18_data_2d_histos(TString setting="pm80", TString dataset="dataset1"){

  // Brief: this function gets specific histogram from DATA for
  // from comissioning deuteron experiment on 2018
  // the arguments are:
  // setting="pm80", "pm580", "pm750"
  // dataset="dataset1", "dataset2", "dataset3"
  

  cout << "calling get_comm18_data_2d_histos()" << endl;
  
  // simc no_rad / rad 2D histo (FSI)
  TH2F *dummy = 0;
  TH2F *H_Pm_vs_thrq_data = 0;

  
  // define generic rootfile name to path 
  TString root_file_path= Form("./commissioning_data/%s/Xsec_%s_lagetfsi_%s.root", setting.Data(), setting.Data(), dataset.Data());

  // if file does not exist, EXIT
  if(gSystem->AccessPathName(root_file_path.Data()) ) { gSystem->Exit(kTRUE);}


  // for each root file, get the desired histogram to be combined
  TFile *data_file =  new TFile(root_file_path.Data(), "READ");


 
  cout << Form("Read data file . . . %s",root_file_path.Data()) << endl;
  
  // Retrieve data histograms
  data_file->cd();
  data_file->GetObject("H_data2DXsec", H_Pm_vs_thrq_data);  


  return H_Pm_vs_thrq_data;

}


void project2d_deut( TH2F *hist2d=0, TH2F *hist2d_corr=0, TString setting="", Bool_t display_plots=0, Bool_t apply_radiative_corr=0 ){

  cout << "calling project2d_deut" << endl;
  
  if(display_plots) {
    //avoid interactive display{
    gROOT->SetBatch(kTRUE);
  }

  /*
    Brief: Projects 2d histograms onto 1d slices in X or Y
    and plots it on a canvas subplot, and also plots relative
    errors on a separate canvas subplot

    For now is specific for deuteron, but can easily be modified for other use
   */

  // set base filename
  TString basename=Form("deut_stats_monitoring_setting_%s_", setting.Data());
  
  // Set basefilename to save projections to root file
  TString ofile="DEUT_OUTPUT/ROOT/" + basename + "output.root";
  
  TFile *fout = new TFile(ofile.Data(),"RECREATE");


  //set output names of projections to be saved to divided canvas
  TString fout_2dHist          = "DEUT_OUTPUT/PDF/" + basename + "2Dhisto.pdf";
  TString fout_projHist        = "DEUT_OUTPUT/PDF/" + basename + "projY.pdf";
  TString fout_projHistErr     = "DEUT_OUTPUT/PDF/" + basename + "projY_relError.pdf";
  TString fout_projsimcRadCorr = "DEUT_OUTPUT/PDF/" + basename + "projY_simcRadCorr.pdf";
  TString fout_projdataRadCorr = "DEUT_OUTPUT/PDF/" + basename + "projY_dataRadCorr.pdf";
  TString fout_projsimcPS      = "DEUT_OUTPUT/PDF/" + basename + "projY_simcPS.pdf";
  TString fout_projdataXsec    = "DEUT_OUTPUT/PDF/" + basename + "projY_dataXsec.pdf";

  // set global title/label sizes
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetLabelSize(.1, "XY");
  gStyle->SetTitleY(1.01); // offset title vertically

  // define external 2d histogrmas
  TH2F *hist2d_Pm_vs_thrq_simc_fsi_norad = 0; //non-radiative fsi simc

  TH2F *hist2d_Pm_vs_thrq_simc_fsi_ratio = 0; //non-radiative / radiative ratio
  TH2F *hist2d_Pm_vs_thrq_simc_fsi_ps = 0;  //phase space
  TH2F *hist2d_Pm_vs_thrq_data_radcorr = 0;  //radiative corr. data
  TH2F *hist2d_Pm_vs_thrq_data_Xsec = 0;  // data cross sections
  TH2F *hist2d_Pm_vs_thrq_simc_fsi_Xsec = 0;  // simc jml fsi cross sections

  TH2F *hist2d_Pm_vs_thrq_simc_pwia_norad = 0; //non-radiative pwia simc
  TH2F *hist2d_Pm_vs_thrq_simc_pwia_ps = 0;  // pwia phase space (should be the same as fsi, since its model independent)
  TH2F *hist2d_Pm_vs_thrq_simc_pwia_Xsec = 0;  // simc jml pwia cross sections

  // define external 2d histograms (from commissioning 2018 deuteron data)
  TH2F *hist2d_Pm_vs_thrq_data_Xsec_pm80;
  TH2F *hist2d_Pm_vs_thrq_data_Xsec_pm580;
  TH2F *hist2d_Pm_vs_thrq_data_Xsec_pm750;
  
  
  cout << "Retrieve histograms get_simc_2d_histos() func. . . " << endl;

  hist2d_Pm_vs_thrq_simc_fsi_norad     = get_simc_2d_histos(setting.Data(), "fsi", "norad");
  hist2d_Pm_vs_thrq_simc_fsi_ratio     = get_simc_2d_histos(setting.Data(), "fsi", "rad_corr_ratio");
  hist2d_Pm_vs_thrq_simc_fsi_ps        = get_simc_2d_histos(setting.Data(), "fsi", "phase_space");

  hist2d_Pm_vs_thrq_simc_pwia_norad     = get_simc_2d_histos(setting.Data(), "pwia", "norad");
  hist2d_Pm_vs_thrq_simc_pwia_ps        = get_simc_2d_histos(setting.Data(), "pwia", "phase_space");

  cout << "Retrieve histograms get_comm18_data_2d_histos() func. . . " << endl;

  hist2d_Pm_vs_thrq_data_Xsec_pm80 = get_comm18_data_2d_histos("pm80","dataset1");
  hist2d_Pm_vs_thrq_data_Xsec_pm580 = get_comm18_data_2d_histos("pm580","dataset2");  // dataset with highest stats is chosen
  hist2d_Pm_vs_thrq_data_Xsec_pm750 = get_comm18_data_2d_histos("pm750","dataset1"); // dataset with highest stats is chosen
  
  //------------------------------
  //Apply radiative corrections to charge-norm histo
  //------------------------------
  hist2d_Pm_vs_thrq_data_radcorr = (TH2F*)hist2d_corr->Clone();
  hist2d_Pm_vs_thrq_data_radcorr->Multiply(hist2d_Pm_vs_thrq_simc_fsi_ratio);
  hist2d_Pm_vs_thrq_data_radcorr->SetLabelSize(.05, "XY");

  //----------------------------
  // divide  radiative-corrected data by phase space to get cross data sections
  //----------------------------

  int xnbins = hist2d->GetXaxis()->GetNbins();
  float xmin = hist2d->GetXaxis()->GetXmin();
  float xmax = hist2d->GetXaxis()->GetXmax();

  int ynbins = hist2d->GetYaxis()->GetNbins();
  float ymin = hist2d->GetYaxis()->GetXmin();
  float ymax = hist2d->GetYaxis()->GetXmax();
  
  hist2d_Pm_vs_thrq_data_Xsec= new TH2F("hist2d_Pm_vs_thrq_data_Xsec", "Data Cross Sections (Online); d^{5}#sigma/d#Omega_{e,p}d#omega [#ub sr^{-2} MeV^{-1}]; #theta_{rq} [deg] ", xnbins, xmin, xmax, ynbins, ymin, ymax);
  hist2d_Pm_vs_thrq_data_Xsec->Divide(hist2d_Pm_vs_thrq_data_radcorr, hist2d_Pm_vs_thrq_simc_fsi_ps);
  hist2d_Pm_vs_thrq_data_Xsec->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d_Pm_vs_thrq_data_Xsec->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");
  hist2d_Pm_vs_thrq_data_Xsec->SetTitle("Data Cross Sections (Online)");
  hist2d_Pm_vs_thrq_data_Xsec->SetLabelSize(.03, "XY");


  //----------------------------
  // divide  non-radiative SIMC by phase space to get cross sections
  //----------------------------
  hist2d_Pm_vs_thrq_simc_fsi_Xsec= new TH2F("hist2d_Pm_vs_thrq_simc_fsi_Xsec", "SIMC JML FSI (Paris) Cross Sections (Online); d^{5}#sigma/d#Omega_{e,p}d#omega [#ub sr^{-2} MeV^{-1}]; #theta_{rq} [deg] ", xnbins, xmin, xmax, ynbins, ymin, ymax);
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->Divide(hist2d_Pm_vs_thrq_simc_fsi_norad, hist2d_Pm_vs_thrq_simc_fsi_ps);
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->SetTitle("SIMC JML FSI (Paris) Cross Sections (Online)");
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->SetLabelSize(.03, "XY");

  hist2d_Pm_vs_thrq_simc_pwia_Xsec= new TH2F("hist2d_Pm_vs_thrq_simc_pwia_Xsec", "SIMC JML PWIA (Paris) Cross Sections (Online); d^{5}#sigma/d#Omega_{e,p}d#omega [#ub sr^{-2} MeV^{-1}]; #theta_{rq} [deg] ", xnbins, xmin, xmax, ynbins, ymin, ymax);
  hist2d_Pm_vs_thrq_simc_pwia_Xsec->Divide(hist2d_Pm_vs_thrq_simc_pwia_norad, hist2d_Pm_vs_thrq_simc_pwia_ps);
  hist2d_Pm_vs_thrq_simc_pwia_Xsec->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d_Pm_vs_thrq_simc_pwia_Xsec->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");
  hist2d_Pm_vs_thrq_simc_pwia_Xsec->SetTitle("SIMC JML PWIA (Paris) Cross Sections (Online)");
  hist2d_Pm_vs_thrq_simc_pwia_Xsec->SetLabelSize(.03, "XY");


  

  cout << "Retrieved 2d histos . . . " << endl;
  
 

  hist2d->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");
  hist2d->SetTitle("data (raw counts)");
  cout << "Set title for hist2d . . . " << endl;

  hist2d_Pm_vs_thrq_simc_fsi_ratio->GetYaxis()->SetTitle("P_{m}, Missing Momentum [GeV/c]");
  hist2d_Pm_vs_thrq_simc_fsi_ratio->GetXaxis()->SetTitle("Recoil Angle, #theta_{rq} [deg]");
  hist2d_Pm_vs_thrq_simc_fsi_ratio->SetTitle("SIMC Y_{norad}/Y_{rad}");
  hist2d_Pm_vs_thrq_simc_fsi_ratio->SetLabelSize(.05, "XY");




  // plot 2d histograms
  TCanvas *c0 = new TCanvas("c0", "", 1500,1500);
  c0->Divide(3,2);

  c0->cd(1);
  gPad->Modified(); gPad->Update();
  hist2d->Draw("contz");
 
  c0->cd(2);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_simc_fsi_ratio->Draw("contz");

  c0->cd(3);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_simc_fsi_ps->Draw("contz");

  c0->cd(4);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_data_radcorr->Draw("contz");
  
  c0->cd(5);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_data_Xsec->Draw("contz");

  c0->cd(6);
  gPad->Modified(); gPad->Update();
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->Draw("contz");
  
  fout->cd();
  hist2d->Write();
  hist2d_Pm_vs_thrq_simc_fsi_ratio->Write();
  hist2d_Pm_vs_thrq_simc_fsi_ps->Write();
  hist2d_Pm_vs_thrq_data_radcorr->Write();
  hist2d_Pm_vs_thrq_data_Xsec->Write();
  hist2d_Pm_vs_thrq_simc_fsi_Xsec->Write();

  // --- define variables for calculative/plotting of relative stats. error on Pmiss (need to reset vector per projection bin) ---
  vector<double> y_val;       // this will be set to 0 (as reference)
  vector<double> y_err;       // relative error on pmiss bin counts
  vector<double> x_val;       // this is the central value of x (pmiss)
  vector<double> x_err;       // this is actually width on x-axis 
  
  double pm_counts;
  double relative_err;
  double pm_center;
  double thrq_center;
  double thrq_width;
  
  // ----------------------------

  
  // --- define canvas subplots (x,y) based on number of bins in hist2d ---
  int yc, xc;
  // if projection is along y, should make a canvas square out of xnbins
  float remainder = sqrt(xnbins) - int(sqrt(xnbins));

  if(remainder>0 && remainder<0.5){
    yc = round(sqrt(xnbins));
    xc = round(sqrt(xnbins))+1;
  }
  else{
    yc = round(sqrt(xnbins));
    xc = round(sqrt(xnbins));
  }
  
  // define canvas for projecting 2d histograms
  TCanvas *c1 = new TCanvas("c1", "Data Missing Momenta ProjY (raw counts)", 1400,1400);
  TCanvas *c2 = new TCanvas("c2", "Data Missing Momenta Relative Error (raw counts)", 1400,1400);
  TCanvas *c3 = new TCanvas("c3", "SIMC Radiative Corrections", 1400,1400);
  TCanvas *c4 = new TCanvas("c4", "DATA Radiative Corrected", 1400,1400);
  TCanvas *c5 = new TCanvas("c5", "SIMC Phase Space", 1400,1400);
  TCanvas *c6 = new TCanvas("c6", "DATA Xsec", 1400,1400);

  c1->cd();
  c1->Divide(yc, xc);
  
  c2->cd();
  c2->Divide(yc, xc);
  
  c3->cd();
  c3->Divide(yc, xc);
  
  c4->cd();
  c4->Divide(yc, xc);

  c5->cd();
  c5->Divide(yc, xc);

  c6->cd();
  c6->Divide(yc, xc);
  //----------------------------


  // define 1d projection histos 
  TH1D *H_dataPm_projY = 0;
  TH1D *H_simcPm_projY_ratio = 0; // ratio of norad/rad projected
  TH1D *H_dataPm_projY_radUnCorr = 0; // charge-normalized data before radiative corrections
  TH1D *H_dataPm_projY_radCorr = 0; // radiative corrected data
  TH1D *H_simcPm_projY_PS = 0; // SIMC phase space projection
  TH1D *H_dataPm_projY_Xsec = 0; // data cross sections (online)
  TH1D *H_simcPm_projY_jmlfsi_Xsec = 0; // simc jmlfsi cross sections (online)
  TH1D *H_simcPm_projY_jmlpwia_Xsec = 0; // simc jmlpwia cross sections (online)

  // 1d projections for commissioning 2018 data cross sections
  TH1D *H_comm_dataPm_projY_Xsec_pm80 = 0; 
  TH1D *H_comm_dataPm_projY_Xsec_pm580 = 0;
  TH1D *H_comm_dataPm_projY_Xsec_pm750 = 0; 

  
  cout << "About to loop over xbins of hist2d . . . " << endl;

  //loop over xbins of hist2d 
  for(int i=1; i<=xnbins; i++){

    // get xbin center value and width
    float thrq_center  = hist2d->GetXaxis()->GetBinCenter(i);
    float thrq_width   = hist2d->GetXaxis()->GetBinWidth(i); 

    //cout << "bin: " << i << ", x-val: " << thrq_center << endl;

    // project hist2d along y-axis (different bins in x)
    H_dataPm_projY         = hist2d->ProjectionY(Form("proj_Pm_thrq%.1f_raw", thrq_center), i, i); // raw data counts (no charge normalized or inefficiency corrected)
    
    H_simcPm_projY_ratio = hist2d_Pm_vs_thrq_simc_fsi_ratio->ProjectionY(Form("proj_ratio_simcPm_thrq%.1f", thrq_center), i, i);


    H_dataPm_projY_radUnCorr = hist2d_corr->ProjectionY(Form("proj_Pm_thrq%.1f_radUnCorr", thrq_center), i, i);
    H_dataPm_projY_radCorr = hist2d_Pm_vs_thrq_data_radcorr->ProjectionY(Form("proj_Pm_thrq%.1f_radCorr", thrq_center), i, i);

    
    H_simcPm_projY_PS = hist2d_Pm_vs_thrq_simc_fsi_ps->ProjectionY(Form("proj_simcPm_thrq%.1f_PS", thrq_center), i, i);

    // data cross sections (#ub sr^{-2} MeV^{-1} (based on SIMC phase space units)
    H_dataPm_projY_Xsec = hist2d_Pm_vs_thrq_data_Xsec->ProjectionY(Form("proj_dataPm_thrq%.1f_Xsec", thrq_center), i, i);
    
    // SIMC cross sections
    H_simcPm_projY_jmlfsi_Xsec = hist2d_Pm_vs_thrq_simc_fsi_Xsec->ProjectionY(Form("proj_simcPm_jmlfsi_thrq%.1f_Xsec", thrq_center), i, i);
    H_simcPm_projY_jmlpwia_Xsec = hist2d_Pm_vs_thrq_simc_pwia_Xsec->ProjectionY(Form("proj_simcPm_jmlpwia_thrq%.1f_Xsec", thrq_center), i, i);

    // commissioning data cross sections
    H_comm_dataPm_projY_Xsec_pm80  = hist2d_Pm_vs_thrq_data_Xsec_pm80->ProjectionY(Form("proj_comm_dataPm_thrq%.1f_Xsec_pm80", thrq_center), i, i);
    H_comm_dataPm_projY_Xsec_pm580 = hist2d_Pm_vs_thrq_data_Xsec_pm580->ProjectionY(Form("proj_comm_dataPm_thrq%.1f_Xsec_pm580", thrq_center), i, i);
    H_comm_dataPm_projY_Xsec_pm750 = hist2d_Pm_vs_thrq_data_Xsec_pm750->ProjectionY(Form("proj_comm_dataPm_thrq%.1f_Xsec_pm750", thrq_center), i, i);

    
    // define integrated counts on projected bin
    float counts = H_dataPm_projY->Integral();

    // reduce divisions of histos (for de-cluttering)
    H_dataPm_projY->GetYaxis()->SetNdivisions(5);
    H_dataPm_projY->GetXaxis()->SetNdivisions(10);
    H_dataPm_projY->GetXaxis()->SetLabelSize(0.1);
    
    H_simcPm_projY_ratio->GetYaxis()->SetNdivisions(5);
    H_simcPm_projY_ratio->GetXaxis()->SetNdivisions(10);
    H_simcPm_projY_ratio->GetXaxis()->SetLabelSize(0.1);
    
    H_dataPm_projY_radUnCorr->GetYaxis()->SetNdivisions(5);
    H_dataPm_projY_radUnCorr->GetXaxis()->SetNdivisions(10);
    H_dataPm_projY_radUnCorr->GetXaxis()->SetLabelSize(0.1);
    
    H_dataPm_projY_radCorr->GetYaxis()->SetNdivisions(5);
    H_dataPm_projY_radCorr->GetXaxis()->SetNdivisions(10);
    H_dataPm_projY_radCorr->GetXaxis()->SetLabelSize(0.1);
    
    H_simcPm_projY_PS->GetYaxis()->SetNdivisions(5);
    H_simcPm_projY_PS->GetXaxis()->SetNdivisions(10);
    H_simcPm_projY_PS->GetXaxis()->SetLabelSize(0.1);

    H_dataPm_projY_Xsec->GetYaxis()->SetNdivisions(5);
    H_dataPm_projY_Xsec->GetXaxis()->SetNdivisions(10);
    H_dataPm_projY_Xsec->GetXaxis()->SetLabelSize(0.1);
  
    H_simcPm_projY_jmlfsi_Xsec->GetYaxis()->SetNdivisions(5);
    H_simcPm_projY_jmlfsi_Xsec->GetXaxis()->SetNdivisions(10);
    H_simcPm_projY_jmlfsi_Xsec->GetXaxis()->SetLabelSize(0.1);

    H_simcPm_projY_jmlpwia_Xsec->GetYaxis()->SetNdivisions(5);
    H_simcPm_projY_jmlpwia_Xsec->GetXaxis()->SetNdivisions(10);
    H_simcPm_projY_jmlpwia_Xsec->GetXaxis()->SetLabelSize(0.1);

    H_comm_dataPm_projY_Xsec_pm80->GetYaxis()->SetNdivisions(5);
    H_comm_dataPm_projY_Xsec_pm80->GetXaxis()->SetNdivisions(10);
    H_comm_dataPm_projY_Xsec_pm80->GetXaxis()->SetLabelSize(0.1);

    H_comm_dataPm_projY_Xsec_pm580->GetYaxis()->SetNdivisions(5);
    H_comm_dataPm_projY_Xsec_pm580->GetXaxis()->SetNdivisions(10);
    H_comm_dataPm_projY_Xsec_pm580->GetXaxis()->SetLabelSize(0.1);

    H_comm_dataPm_projY_Xsec_pm750->GetYaxis()->SetNdivisions(5);
    H_comm_dataPm_projY_Xsec_pm750->GetXaxis()->SetNdivisions(10);
    H_comm_dataPm_projY_Xsec_pm750->GetXaxis()->SetLabelSize(0.1);

    
    
    //cout << "thrq_center, counts (v1) = " << thrq_center << ", " << counts << endl;
    H_dataPm_projY->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_dataPm_projY->SetTitleSize(10);
    
    H_simcPm_projY_ratio->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f", thrq_center, thrq_width/2.));
    H_simcPm_projY_ratio->SetTitleSize(10);
    
    H_dataPm_projY_radUnCorr->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_dataPm_projY_radUnCorr->SetTitleSize(10);
    
    H_dataPm_projY_radCorr->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_dataPm_projY_radCorr->SetTitleSize(10);
    
    H_simcPm_projY_PS->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (phase space)", thrq_center, thrq_width/2.));
    H_simcPm_projY_PS->SetTitleSize(10);
    
    H_dataPm_projY_Xsec->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_dataPm_projY_Xsec->SetTitleSize(10);
    
    H_simcPm_projY_jmlfsi_Xsec->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_simcPm_projY_jmlfsi_Xsec->SetTitleSize(10);

    H_simcPm_projY_jmlpwia_Xsec->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    H_simcPm_projY_jmlpwia_Xsec->SetTitleSize(10);
        
    H_comm_dataPm_projY_Xsec_pm80->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f", thrq_center, thrq_width/2.));
    H_comm_dataPm_projY_Xsec_pm80->SetTitleSize(10);

    H_comm_dataPm_projY_Xsec_pm580->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f", thrq_center, thrq_width/2.));
    H_comm_dataPm_projY_Xsec_pm580->SetTitleSize(10);

    H_comm_dataPm_projY_Xsec_pm750->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f", thrq_center, thrq_width/2.));
    H_comm_dataPm_projY_Xsec_pm750->SetTitleSize(10);
    
    
    //cout << Form("bin#: %d,  x-center: %.1f, counts: %.3f", i, thrq_center, counts) << endl;

    
    //-----------------------------------
    // Compute Relative Error on Counts
    //-----------------------------------
    
    // For each H_dataPm_projY hsitogram, get bin content to calculate its error and plot relative error
    y_val.clear();
    y_err.clear();
    x_val.clear();
    x_err.clear();
    
    pm_counts = 0;
    relative_err = 0;
    
    // set statistical lower (inner) limit for guidance during online data-taking
    int inner_stats = 15;  // +/- 15 %
  
    // loop over all bins of Pmiss (H_dataPm_projY histogram) 
    for(int i =1; i<=H_dataPm_projY->GetNbinsX(); i++){

      // get bin content
      pm_counts = H_dataPm_projY->GetBinContent(i);

      // get bin center
      pm_center = H_dataPm_projY->GetBinCenter(i);
      
      // calculate relative error of bin
      relative_err = (sqrt(pm_counts) / pm_counts ) * 100.; // in %

      // push values to vector
      y_val.push_back(0);
      y_err.push_back( relative_err );
      x_val.push_back(pm_center);
      x_err.push_back(0); // no need to set this width
      
      
    }

    // change to TFile for writing to root
    fout->cd();

    // at the end, should have vector of length N for plotting
    int n=H_dataPm_projY->GetNbinsX();
    
    TGraphErrors *gr = new TGraphErrors(n, &x_val[0], &y_val[0], &x_err[0], &y_err[0]);
    gr->SetTitle(Form("proj_Pm_thrq%.1f_relErr", thrq_center));
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerSize(0.);
    gr->SetMarkerStyle(21);

    TLine *lo_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), -inner_stats, *(max_element(x_val.begin(), x_val.end())) , -inner_stats);
    TLine *up_limit = new TLine( *(min_element(x_val.begin(), x_val.end())), inner_stats, *(max_element(x_val.begin(), x_val.end())) , inner_stats);

    lo_limit->SetLineColor(kRed);
    lo_limit->SetLineStyle(1);
    lo_limit->SetLineWidth(1);

    up_limit->SetLineColor(kRed);
    up_limit->SetLineStyle(1);
    up_limit->SetLineWidth(1);
    
  
    
    //cout << "thrq_center, counts (v2) = " << thrq_center << ", " << counts << endl;
    gr->SetTitle(Form("#theta_{rq} = %.1f#pm%.1f (N=%.1f)", thrq_center, thrq_width/2., counts));
    
    // reduce divisions of histos (for de-cluttering)
    gr->GetYaxis()->SetNdivisions(5);
    gr->GetXaxis()->SetNdivisions(10);

    gr->Write();

    //---------------------------------------------------
        
    /*
    c1->cd(i);
    gPad->Modified();
    gPad->Update();

    H_dataPm_projY->SetMarkerStyle(kFullCircle);
    H_dataPm_projY->SetMarkerSize(1);
    H_dataPm_projY->SetMarkerColor(kBlack);
    H_dataPm_projY->SetLineColor(kBlack);
    H_dataPm_projY->Draw("histE0");
    H_dataPm_projY->Write();
   
    
    c3->cd(i);
    gPad->Modified();
    gPad->Update();    

    H_simcPm_projY_ratio->SetMarkerStyle(kFullCircle);
    H_simcPm_projY_ratio->SetMarkerSize(1);
    H_simcPm_projY_ratio->SetMarkerColor(kBlack);
    H_simcPm_projY_ratio->SetLineColor(kBlack);

    H_simcPm_projY_ratio->Draw("PE0");
    H_simcPm_projY_ratio->Write();
    
    
    c4->cd(i);
    gPad->Modified();
    gPad->Update();

    H_dataPm_projY_radUnCorr->SetMarkerStyle(kFullCircle);
    H_dataPm_projY_radCorr->SetMarkerStyle(kFullCircle);
    
    H_dataPm_projY_radUnCorr->SetMarkerSize(1);
    H_dataPm_projY_radCorr->SetMarkerSize(1);
    H_dataPm_projY_radUnCorr->SetMarkerColor(kBlue);
    H_dataPm_projY_radCorr->SetMarkerColor(kRed);
    H_dataPm_projY_radUnCorr->SetLineColor(kBlue);
    H_dataPm_projY_radCorr->SetLineColor(kRed);
    
    H_dataPm_projY_radCorr->Draw("histE0");
    H_dataPm_projY_radUnCorr->Draw("histE0same");
    
 

    // add legend
    if(i==1){
      auto legend2 = new TLegend(0.5,0.2,0.9,0.85);
      legend2->AddEntry("H_dataPm_projY_radUnCorr","no_rad_corr","%s");
      legend2->SetBorderSize(0);
      legend2->SetTextSize(0.08);
      legend2->SetTextColor(kBlue);
      legend2->Draw();

      auto legend3 = new TLegend(0.5,0.6,0.9,0.85);
      legend3->AddEntry("H_dataPm_projY_radCorr","rad_corr","%s");
      legend3->SetBorderSize(0);
      legend3->SetTextSize(0.08);
      legend3->SetTextColor(kRed);
      legend3->Draw("same");
    }

    H_dataPm_projY_radUnCorr->Write();
    H_dataPm_projY_radCorr->Write();

    */ 
    c5->cd(i);
    H_simcPm_projY_PS->SetMarkerStyle(kFullCircle);
    H_simcPm_projY_PS->SetMarkerSize(1);
    H_simcPm_projY_PS->SetMarkerColor(kBlack);
    H_simcPm_projY_PS->SetLineColor(kBlack);

    H_simcPm_projY_PS->Draw("histE0");
    H_simcPm_projY_PS->Write();
    /*
    //----------- online data cross sections -----------
    c6->cd(i);
    gPad->SetLogy();
    gPad->Modified();
    gPad->Update();


    // ---------- PLOT THEORY---------
    
    H_simcPm_projY_jmlfsi_Xsec->SetMarkerStyle(kFullStar);
    H_simcPm_projY_jmlfsi_Xsec->SetMarkerSize(1.3);
    H_simcPm_projY_jmlfsi_Xsec->SetMarkerColor(kGreen+2);
    H_simcPm_projY_jmlfsi_Xsec->Draw("PLC");
  
    H_simcPm_projY_jmlpwia_Xsec->SetMarkerStyle(kFullStar);
    H_simcPm_projY_jmlpwia_Xsec->SetMarkerSize(1.3);
    H_simcPm_projY_jmlpwia_Xsec->SetMarkerColor(kBlue+2);   
    H_simcPm_projY_jmlpwia_Xsec->Draw("PLCsame");
   
    H_simcPm_projY_jmlfsi_Xsec->Write();
    H_simcPm_projY_jmlpwia_Xsec->Write();

    // check if plotting online pm120 MeV setting, then compare it to the pm80 from commissioning
    if(setting=="pm120"){
      H_comm_dataPm_projY_Xsec_pm80->SetMarkerStyle(kFullTriangleUp);
      H_comm_dataPm_projY_Xsec_pm80->SetMarkerSize(1);
      H_comm_dataPm_projY_Xsec_pm80->SetMarkerColor(kRed);
      H_comm_dataPm_projY_Xsec_pm80->SetLineColor(kRed);
      H_comm_dataPm_projY_Xsec_pm80->Draw("PE0same");
      H_comm_dataPm_projY_Xsec_pm80->Write();
    }

    if(setting=="pm580"){
      H_comm_dataPm_projY_Xsec_pm580->SetMarkerStyle(kFullTriangleUp);
      H_comm_dataPm_projY_Xsec_pm580->SetMarkerSize(1);
      H_comm_dataPm_projY_Xsec_pm580->SetMarkerColor(kOrange-3);
      H_comm_dataPm_projY_Xsec_pm580->SetLineColor(kOrange-3);   
      H_comm_dataPm_projY_Xsec_pm580->Draw("PE0same");
      H_comm_dataPm_projY_Xsec_pm580->Write();
    }

    if(setting=="pm800" || setting=="pm900"){
      H_comm_dataPm_projY_Xsec_pm580->SetMarkerStyle(kFullTriangleUp);
      H_comm_dataPm_projY_Xsec_pm580->SetMarkerSize(1);
      H_comm_dataPm_projY_Xsec_pm580->SetMarkerColor(kOrange-3);
      H_comm_dataPm_projY_Xsec_pm580->SetLineColor(kOrange-3);
      H_comm_dataPm_projY_Xsec_pm580->Draw("PE0same");
      H_comm_dataPm_projY_Xsec_pm580->Write();

      H_comm_dataPm_projY_Xsec_pm750->SetMarkerStyle(kFullTriangleUp);
      H_comm_dataPm_projY_Xsec_pm750->SetMarkerSize(1);
      H_comm_dataPm_projY_Xsec_pm750->SetMarkerColor(kMagenta-7);
      H_comm_dataPm_projY_Xsec_pm750->SetLineColor(kMagenta-7);
      H_comm_dataPm_projY_Xsec_pm750->Draw("PE0same");
      H_comm_dataPm_projY_Xsec_pm750->Write();
    }
    
    // add legend
    if(i==1){
      auto legend4 = new TLegend(0.2,0.7, 0.4,0.9);
      legend4->AddEntry("H_dataPm_projY_Xsec","d^{5}#sigma/d#Omega_{e,p}d#omega [#mub sr^{-2} MeV^{-1}]","%s");
      legend4->SetBorderSize(0);
      legend4->SetTextSize(0.09);
      legend4->Draw();

      auto leg_data = new TLegend(0.3,0.6,0.5,0.7);   
      leg_data->AddEntry("H_dataPm_projY_Xsec","data (online)","%s");    
      leg_data->SetBorderSize(0);   
      leg_data->SetTextSize(0.09); 
      leg_data->Draw("same");

      auto leg_pwia = new TLegend(0.3,0.5,0.5,0.6);
      leg_pwia->AddEntry("H_simcPm_projY_jmlpwia_Xsec","JML PWIA","%s");
      leg_pwia->SetBorderSize(0);
      leg_pwia->SetTextSize(0.08);
      leg_pwia->SetTextColor(kBlue+2);
      leg_pwia->Draw();
      
      auto leg_fsi = new TLegend(0.3,0.4,0.5,0.5);
      leg_fsi->AddEntry("H_simcPm_projY_jmlfsi_Xsec","JML FSI","%s");
      leg_fsi->SetBorderSize(0);
      leg_fsi->SetTextSize(0.08);
      leg_fsi->SetTextColor(kGreen+2);
      leg_fsi->Draw("same");

      if(setting=="pm120"){
	auto leg_pm80 = new TLegend(0.3,0.3,0.5,0.4);
	leg_pm80->AddEntry("H_comm_dataPm_projY_Xsec_pm80","data 80 MeV/c (2018)","%s");
	leg_pm80->SetBorderSize(0);
	leg_pm80->SetTextSize(0.08);
	leg_pm80->SetTextColor(kRed);
	leg_pm80->Draw("same");
      }

      if(setting=="pm580"){
	auto leg_pm580 = new TLegend(0.3,0.3,0.5,0.4);
	leg_pm580->AddEntry("H_comm_dataPm_projY_Xsec_pm580","data 580 MeV/c (2018)","%s");
	leg_pm580->SetBorderSize(0);
	leg_pm580->SetTextSize(0.08);
	leg_pm580->SetTextColor(kOrange-3);
	leg_pm580->Draw("same");
      }

      if(setting=="pm800" || setting=="pm900"){

	auto leg_pm580 = new TLegend(0.3,0.3,0.5,0.4);
	leg_pm580->AddEntry("H_comm_dataPm_projY_Xsec_pm580","data 580 MeV/c (2018)","%s");
	leg_pm580->SetBorderSize(0);
	leg_pm580->SetTextSize(0.08);
	leg_pm580->SetTextColor(kOrange-3);
	leg_pm580->Draw("same");
	
	auto leg_pm750 = new TLegend(0.3,0.2,0.5,0.3);
	leg_pm750->AddEntry("H_comm_dataPm_projY_Xsec_pm750","data 750 MeV/c (2018)","%s");
	leg_pm750->SetBorderSize(0);
	leg_pm750->SetTextSize(0.08);
	leg_pm750->SetTextColor(kMagenta-7);
	leg_pm750->Draw("same");

      }
      
    }

    H_dataPm_projY_Xsec->SetMarkerStyle(kFullCircle);
    H_dataPm_projY_Xsec->SetMarkerSize(1);
    H_dataPm_projY_Xsec->SetMarkerColor(kBlack);
    H_dataPm_projY_Xsec->SetLineColor(kBlack);
    H_dataPm_projY_Xsec->Draw("PE0same");
    H_dataPm_projY_Xsec->Write();


    //------------------------------------------------------------------
   
    c2->cd(i); 
   
    // draw to graph
    gr->Draw("AP");
    lo_limit->Draw();
    up_limit->Draw();
     
    // add legend
    if(i==1){
      auto legend = new TLegend(0.5,0.6,0.9,0.8);
      legend->AddEntry("gr",Form("#pm %d %%",inner_stats),"%d");
      legend->SetBorderSize(0);
      legend->SetTextSize(0.15);
      legend->SetTextColor(kRed);
      legend->Draw();
    }
    
    */
  } // end loop over 2D xbins [th_rq]
  
  
  // save canvas
  gStyle->SetOptStat(0);
  c0->SaveAs( fout_2dHist.Data()      );
  c1->SaveAs( fout_projHist.Data()    );
  c2->SaveAs( fout_projHistErr.Data() );
  c3->SaveAs( fout_projsimcRadCorr.Data() );
  c4->SaveAs( fout_projdataRadCorr.Data() );
  c5->SaveAs( fout_projsimcPS.Data() );
  c6->SaveAs( fout_projdataXsec.Data() );

  if(display_plots){
    
    //open plots with evince or any other viewer
    
    gSystem->Exec(Form("evince %s", fout_2dHist.Data() ));
    gSystem->Exec(Form("evince %s", fout_projHist.Data() ));
    gSystem->Exec(Form("evince %s", fout_projHistErr.Data() )); 
    gSystem->Exec(Form("evince %s",  fout_projsimcRadCorr.Data() ));
    gSystem->Exec(Form("evince %s",  fout_projdataRadCorr.Data() ));
    gSystem->Exec(Form("evince %s",  fout_projsimcPS.Data() ));
    gSystem->Exec(Form("evince %s",  fout_projdataXsec.Data() ));


    /*
    gSystem->Exec(Form("open %s", fout_2dHist.Data() ));
    gSystem->Exec(Form("open %s", fout_projHist.Data() ));
    gSystem->Exec(Form("open %s", fout_projHistErr.Data() )); 
    gSystem->Exec(Form("open %s",  fout_projsimcRadCorr.Data() ));
    */
  }
  
}


void project2d_online_develop() {

  int run_min;
  int run_max;
  TString pm_setting;

  cout << "\n Please enter Pmiss setting (e.g. pm120, pm580, pm800, pm900): ";
  cin >> pm_setting;    

  
  if(pm_setting=="pm120"){
    run_min=20871;  
    run_max=20872;  
  }    
  
  if(pm_setting=="pm580"){   
    run_min=20873;
    run_max=20883;   
  }

  if(pm_setting=="pm800"){
    run_min=20886; 
    run_max=20956;
  }
  
  if(pm_setting=="pm900"){
    run_min=20958; 
    run_max=21009;    
  }
  /*
  cout << "" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "" << endl;
  cout << "function: project2d_online(int run_min, int run_max, TString pm_setting)" << endl;
  cout << "" << endl;
  cout << "Brief: This function plots the following: \n"
  "(1) the combined 2d histogram Pmiss vs. th_rq \n"
  "(2) projection of Pmiss (y-axis) in th_rq bins (x-axis)   \n"
  "(3) relative statistical errors (sqrt[N]/[N]) of projected bins\n (for monitoring statistical goals per thrq bin)" << endl;
  cout << "" << endl;
  cout <<  "**NOTE** : if a run range is over a single or multiple settings, \n "
    "name it accordintly e.g. pm_120 (if single), pm_120_580 (multiple settings) etc. \n "
    "This will be used to name the output .pdf \n" << endl;
  cout << "" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "" << endl;
  cout << "Please enter minimum run in range (e.g., 3288): ";
  cin >> run_min;
  cout << "\n Please enter maximum run in range (e.g., 3377): ";
  cin >> run_max;
  cout << "\n Please enter Pmiss setting (e.g. pm_120, pm_120_580): ";
  cin >> pm_setting;
  cout << "" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  */
  // Brief: this function calls two functions described below. 

  
  // funct1: retrieve the combined 2d histos on a given run range. The root file is predefined in the function,
  // so  all the user needs is to input the path where the histogram is on the .root file
  TH2F * myhist2d_corr = combine_2dhistos(run_min, run_max, "randSub_plots/H_Pm_vs_thrq_rand_sub", true); //corrected for ineff. and charge
  TH2F * myhist2d = combine_2dhistos(run_min, run_max, "randSub_plots/H_Pm_vs_thrq_rand_sub", false); // raw counts

  //TCanvas *c1 = new TCanvas("c1");
  //c1->Divide(1,2);
  //c1->cd(1);
  //myhist2d->Draw("colz");
  //c1->cd(2);
  //myhist2d_corr->Draw("colz");

  // func2: projectes the 2d histo (slices of x-bins along y-axis) onto 1d bins, the pm_setting is just for histogram naming purposes
  // and should be consistent with the histogram range chosen
  project2d_deut( myhist2d, myhist2d_corr, pm_setting, false, true ); // the bool flag is to display the plots (otherwise, they will be saved)

}

//#endif
