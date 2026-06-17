void get_boiling(Int_t runNum, Double_t &set_current, Double_t &avg_current_bcm4a, Double_t &scl_Yield_bcm4a,Double_t &scl_err_bcm4a, 
		 Double_t &notrk_Yield_bcm4a, Double_t &notrk_err_bcm4a, Double_t &trk_Yield_bcm4a, Double_t &trk_err_bcm4a,  
		 Double_t &cpuLT_bcm4a, Double_t &tLT_bcm4a, Double_t &Qtot_bcm4a, Double_t &pS1XscalerRate_bcm4a, Double_t &pTRIG1scalerRate_bcm4a, 
		 Double_t &pTRIG2scalerRate_bcm4a, Double_t &pTRIG3scalerRate_bcm4a, Double_t &pEDTMscalerRate_bcm4a, Double_t &eTrkEff_bcm4a,
		 Double_t &avg_current_bcm4b, Double_t &scl_Yield_bcm4b,Double_t &scl_err_bcm4b, 
		 Double_t &notrk_Yield_bcm4b, Double_t &notrk_err_bcm4b, Double_t &trk_Yield_bcm4b, Double_t &trk_err_bcm4b,  
		 Double_t &cpuLT_bcm4b, Double_t &tLT_bcm4b, Double_t &Qtot_bcm4b, Double_t &pS1XscalerRate_bcm4b, Double_t &pTRIG1scalerRate_bcm4b, 
		 Double_t &pTRIG2scalerRate_bcm4b, Double_t &pTRIG3scalerRate_bcm4b, Double_t &pEDTMscalerRate_bcm4b, Double_t &eTrkEff_bcm4b
		 )
{

  TH1F::SetDefaultSumw2(kTRUE);

  //gROOT->SetBatch(kTRUE);

  //Define Constants
  //Int_t runNum = 2093;
  Double_t current_thrs_bcm4a = 3.;
  Double_t current_thrs_bcm4b = 3.;

  Double_t pS;   //Pre-Scale FActor


  //SHMS LD2 Boiling Study
  if(runNum==21049){
    set_current = 10.;
    pS = 1.;
    //current_thrs_bcm4a = 5.; 
    //current_thrs_bcm4b = 5.; 
    //eTrkEff = 0.9931;  
  }
  if(runNum==21050){
    set_current = 15.;
    pS = 1.;
    //eTrkEff = 0.9932;  
  }
  if(runNum==21051){
    set_current = 25.;
    pS = 1.;
    //eTrkEff = 0.9933;  
  }
  if(runNum==21052){
    set_current = 35.;
    pS = 2.;
    //eTrkEff = 0.9929;  
  }
  if(runNum==21053){
    set_current = 40.;
    pS = 2.;
    //eTrkEff = 0.9932;  
  }
  if(runNum==21054){
    set_current = 40.;
    pS = 2.;
    //eTrkEff = 0.9932;  
  }
  if(runNum==21055){
    set_current = 45.;
    pS = 2.;
    //eTrkEff = 0.9933;  
  }
  if(runNum==21056){
    set_current = 45.;
    pS = 3.;
    //eTrkEff = 0.9933;  
  }
  if(runNum==21057){
    set_current = 65.;
    pS = 3.;
    //eTrkEff = 0.9933;  
  }
  if(runNum==21058){
    set_current = 70.;
    pS = 3.;
    //eTrkEff = 0.9933;  
  }
  
  //SHMS Carbon Boiling Study
  if(runNum==21059){
    set_current = 15.;
    pS = 1.;
    //eTrkEff = 0.9936;  
  }
  if(runNum==21060){
    set_current = 30.;
    pS = 1.;
    //eTrkEff = 0.9930;  
  }
  if(runNum==21061){
    set_current = 40.;
    pS = 1.;
    //eTrkEff = 0.9936;  
  }
  if(runNum==21062){
    set_current = 55.;
    pS = 1.;
    //eTrkEff = 0.9936;  
  }
  if(runNum==21063){
    set_current = 70.;
    pS = 1.;
    //eTrkEff = 0.9936;  
  }
  //SHMS Al. dummy
  if(runNum==2096){
    set_current = 40;
    pS = 1.;
    //eTrkEff = 0.9936;  
  }

  //SHMS LH2 Boiling Study
  if(runNum==2075){
    set_current = 80.;
    pS = 1.;
    current_thrs_bcm4a = 5.; 
    current_thrs_bcm4b = 5.; 
    //eTrkEff = 0.9932;
  }
  if(runNum==2076){
    set_current = 70.;
    pS = 1.;
    //eTrkEff = 0.9932;
  }
  if(runNum==2078){
    set_current =  10.;
    pS = 1.;
    //eTrkEff = 0.9935;
  }
  if(runNum==2080){
    set_current = 10.;
    pS = 1.;
    //eTrkEff = 0.9930;
  }
  if(runNum==2081){
    set_current = 20.;
    pS = 1.;
    //eTrkEff = 0.9934;
  }
  if(runNum==2082){
    set_current = 35.;
    pS = 1.;
    //eTrkEff = 0.9952;
  }
  if(runNum==2083){
    set_current = 35.;
    pS = 1.;
    //eTrkEff = 0.9933;
  }
  if(runNum==2084){
    set_current = 45.;
    pS = 1.;
    //eTrkEff = 0.9924;
  }
  if(runNum==2085){
    set_current = 55.;
    pS = 1.;
    //eTrkEff = 0.9932;
  }


  //Created and Read Rootfile
  TString filename=Form("/media/gvill/Gema's T7/ROOTfiles/pass_2/deut_replay_prod_%d_-1.root", runNum);
  
  TFile *file = new TFile(filename);
  
  TTree *tdata = (TTree*) file->Get("T");    //data TTree
  TTree *tscal = (TTree*) file->Get("TSP");  //Scaler TTree


  //Get Scaler Leafs
  Double_t Scal_evNum;    //event number associated with scaler reads
  tscal->SetBranchAddress("evNumber", &Scal_evNum);
  Double_t  Scal_BCM4A_charge;
  tscal->SetBranchAddress("P.BCM4A.scalerCharge",&Scal_BCM4A_charge);
  Double_t  Scal_BCM4A_current;
  tscal->SetBranchAddress("P.BCM4A.scalerCurrent",&Scal_BCM4A_current);  
  Double_t  Scal_BCM4B_charge;
  tscal->SetBranchAddress("P.BCM4B.scalerCharge",&Scal_BCM4B_charge);
  Double_t  Scal_BCM4B_current;
  tscal->SetBranchAddress("P.BCM4B.scalerCurrent",&Scal_BCM4B_current);  
  Double_t  Scal_time;
  tscal->SetBranchAddress("P.1MHz.scalerTime",&Scal_time);
  Double_t  pS1X_scaler;
  tscal->SetBranchAddress("P.S1X.scaler",&pS1X_scaler);         //SHMS S1X 
  Double_t  pTRIG1_scaler;
  tscal->SetBranchAddress("P.pTRIG1.scaler",&pTRIG1_scaler);    //SHMS 3/4
  Double_t  pTRIG2_scaler;  
  tscal->SetBranchAddress("P.pTRIG2.scaler",&pTRIG2_scaler);    //SHMS EL-REAL
  Double_t  pTRIG3_scaler; 
  tscal->SetBranchAddress("P.pTRIG3.scaler",&pTRIG3_scaler);    //SHMS EL-CLEAN
  Double_t  pEDTM_scaler; 
  tscal->SetBranchAddress("P.EDTM.scaler",&pEDTM_scaler);    //SHMS EDTM 
  

  //Defive Quantities To Store Previous Reads and cumulative quantities
  Double_t prev_time = 0.;
  Double_t prev_charge_bcm4a = 0.;
  Double_t prev_charge_bcm4b = 0.;
  Double_t prev_ps1x_scaler = 0;
  Double_t prev_ptrig1_scaler = 0;
  Double_t prev_ptrig2_scaler = 0;
  Double_t prev_ptrig3_scaler = 0;
  Double_t prev_pedtm_scaler = 0;

  Double_t total_time = 0.;
  Double_t total_charge_bcm4a = 0.;
  Double_t total_charge_bcm4b = 0.;
  Double_t total_ps1x_scaler = 0;
  Double_t total_ptrig1_scaler = 0;
  Double_t total_ptrig2_scaler = 0;
  Double_t total_ptrig3_scaler = 0;
  Double_t total_pedtm_scaler = 0;

  //BCM4A Cut variables
  Double_t total_time_bcm4a_cut = 0.;
  Double_t total_charge_bcm4a_cut = 0.;
  Double_t total_ps1x_scaler_bcm4a_cut = 0;
  Double_t total_ptrig1_scaler_bcm4a_cut = 0;
  Double_t total_ptrig2_scaler_bcm4a_cut = 0;
  Double_t total_ptrig3_scaler_bcm4a_cut = 0;
  Double_t total_pedtm_scaler_bcm4a_cut = 0;

  //BCM4B Cut variables
  Double_t total_time_bcm4b_cut = 0.;
  Double_t total_charge_bcm4b_cut = 0.;
  Double_t total_ps1x_scaler_bcm4b_cut = 0;
  Double_t total_ptrig1_scaler_bcm4b_cut = 0;
  Double_t total_ptrig2_scaler_bcm4b_cut = 0;
  Double_t total_ptrig3_scaler_bcm4b_cut = 0;
  Double_t total_pedtm_scaler_bcm4b_cut = 0;

  //Loop Over Scaler Reads 
  Long64_t scal_entries = tscal->GetEntries();

  Int_t evt_flag_bcm4a[scal_entries];             //Store Flag [0 or 1], to know if scaler read passed current cut (1) or not (0)
  Int_t evt_flag_bcm4b[scal_entries];             //Store Flag [0 or 1], to know if scaler read passed current cut (1) or not (0)

  Int_t scal_evt_num[scal_entries];    //Store Event Associated with Scaler Read
 
  for (int i = 0; i < scal_entries; i++) {
    
    //**NOTE: Each scaler read is associated with as specific event number
    //        as (scaler read 1-> event 1000,  scaler read 2 -> event 2300, ...)
    //        This means events up to 1000 correspond to scaler read 1, ...
   
    tscal->GetEntry(i);
       
    
    evt_flag_bcm4a[i] = 0;
    evt_flag_bcm4b[i] = 0;
    
    scal_evt_num[i] = Scal_evNum; //store events associated with scaler reads (to be used in data tree)

    //Store Cumulative Quantities
    total_charge_bcm4a = Scal_BCM4A_charge;
    total_charge_bcm4b = Scal_BCM4B_charge;
    total_time = Scal_time;
    total_ps1x_scaler = pS1X_scaler;
    total_ptrig1_scaler = pTRIG1_scaler;
    total_ptrig2_scaler = pTRIG2_scaler;
    total_ptrig3_scaler = pTRIG3_scaler;
    total_pedtm_scaler = pEDTM_scaler;

    //Check If BCM4A Beam Current in Between Reads is Over Threshold
    if(abs(Scal_BCM4A_current-set_current)<current_thrs_bcm4a)
      {

	//Turn Event Flag ON, if beam current is within threshold
	evt_flag_bcm4a[i] = 1;

	//Store Quantities that Passed the Current Threshold
	total_time_bcm4a_cut = total_time_bcm4a_cut + (Scal_time - prev_time);
	total_charge_bcm4a_cut = total_charge_bcm4a_cut + (Scal_BCM4A_charge - prev_charge_bcm4a);  
	total_ps1x_scaler_bcm4a_cut = total_ps1x_scaler_bcm4a_cut + (pS1X_scaler-prev_ps1x_scaler);
	total_ptrig1_scaler_bcm4a_cut = total_ptrig1_scaler_bcm4a_cut + (pTRIG1_scaler-prev_ptrig1_scaler);
	total_ptrig2_scaler_bcm4a_cut = total_ptrig2_scaler_bcm4a_cut + (pTRIG2_scaler-prev_ptrig2_scaler);
	total_ptrig3_scaler_bcm4a_cut = total_ptrig3_scaler_bcm4a_cut + (pTRIG3_scaler-prev_ptrig3_scaler);
	total_pedtm_scaler_bcm4a_cut = total_pedtm_scaler_bcm4a_cut + (pEDTM_scaler - prev_pedtm_scaler);
      }
    
    //Check If BCM4B Beam Current in Between Reads is Over Threshold
    if(abs(Scal_BCM4B_current-set_current)<current_thrs_bcm4b)
      {
	
	//Turn Event Flag ON, if beam current is within threshold
	evt_flag_bcm4b[i] = 1;

	//Store Quantities that Passed the Current Threshold
	total_time_bcm4b_cut = total_time_bcm4b_cut + (Scal_time - prev_time);
	total_charge_bcm4b_cut = total_charge_bcm4b_cut + (Scal_BCM4B_charge - prev_charge_bcm4b);  
	total_ps1x_scaler_bcm4b_cut = total_ps1x_scaler_bcm4b_cut + (pS1X_scaler-prev_ps1x_scaler);
	total_ptrig1_scaler_bcm4b_cut = total_ptrig1_scaler_bcm4b_cut + (pTRIG1_scaler-prev_ptrig1_scaler);
	total_ptrig2_scaler_bcm4b_cut = total_ptrig2_scaler_bcm4b_cut + (pTRIG2_scaler-prev_ptrig2_scaler);
	total_ptrig3_scaler_bcm4b_cut = total_ptrig3_scaler_bcm4b_cut + (pTRIG3_scaler-prev_ptrig3_scaler);
	total_pedtm_scaler_bcm4b_cut = total_pedtm_scaler_bcm4b_cut + (pEDTM_scaler - prev_pedtm_scaler);
      }

    //Previous Scaler Reads (Necessary to Take Average between S-1 and S scaler reads, to get values in between)
    prev_time = Scal_time;
    prev_charge_bcm4a = Scal_BCM4A_charge;
    prev_charge_bcm4b = Scal_BCM4B_charge;
    prev_ps1x_scaler = pS1X_scaler;
    prev_ptrig1_scaler = pTRIG1_scaler;
    prev_ptrig2_scaler = pTRIG2_scaler;
    prev_ptrig3_scaler = pTRIG3_scaler;
    prev_pedtm_scaler = pEDTM_scaler;

  } //end scaler read loop

  //Calculate Rates if Current Cut Passed
  //BCM4A
  pS1XscalerRate_bcm4a = total_ps1x_scaler_bcm4a_cut / total_time_bcm4a_cut;
  pTRIG1scalerRate_bcm4a = total_ptrig1_scaler_bcm4a_cut / total_time_bcm4a_cut;
  pTRIG2scalerRate_bcm4a = total_ptrig2_scaler_bcm4a_cut / total_time_bcm4a_cut;
  pTRIG3scalerRate_bcm4a = total_ptrig3_scaler_bcm4a_cut / total_time_bcm4a_cut;
  pEDTMscalerRate_bcm4a =  total_pedtm_scaler_bcm4a_cut / total_time_bcm4a_cut;
  //BCM4B
  pS1XscalerRate_bcm4b = total_ps1x_scaler_bcm4b_cut / total_time_bcm4b_cut;
  pTRIG1scalerRate_bcm4b = total_ptrig1_scaler_bcm4b_cut / total_time_bcm4b_cut;
  pTRIG2scalerRate_bcm4b = total_ptrig2_scaler_bcm4b_cut / total_time_bcm4b_cut;
  pTRIG3scalerRate_bcm4b = total_ptrig3_scaler_bcm4b_cut / total_time_bcm4b_cut;
  pEDTMscalerRate_bcm4b =  total_pedtm_scaler_bcm4b_cut / total_time_bcm4b_cut;

  

  //*******************
  //      

  //Get Data Leafs
  Double_t gevtyp;
  tdata->SetBranchAddress("g.evtyp",&gevtyp);
  Double_t gevnum;
  tdata->SetBranchAddress("g.evnum",&gevnum);
  
  
  Double_t pTRIG1_tdcTimeRaw;
  tdata->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw",&pTRIG1_tdcTimeRaw);
  Double_t pTRIG2_tdcTimeRaw;
  tdata->SetBranchAddress("T.coin.pTRIG2_ROC2_tdcTimeRaw",&pTRIG2_tdcTimeRaw);
  Double_t pTRIG3_tdcTimeRaw;
  tdata->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTimeRaw",&pTRIG3_tdcTimeRaw);
  Double_t pEDTM_tdcTimeRaw;
  tdata->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&pEDTM_tdcTimeRaw);
  Double_t pHGCER_npeSum; 
  tdata->SetBranchAddress("P.hgcer.npeSum",&pHGCER_npeSum);
  Double_t pNGCER_npeSum; 
  tdata->SetBranchAddress("P.ngcer.npeSum",&pNGCER_npeSum);
  Double_t pCAL_etotnorm; 
  tdata->SetBranchAddress("P.cal.etotnorm",&pCAL_etotnorm);    
  Double_t pCAL_etottracknorm; 
  tdata->SetBranchAddress("P.cal.etottracknorm",&pCAL_etottracknorm);    
  Double_t pDelta;
  tdata->SetBranchAddress("P.gtr.dp",&pDelta);    
  Double_t pBetanotrk;
  tdata->SetBranchAddress("P.hod.betanotrack",&pBetanotrk);    
  Double_t pGoodScinHit;
  tdata->SetBranchAddress("P.hod.goodscinhit",&pGoodScinHit);    
  Double_t pdc_ntrack;
  tdata->SetBranchAddress("P.dc.ntrack",&pdc_ntrack);    

  //Create Histograms 
  TH1F *H_pHGCER_npeSum = new TH1F("H_pHGCER_npeSum", "SHMS HGCer NpeSum", 100,0.1,10);
  TH1F *H_pCAL_etotnorm = new TH1F("H_pCAL_etotnorm", "SHMS Cal EtotNorm", 100,0.1,1.5);
  TH1F *H_pCAL_etottracknorm = new TH1F("H_pCAL_etottracknorm", "SHMS Cal EtotTrackNorm", 100,0.1,1.5);
  TH1F *H_pDelta = new TH1F("H_pDelta", "SHMS Delta", 100,-15,25);
  
  TH1F *H_pHGCER_npeSum_cut = new TH1F("H_pHGCER_npeSum_cut", "SHMS HGCer NpeSum", 100,0.1,10);
  TH1F *H_pCAL_etotnorm_cut = new TH1F("H_pCAL_etotnorm_cut", "SHMS Cal EtotNorm", 100,0.1,1.5);
  TH1F *H_pCAL_etottracknorm_cut = new TH1F("H_pCAL_etottracknorm_cut", "SHMS Cal EtotTrackNorm", 100,0.1,1.5);
  TH1F *H_pDelta_cut = new TH1F("H_pDelta_cut", "SHMS Delta_cut", 100,-15,25);

  //BCM4A
  TH1F *H_pCAL_etotnorm_ntrkY_4a = new TH1F("H_pCAL_etotnorm_ntrkY_4a", "SHMS Cal. E_{tot} Norm, Yield_{ntrk}", 100, 0.1, 2.0);
  TH1F *H_pDelta_trkY_4a = new TH1F("H_pDelta_trkY_4a", "SHMS #delta, Yield_{trk}", 100, -15, 25);
  //BCM4B
  TH1F *H_pCAL_etotnorm_ntrkY_4b = new TH1F("H_pCAL_etotnorm_ntrkY_4b", "SHMS Cal. E_{tot} Norm, Yield_{ntrk}", 100, 0.1, 2.0);
  TH1F *H_pDelta_trkY_4b = new TH1F("H_pDelta_trkY_4b", "SHMS #delta, Yield_{trk}", 100, -15, 25);
  
  //Set-Up Tdc Counters
  Double_t total_ptrig1_accp = 0;
  Double_t total_ptrig2_accp = 0;
  Double_t total_ptrig3_accp = 0;
  Double_t total_pedtm_accp = 0;
  //Count ONLY for BCM4A Cut
  Double_t total_ptrig1_accp_cut_4a = 0;
  Double_t total_ptrig2_accp_cut_4a = 0;
  Double_t total_ptrig3_accp_cut_4a = 0;
  Double_t total_pedtm_accp_cut_4a = 0;
  //Count ONLY for BCM4B Cut
  Double_t total_ptrig1_accp_cut_4b = 0;
  Double_t total_ptrig2_accp_cut_4b = 0;
  Double_t total_ptrig3_accp_cut_4b = 0;
  Double_t total_pedtm_accp_cut_4b = 0;

  //Set Up Yield Counter
  //BCM4A
  Double_t Y_ntrk_4a;
  Double_t Y_trk_4a;
  //BCM4B
  Double_t Y_ntrk_4b;
  Double_t Y_trk_4b;


  //Tracking Efficiency Counter
  //BCM4A
  Double_t e_did_4a=0.;
  Double_t e_should_4a=0.;
  Double_t eTrkEff_err_4a;
  //BCM4B
  Double_t e_did_4b=0.;
  Double_t e_should_4b=0.;
  Double_t eTrkEff_err_4b;
  
  //Loop Over Data Events 
  Int_t scal_read = 0; //scaler read counter

  Long64_t data_entries = tdata->GetEntries();
  //Double_t gevnum[data_entries];
  
  //for (int j = 0; j < data_entries; j++){
  //    gevnum[j] = j+1;
  //}
  
  
    for (int i = 0; i < data_entries; i++) {

      
      tdata->GetEntry(i);

     

      //Define Cuts
      Bool_t c_noedtm = pEDTM_tdcTimeRaw==0;
      Bool_t c_edtm = pEDTM_tdcTimeRaw>0;
      Bool_t c_ptrig2 = pTRIG2_tdcTimeRaw>0;
      Bool_t c_cerNpesum =  pHGCER_npeSum>0.7;
      Bool_t c_pdelta = pDelta>-10.&&pDelta<22.;
      Bool_t c_etotnorm = pCAL_etotnorm>0.6;
      Bool_t c_etottrknorm = pCAL_etottracknorm>0.6;
      //e- tracking efficiency 
      Bool_t good_elec_should = pGoodScinHit==1.&&pBetanotrk>0.85&&pBetanotrk<1.5&&pCAL_etotnorm>0.6&&pHGCER_npeSum>0.7;
      Bool_t good_elec_did = good_elec_should&&pdc_ntrack>0.;


      //Count Accepted EDTM events
      if(c_edtm){
	total_pedtm_accp++;
      }
      //Count Accepted pTRIG2 events
      if(c_ptrig2){
	total_ptrig2_accp++;
      }
      
      //--------Check If BCM4A Current is within limits---------
      if(evt_flag_bcm4a[scal_read]==1){
       

	//Count Accepted EDTM events
	if(c_edtm){
	  total_pedtm_accp_cut_4a++;
	}

	//Count Accepted pTRIG2 events (without EDTM)
	if(c_ptrig2&&c_noedtm){
	  total_ptrig2_accp_cut_4a++;
	}
      
	if(c_noedtm) {

	  /*
	  //Fill DATA Histograms
	  H_pHGCER_npeSum->Fill(pHGCER_npeSum);
	  H_pCAL_etotnorm->Fill(pCAL_etotnorm);
	  H_pCAL_etottracknorm->Fill(pCAL_etottracknorm);
	  H_pDelta->Fill(pDelta);
	  //Fill With Cuts 
	  if(c_cerNpesum){ H_pHGCER_npeSum_cut->Fill(pHGCER_npeSum);}
	  if(c_etotnorm){H_pCAL_etotnorm_cut->Fill(pCAL_etotnorm);}
	  if(c_etottrknorm){H_pCAL_etottracknorm_cut->Fill(pCAL_etottracknorm);}
	  if(c_pdelta){H_pDelta_cut->Fill(pDelta);}
	  */

	  //Make Cuts to Count Yield
	  //---Get NON-TRACKING Yield
	  if(c_cerNpesum&&c_etotnorm)
	    {
	      H_pCAL_etotnorm_ntrkY_4a->Fill(pCAL_etotnorm);
	      Y_ntrk_4a++;
	  
	    }
	  
	  if(c_etottrknorm&&c_pdelta)
	    {
	      H_pDelta_trkY_4a->Fill(pDelta);
	      Y_trk_4a++;
	    }

	  //Calculate Electron Tracking Efficiency
	  if(good_elec_did){
	    e_did_4a++;
	  }
	  if(good_elec_should){
	    e_should_4a++;
	  }

	} //end cut on NO EDTM


      } 
      //-------------End BCM4A Current Cut----------------
     

      //--------Check If BCM4B Current is within limits---------
      if(evt_flag_bcm4b[scal_read]==1){
       

	//Count Accepted EDTM events
	if(c_edtm){
	  total_pedtm_accp_cut_4b++;
	}

	//Count Accepted pTRIG2 events (without EDTM)
	if(c_ptrig2&&c_noedtm){
	  total_ptrig2_accp_cut_4b++;
	}
      
	if(c_noedtm) {

	  
	  //Fill DATA Histograms
	  H_pHGCER_npeSum->Fill(pHGCER_npeSum);
	  H_pCAL_etotnorm->Fill(pCAL_etotnorm);
	  H_pCAL_etottracknorm->Fill(pCAL_etottracknorm);
	  H_pDelta->Fill(pDelta);
	  //Fill With Cuts 
	  if(c_cerNpesum){ H_pHGCER_npeSum_cut->Fill(pHGCER_npeSum);}
	  if(c_etotnorm){H_pCAL_etotnorm_cut->Fill(pCAL_etotnorm);}
	  if(c_etottrknorm){H_pCAL_etottracknorm_cut->Fill(pCAL_etottracknorm);}
	  if(c_pdelta){H_pDelta_cut->Fill(pDelta);}
	  

	  //Make Cuts to Count Yield
	  //---Get NON-TRACKING Yield
	  if(c_cerNpesum&&c_etotnorm)
	    {
	      H_pCAL_etotnorm_ntrkY_4b->Fill(pCAL_etotnorm);
	      Y_ntrk_4b++;
	  
	    }
	  
	  if(c_etottrknorm&&c_pdelta)
	    {
	      H_pDelta_trkY_4b->Fill(pDelta);
	      Y_trk_4b++;
	    }

	  //Calculate Electron Tracking Efficiency
	  if(good_elec_did){
	    e_did_4b++;
	  }
	  if(good_elec_should){
	    e_should_4b++;
	  }

	} //end cut on NO EDTM


      } 
      //-------------End BCM4B Current Cut----------------



      //Increment Scaler Read if event == scaler_evt_pperlimit for that scaler read
      if(gevnum==scal_evt_num[scal_read]){
	scal_read++;
      }
  
    } //end data loop
    
    /*
    cout << "Total Accepted EDTM (CUT)= " << total_pedtm_accp_cut_4a << endl;
    cout << "Total Accepted pTRIG2 (CUT) = " << total_ptrig2_accp_cut_4a << endl;
    cout << "Total Accepted EDTM = " << total_pedtm_accp << endl;
    cout << "Total Accepted pTRIG2 = " << total_ptrig2_accp << endl;
    */


    //Draw Histograms
    TCanvas *c1 = new TCanvas("c1", "", 1500,1500);
    c1->Divide(2,2);
    c1->cd(1);
    H_pHGCER_npeSum->Draw();
    H_pHGCER_npeSum_cut->SetLineColor(kRed);
    H_pHGCER_npeSum_cut->Draw("sames");
    c1->cd(2);
    H_pCAL_etotnorm->Draw();
    H_pCAL_etotnorm_cut->SetLineColor(kRed);
    H_pCAL_etotnorm_cut->Draw("sames");
    c1->cd(3);
    H_pCAL_etottracknorm->Draw();
    H_pCAL_etottracknorm_cut->SetLineColor(kRed);
    H_pCAL_etottracknorm_cut->Draw("sames");
    c1->cd(4);
    H_pDelta->Draw();
    H_pDelta_cut->SetLineColor(kRed);
    H_pDelta_cut->Draw("sames");



    
    //Calculate  Live Time
    //BCM4A
    cpuLT_bcm4a =  total_ptrig2_accp_cut_4a / (total_ptrig2_scaler_bcm4a_cut-total_pedtm_scaler_bcm4a_cut);
    tLT_bcm4a =  total_pedtm_accp_cut_4a / total_pedtm_scaler_bcm4a_cut;
    //BCM4B
    cpuLT_bcm4b =  total_ptrig2_accp_cut_4b / (total_ptrig2_scaler_bcm4b_cut-total_pedtm_scaler_bcm4b_cut);
    tLT_bcm4b =  total_pedtm_accp_cut_4b / total_pedtm_scaler_bcm4b_cut;
    
    //Calculate Tracking Efficiency
    //BCM4A
    eTrkEff_bcm4a = e_did_4a / e_should_4a;
    eTrkEff_err_4a = sqrt(e_should_4a-e_did_4a) / e_should_4a;
    //BCM4B
    eTrkEff_bcm4b = e_did_4b / e_should_4b;
    eTrkEff_err_4b = sqrt(e_should_4b-e_did_4b) / e_should_4b;
    
    //Avg Current
    avg_current_bcm4a = total_charge_bcm4a_cut / total_time_bcm4a_cut; //uA
    avg_current_bcm4b = total_charge_bcm4b_cut / total_time_bcm4b_cut; //uA

    total_charge_bcm4a_cut = total_charge_bcm4a_cut/1000.;  //COnvert from uC to mC
    Qtot_bcm4a = total_charge_bcm4a_cut;

    total_charge_bcm4b_cut = total_charge_bcm4b_cut/1000.;  //COnvert from uC to mC
    Qtot_bcm4b = total_charge_bcm4b_cut;
    
    //Calculate Normalized Yields
    //Yield involving ONLY scaler quantities  ---> Y / Q
    Double_t notrk_Yield;//, notrk_Yield_bcm4a, notrk_err_bcm4a;    //Yield counting data variables, NO TRACKING ----> Y / (Q*tLT_bcm4a)
    Double_t trk_Yield;//, trk_Yield_bcm4a, trk_err_bcm4a;      //Yield counting data variables,  TRACKING ------> Y / (Q*tLT_bcm4a*trkEff)
    
    //-----SCALER YIELD-----
    //BCM4A
    scl_Yield_bcm4a = (total_ptrig2_scaler_bcm4a_cut-total_pedtm_scaler_bcm4a_cut) / (total_charge_bcm4a_cut);
    scl_err_bcm4a = sqrt(total_ptrig2_scaler_bcm4a_cut) / total_charge_bcm4a_cut;
    //BCM4B
    scl_Yield_bcm4b = (total_ptrig2_scaler_bcm4b_cut-total_pedtm_scaler_bcm4b_cut) / (total_charge_bcm4b_cut);
    scl_err_bcm4b = sqrt(total_ptrig2_scaler_bcm4b_cut) / total_charge_bcm4b_cut;

    
    //-----NON-TRACKING Yield---
    notrk_Yield = Y_ntrk_4a / (total_charge_bcm4a_cut*cpuLT_bcm4a);

    //BCM4A
    H_pCAL_etotnorm_ntrkY_4a->Scale(1./(total_charge_bcm4a_cut*cpuLT_bcm4a));
    notrk_Yield_bcm4a = H_pCAL_etotnorm_ntrkY_4a->IntegralAndError(H_pCAL_etotnorm_ntrkY_4a->FindBin(0.6), H_pCAL_etotnorm_ntrkY_4a->FindBin(5.), notrk_err_bcm4a);
    //BCM4B
    H_pCAL_etotnorm_ntrkY_4b->Scale(1./(total_charge_bcm4b_cut*cpuLT_bcm4b));
    notrk_Yield_bcm4b = H_pCAL_etotnorm_ntrkY_4b->IntegralAndError(H_pCAL_etotnorm_ntrkY_4b->FindBin(0.6), H_pCAL_etotnorm_ntrkY_4b->FindBin(5.), notrk_err_bcm4b);


    
    //-----TRACKING Yield---
    trk_Yield = Y_trk_4a / (total_charge_bcm4a_cut*cpuLT_bcm4a*eTrkEff_bcm4a);
    //BCM4A
    H_pDelta_trkY_4a->Scale(1./(total_charge_bcm4a_cut*cpuLT_bcm4a*eTrkEff_bcm4a));
    trk_Yield_bcm4a = H_pDelta_trkY_4a->IntegralAndError(H_pDelta_trkY_4a->FindBin(-8.), H_pDelta_trkY_4a->FindBin(8.), trk_err_bcm4a);    
    //BCM4B
    H_pDelta_trkY_4b->Scale(1./(total_charge_bcm4b_cut*cpuLT_bcm4b*eTrkEff_bcm4a));
    trk_Yield_bcm4b = H_pDelta_trkY_4b->IntegralAndError(H_pDelta_trkY_4b->FindBin(-8.), H_pDelta_trkY_4b->FindBin(8.), trk_err_bcm4b);  
}


void run_boiling_shms(string target)
{

  //User Input: target-> "Al", "C12", "LH2",  "LD2"

  Int_t runNUM;
  string line;

  //Read FileName
  TString filename = Form("shms_%s.dat", target.c_str());
  ifstream ifs;
  ifs.open(filename);
  //int C12_runlist[5] = {21059,21060,21061,21062,21063};

  //Create Output File to write results
  ofstream ofs;
  ofs.open(Form("shms_%s_yield.data", target.c_str()));
  ofs << "#!Run[i,0]/    set_current[f,1]/   avg_current_bcm4a[f,2]/     sclY4a[f,3]/     sclY_err4a[f,4]/    ntrkY4a[f,5]/    ntrkY_err4a[f,6]/    trkY4a[f,7]/    trkY_err4a[f,8]/    cpuLT_bcm4a[f,9]/     tLT_bcm4a[f,10]/    Qtot_bcm4a[f,11]/    ps1x_rate_bcm4a[f,12]/    ptrig1_rate_bcm4a[f,13]/     ptrig2_rate_bcm4a[f,14]/     ptrig3_rate_bcm4a[f,15]/    pedtm_rate_bcm4a[f,16]/   p_eTrkEff_bcm4a[f,17]/      avg_current_bcm4b[f,18]/     sclY4b[f,19]/     sclY_err4b[f,20]/    ntrkY4b[f,21]/    ntrkY_err4b[f,22]/    trkY4b[f,23]/    trkY_err4b[f,24]/    cpuLT_bcm4b[f,25]/     tLT_bcm4b[f,26]/    Qtot_bcm4b[f,27]/     ps1x_rate_bcm4b[f,28]/      ptrig1_rate_bcm4b[f,29]/     ptrig2_rate_bcm4b[f,30]/     ptrig3_rate_bcm4b[f,31]/    pedtm_rate_bcm4b[f,32]/   p_eTrkEff_bcm4b[f,33]/" << endl;
    
  
  //Variables to be passed to function
  Double_t set_current;
  Double_t sclY4a, sclY4a_err;
  Double_t ntrkY4a, ntrkY4a_err;
  Double_t trkY4a, trkY4a_err;
  Double_t avg_current_bcm4a;
  Double_t cpuLT_bcm4a;
  Double_t tLT_bcm4a;
  Double_t Qtot_bcm4a;
  Double_t ptrig1_rate_bcm4a;
  Double_t ptrig2_rate_bcm4a;
  Double_t ptrig3_rate_bcm4a;
  Double_t pedtm_rate_bcm4a;
  Double_t p_eTrkEff_bcm4a;
  Double_t ps1x_rate_bcm4a;

  //BCM4B
  Double_t sclY4b, sclY4b_err;
  Double_t ntrkY4b, ntrkY4b_err;
  Double_t trkY4b, trkY4b_err;
  Double_t avg_current_bcm4b;
  Double_t cpuLT_bcm4b;
  Double_t tLT_bcm4b;
  Double_t Qtot_bcm4b;
  Double_t ptrig1_rate_bcm4b;
  Double_t ptrig2_rate_bcm4b;
  Double_t ptrig3_rate_bcm4b;
  Double_t pedtm_rate_bcm4b;
  Double_t p_eTrkEff_bcm4b;
  Double_t ps1x_rate_bcm4b;
  
  //Read run txt file
  while (getline(ifs, line)){
  
  //for (int i = 0; i < 5; i++) { 
    //convert run from string to int
    runNUM = stoi(line);  
    //runNUM = C12_runlist[i];
    
    get_boiling(runNUM, set_current,  avg_current_bcm4a, sclY4a, sclY4a_err, ntrkY4a, ntrkY4a_err, trkY4a, trkY4a_err,
		cpuLT_bcm4a, tLT_bcm4a, Qtot_bcm4a, ps1x_rate_bcm4a, ptrig1_rate_bcm4a, ptrig2_rate_bcm4a, ptrig3_rate_bcm4a, pedtm_rate_bcm4a, p_eTrkEff_bcm4a,
		avg_current_bcm4b, sclY4b, sclY4b_err, ntrkY4b, ntrkY4b_err, trkY4b, trkY4b_err,
		cpuLT_bcm4b, tLT_bcm4b, Qtot_bcm4b, ps1x_rate_bcm4b, ptrig1_rate_bcm4b, ptrig2_rate_bcm4b, ptrig3_rate_bcm4b, pedtm_rate_bcm4b, p_eTrkEff_bcm4b);
    
    //Write to file
    ofs << Form("%i     %f     %f      %f     %f     %f     %f     %f     %f     %f     %f     %f     %f    %f    %f     %f     %f    %f    %f    %f     %f     %f     %f     %f     %f     %f     %f     %f    %f    %f     %f     %f    %f    %f", runNUM, set_current, avg_current_bcm4a, sclY4a, sclY4a_err, ntrkY4a, ntrkY4a_err, trkY4a, trkY4a_err, cpuLT_bcm4a, tLT_bcm4a, Qtot_bcm4a, ps1x_rate_bcm4a, ptrig1_rate_bcm4a, ptrig2_rate_bcm4a, ptrig3_rate_bcm4a, pedtm_rate_bcm4a, p_eTrkEff_bcm4a, avg_current_bcm4b, sclY4b, sclY4b_err, ntrkY4b, ntrkY4b_err, trkY4b, trkY4b_err, cpuLT_bcm4b, tLT_bcm4b, Qtot_bcm4b, ps1x_rate_bcm4b, ptrig1_rate_bcm4b, ptrig2_rate_bcm4b, ptrig3_rate_bcm4b, pedtm_rate_bcm4b, p_eTrkEff_bcm4b) << endl;



  }
  
  ifs.close();
  ofs.close();
  
}
