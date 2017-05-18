#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

using namespace TMath;



void TryEff()
{
  const double ETAMIN  =      1.24;  // Minimum GEN / trigger eta to consider              
  const double ETAMAX  =      2.40;  // Maximum GEN / trigger eta to consider
  const double BIT     =      0.01;  // Shift to trigger pT to get it off bin edge
  const double PTMIN   =       0.0;
  const double PTMAX   =      50.0;
  const int    PTBINS  =        50;


  TString in_file_name = "PtRegression_AWB_v0_17_01_05_no_wgt_inv_pt_orig_deriv_vars_10k.root";
  
  TFile *f = new TFile(in_file_name);
  if(!f) cout << "No data file!"<<endl;

  TTree *t_test = (TTree*)f->Get("dataset/TestTree");
  Float_t bGEN_pt=0, bGEN_eta=0, bEMTF_eta=0, bTRG_pt=0;
  t_test->SetBranchAddress("GEN_pt",&bGEN_pt);
  t_test->SetBranchAddress("GEN_eta",&bGEN_eta);
  t_test->SetBranchAddress("EMTF_eta",&bEMTF_eta);
  //trigger pt = 1/pt
  t_test->SetBranchAddress("BDTG_default",&bTRG_pt);


  cout<<"Nentries: "<<t_test->GetEntries()<<endl;

  //********* HISTO **********//
  TH1F *h_GENpt = new TH1F("h_GENpt","GEN pt",100,0,1000);
  TH1F *h_TRGpt = new TH1F("h_TRGpt","TRG pt",100,0,1000);
  TH2F *h_Dev = new TH2F("h_Dev",";GEN pt (GeV); TRGpt - GENpt (GeV)",100,0,100,100,-10,10);
  TH1D *h_count =  new TH1D( "h_count", "  h_count ", PTBINS, PTMIN, PTMAX );
  TH2D *h_pt = new TH2D( "h_pt", "h_pt", 1, 10, 200, PTBINS, PTMIN, PTMAX);
  TH2D *h_pt_1 = new TH2D( "h_pt_1", "h_pt", PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX);
  TH2D *h_eff = new TH2D("h_eff","h_eff",PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX);
  TH2D *h_eff_eta = new TH2D("h_eff_eta","h_eff_eta",PTBINS, ETAMIN, ETAMAX, 50,0,1);
  TH1D *h_eff1D = new TH1D("h_eff1D","h_eff1D",100,0,2);

  //////////////////////////////////


  for(Int_t iEv=0; iEv<t_test->GetEntries();iEv++)
    {
      t_test->GetEntry(iEv);

      Double_t GEN_pt = (Double_t)bGEN_pt;
      Double_t GEN_eta = (Double_t)bGEN_eta;
      Double_t EMTF_eta = (Double_t)bEMTF_eta;
      Double_t TRG_pt = (Double_t)bTRG_pt;

      if(fabs(EMTF_eta)<ETAMIN) continue;
      if(fabs(EMTF_eta)>ETAMAX) continue;
      if(TRG_pt<0) continue;

      TRG_pt = 1/TRG_pt;
      TRG_pt += BIT;

      h_GENpt->Fill(GEN_pt);
      h_TRGpt->Fill(TRG_pt);
      h_Dev->Fill(GEN_pt,TRG_pt-GEN_pt);

      // Fill counts from ZeroBias events
      if (GEN_eta < -10)
	h_count->Fill(TRG_pt);

      if ( TRG_pt < 0 ) continue;
      if ( fabs(GEN_eta) < ETAMIN ) continue;
      if ( fabs(GEN_eta) > ETAMAX ) continue;
      TRG_pt = min(PTMAX - BIT, TRG_pt);
      GEN_pt = min(PTMAX - BIT, GEN_pt);

      h_pt->Fill( GEN_pt, TRG_pt );

      //Old histo
      h_pt_1->Fill( GEN_pt, TRG_pt );

      float num=0, den = 0, eff=0;

      for (int iBin = 1; iBin <= 1; iBin++) {
	for (int jBin = PTBINS; jBin >= 1; jBin--) { // Loop over trigger pT bins (y-axis)
	  num= h_pt->Integral(iBin, iBin, jBin, PTBINS); // Events passing trigger pT cut
	  den = h_pt->Integral(iBin, iBin,    1, PTBINS); // All events in GEN pT bin
	  eff = num / den;
	  h_eff1D->Fill(eff);
	}
      }
      

      // Compute 2D efficiency histograms and rate at efficiency threshold histograms
      for (int iBin = 1; iBin <= PTBINS; iBin++) { // Loop over GEN pT bins (x-axis)
	for (int jBin = PTBINS; jBin >= 1; jBin--) { // Loop over trigger pT bins (y-axis)
	  float num1 = h_pt_1->Integral(iBin, iBin, jBin, PTBINS); // Events passing trigger pT cut
	  float den1 = h_pt_1->Integral(iBin, iBin,    1, PTBINS); // All events in GEN pT bin
	  float eff1 = num1 / den1;
	  h_eff_eta->Fill(GEN_eta,eff1);
	  h_eff->SetBinContent(iBin, jBin, eff1);
	  h_eff->SetBinError  (iBin, jBin, eff1 * sqrt( (1/num1) + (1/den1) ) );
	}//Loop over y-axis
      }//Loop over x-axis

      //for(Int_t iBin =0;iBin<PTBINS;iBin++)
	

    }//end loop on tree

  //////////// DRAW HISTO ///////////////////
  TCanvas *c = new TCanvas("c","",1000,500);
  c->Divide(2,1);
  c->cd(1);
  h_GENpt->Draw("hist");
  c->cd(2);
  h_TRGpt->Draw("hist");

  TCanvas *c1 = new TCanvas("c1","",500,500);
  h_Dev->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","",1000,500);
  c2->Divide(2,1);
  c2->cd(1);
  h_count->Draw("hist");
  c2->cd(2);
  h_pt->Draw("colz");

  TCanvas *c3  = new TCanvas("c3","",500,500);
  h_eff->Draw("colz");

  TCanvas *c4  = new TCanvas("c4","",500,500);
  h_eff1D->Draw("hist");

  TCanvas *c5  = new TCanvas("c5","",500,500);
  h_eff_eta->Draw("hist");



}//THE End


