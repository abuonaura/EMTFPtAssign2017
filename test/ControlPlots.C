#include "TFile.h"
#include "TTree.h"

using namespace TMath;

const int MAX_EVT =   1000000;
const int REPORT_EVT = 1000;
const double BIT = 0.000001; // Tiny value or offset
const double PI = 3.14159265359;

void ControlPlots()
{
  //HISTOS
  TH1F *h_eta_all_CSC = new TH1F("h_eta_all_CSC","",50,1,2);
  TH1F *h_eta_4hits_CSC = new TH1F("h_eta_4hits_CSC","",50,1,2);
  TH1F *h_eta_all_RPC = new TH1F("h_eta_all_RPC","",50,1,2);
  TH1F *h_eta_4hits_RPC = new TH1F("h_eta_4hits_RPC","",50,1,2);

  TH1F *h_nHits = new TH1F("h_nHits","h_nHits",5,0,5);
  TH1F *h_dphi12_1 = new TH1F("h_dphi12_1","",60,-1,0.2);
  TH1F *h_dphi12_2 = new TH1F("h_dphi12_2","",60,-1,0.2);
  TH1F *h_dphi12_3 = new TH1F("h_dphi12_3","",60,-1,0.2);

  TH1F *h_dphi12_1a = new TH1F("h_dphi12_1a","",60,-1,0.2);
  TH1F *h_dphi12_2a = new TH1F("h_dphi12_2a","",60,-1,0.2);
  TH1F *h_dphi12_3a = new TH1F("h_dphi12_3a","",60,-1,0.2);

  TH1F *h_dR_CSC = new TH1F("h_dR_CSC","",100,0,7);
  TH1F *h_dR_CSC_1 = new TH1F("h_dR_CSC_1","",50,0,1);
  TH1F *h_dR_min_CSC = new TH1F("h_dR_min_CSC","",50,0,1);
  TH1F *h_dR_RPC = new TH1F("h_dR_RPC","",100,0,7);
  TH1F *h_dR_RPC_1 = new TH1F("h_dR_RPC_1","",50,0,1);
  TH1F *h_dR_min_RPC = new TH1F("h_dR_min_RPC","",50,0,1);

  TH1F *h_eta_gen_CSC = new TH1F("h_eta_gen_CSC","CSC only;#eta (rad);",100,0,7);
  TH1F *h_eta_emtf_CSC = new TH1F("h_eta_emtf_CSC","CSC only;#eta (rad);",100,0,7);
  TH1F *h_eta_gen_RPC = new TH1F("h_eta_gen_RPC","CSC + RPC;#eta (rad);",100,0,7);
  TH1F *h_eta_emtf_RPC = new TH1F("h_eta_emtf_RPC","CSC + RPC;#eta (rad);",100,0,7);

  TH1F *h_phi_gen_CSC = new TH1F("h_phi_gen_CSC","CSC only;#phi (rad)",100,0,7);
  TH1F *h_phi_emtf_CSC = new TH1F("h_phi_emtf_CSC","CSC only;#phi (rad);",100,0,7);
  TH1F *h_phi_gen_RPC = new TH1F("h_phi_gen_RPC","CSC + RPC;#phi (rad);",100,0,7);
  TH1F *h_phi_emtf_RPC = new TH1F("h_phi_emtf_RPC","CSC + RPC;#phi (rad);",100,0,7);

  TH1F *h_pt_gen_CSC = new TH1F("h_pt_gen_CSC","CSC only; p_{t} (GeV/c);",100,0,200);
  TH1F *h_pt_emtf_CSC = new TH1F("h_pt_emtf_CSC","CSC only; p_{t} (GeV/c);",100,0,200);
  TH1F *h_pt_gen_RPC = new TH1F("h_pt_gen_RPC","CSC + RPC; p_{t} (GeV/c);",100,0,200);
  TH1F *h_pt_emtf_RPC = new TH1F("h_pt_emtf_RPC","CSC + RPC; p_{t} (GeV/c);",100,0,200);

  TH1F *h_dEta = new TH1F("h_dEta","",50,-0.2,0.2);
  TH1F *h_dPhi = new TH1F("h_dPhi","",50,-0.5,0.5);
  TH1F *h_dPt = new TH1F("h_dPt","",50,-1,1);

  //Read input file
  TFile *ifile;
  TChain *regChain;
  for(Int_t j=0;j<2;j++)
    {
      cout<< j<<endl;
      
      if(j==0)
	ifile = new TFile("EMTF_MC_NTuple_SingleMu_noRPC_300k.root","READ");
      if(j==1)
	ifile = new TFile("EMTF_MC_NTuple_SingleMu_RPC_300k.root","READ");
      if (!ifile) {
      cout << "ERROR: could not open data file!" <<endl;
      	exit(1);
      }
      if(ifile) cout<<"File found"<<endl;
      
      // Have to use TChain for both SetBranchAddress and GetEntry to work
      if(j==0)
	{
	  regChain = new TChain("ntuple/tree");
	  regChain->Add("EMTF_MC_NTuple_SingleMu_noRPC_300k.root");
	}
      if(j==1)
      	{
	  regChain = new TChain("ntuple/tree");
	  regChain->Add("EMTF_MC_NTuple_SingleMu_RPC_300k.root");
	  }
      
      TBranch *muon_br = regChain->GetBranch("muon");
      TBranch *hit_br   = regChain->GetBranch("hit");
      TBranch *track_br = regChain->GetBranch("track");
      
      
      for (UInt_t iEvt = 0; iEvt < regChain->GetEntries(); iEvt++) {if (iEvt > MAX_EVT) break;
	if ( (iEvt % REPORT_EVT) == 0 )
	  cout << "Looking at event " << iEvt <<endl;
	regChain->GetEntry(iEvt);
	
	UInt_t nMuons  = (muon_br->GetLeaf("nMuons"))->GetValue();
	UInt_t nTracks = (track_br->GetLeaf("nTracks"))->GetValue();
	//	cout << "There are " << nMuons << " GEN muons and " << nTracks << " EMTF tracks\n" << endl;

	for (UInt_t iMu = 0; iMu < nMuons; iMu++) 
	  {
	    Double_t mu_eta = (muon_br->GetLeaf("eta"))->GetValue(iMu);
	    Double_t mu_phi = (muon_br->GetLeaf("phi"))->GetValue(iMu);
	    mu_phi = PI*mu_phi/180.;
	    Double_t mu_pt  = (muon_br->GetLeaf("pt"))->GetValue(iMu);
	    Int_t mu_charge    = (muon_br->GetLeaf("charge"))->GetValue(iMu);

	    if(j==0)
	      {
		h_eta_gen_CSC->Fill(mu_eta);
		h_phi_gen_CSC->Fill(mu_phi);
		h_pt_gen_CSC->Fill(mu_pt);
	      }

	    if(j==1)
              {
                h_eta_gen_RPC->Fill(mu_eta);
                h_phi_gen_RPC->Fill(mu_phi);
                h_pt_gen_RPC->Fill(mu_pt);
              }


	    if (mu_pt<10 || mu_pt>200) continue;
	    if ( fabs( mu_eta ) < 1.10 ) continue;
	    if ( fabs( mu_eta ) > 2.50 ) continue;
	
	    Double_t dR[nTracks];

	    //I look for the track which has the lowest dR and I save its index
	    Double_t dRmin = 1E10;
	    Int_t iTrk_ok=0;	
	    for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) 
	      {
		Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk);
		Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk);
		Double_t trk_theta = (track_br->GetLeaf("theta"))->GetValue(iTrk);
		Double_t trk_phi   = (track_br->GetLeaf("phi"))->GetValue(iTrk);
		trk_phi = PI*trk_phi/180.;

		Double_t dEta = trk_eta-mu_eta;
		Double_t dPhi = trk_phi-mu_phi;
		Double_t dPt = trk_pt-mu_pt;

		dR[iTrk]= Sqrt(dEta*dEta+dPhi*dPhi);
		if(j==0)
		  {
		    h_dR_CSC->Fill(dR[iTrk]);
		    h_dR_CSC_1->Fill(dR[iTrk]);
		  }
		if(j==1)
		  {
		    h_dR_RPC->Fill(dR[iTrk]);
		    h_dR_RPC_1->Fill(dR[iTrk]);
		  }
		//cout<<" dR: "<<dR[iTrk]<<endl;
		if(dR[iTrk]<dRmin)
		  {
		    dRmin = dR[iTrk];
		    iTrk_ok = iTrk;
		  }
	      }
	    if(dRmin<PI)
	      {

		//	cout<<" dRmin: "<<dRmin<<endl
		Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk_ok);
		Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk_ok);
		Double_t trk_theta = (track_br->GetLeaf("theta"))->GetValue(iTrk_ok);
		Double_t trk_phi   = (track_br->GetLeaf("phi"))->GetValue(iTrk_ok);
		trk_phi = PI*trk_phi/180.;
		Int_t trk_nHit    = (track_br->GetLeaf("nHits"))->GetValue(iTrk_ok);
		Int_t trk_charge    = (track_br->GetLeaf("charge"))->GetValue(iTrk_ok);
		Double_t phi1 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 0);
		Double_t phi2 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 1);
		Double_t phi3 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 2);
		Double_t phi4 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 3);
		
		Double_t dEta = trk_eta-mu_eta;
		Double_t dPhi = trk_phi-mu_phi;
		Double_t dPt = trk_pt-mu_pt;

		Double_t dPhi12 = acos( cos( (phi2 - phi1)*(PI/180.) ) );
		dPhi12 *= ( sin( (phi2 - phi1)*(PI/180.) ) / max( BIT, abs( sin( (phi2 - phi1)*(PI/180.) ) ) ) );
		dPhi12 *= (180./PI);
		Double_t dPhi13 = acos( cos( (phi3 - phi1)*(PI/180.) ) );
		dPhi13 *= ( sin( (phi3 - phi1)*(PI/180.) ) / max( BIT, abs( sin( (phi3 - phi1)*(PI/180.) ) ) ) );
		dPhi13 *= (180./PI);
		Double_t dPhi14 = acos( cos( (phi4 - phi1)*(PI/180.) ) );
		dPhi14 *= ( sin( (phi4 - phi1)*(PI/180.) ) / max( BIT, abs( sin( (phi4 - phi1)*(PI/180.) ) ) ) );
		dPhi14 *= (180./PI);
		Double_t dPhi23 = acos( cos( (phi3 - phi2)*(PI/180.) ) );
		dPhi23 *= ( sin( (phi3 - phi2)*(PI/180.) ) / max( BIT, abs( sin( (phi3 - phi2)*(PI/180.) ) ) ) );
		dPhi23 *= (180./PI);
		Double_t dPhi24 = acos( cos( (phi4 - phi2)*(PI/180.) ) );
		dPhi24 *= ( sin( (phi4 - phi2)*(PI/180.) ) / max( BIT, abs( sin( (phi4 - phi2)*(PI/180.) ) ) ) );
		dPhi24 *= (180./PI);
		Double_t dPhi34 = acos( cos( (phi4 - phi3)*(PI/180.) ) );
		dPhi34 *= ( sin( (phi4 - phi3)*(PI/180.) ) / max( BIT, abs( sin( (phi4 - phi3)*(PI/180.) ) ) ) );
		dPhi34 *= (180./PI);
		
		if(j==0)
		  {
		    h_dR_min_CSC->Fill(dRmin);
		    
		    h_eta_emtf_CSC->Fill(trk_eta);
		    h_phi_emtf_CSC->Fill(trk_phi);
		    h_pt_emtf_CSC->Fill(trk_pt);

		    h_dEta->Fill(dEta/mu_eta);
		    h_dPt->Fill(dPt/mu_pt);
		    h_dPhi->Fill(dPhi/mu_phi);
		    
		    
		    //cout<<trk_nHit<<endl;
		    h_nHits->Fill(trk_nHit);
		    
		    /*
		      Double_t dPhi12 = phi2 -phi1;
		      Double_t dPhi32 = phi3 -phi2;
		      Double_t dPhi43 = phi4 -phi3;
		    */

		    
		    if(mu_eta>1.2&&mu_eta<1.55)
		      {
			if(mu_pt>15&&mu_pt<30)
			  h_dphi12_1->Fill(dPhi12*trk_charge);
			if(mu_pt>30&&mu_pt<60)
			  h_dphi12_2->Fill(dPhi12*trk_charge);
			if(mu_pt>120&&mu_pt<250)
			  h_dphi12_3->Fill(dPhi12*trk_charge);
		      }
		    
		    if(mu_eta>2.1&&mu_eta<2.4)
		      {
			if(mu_pt>15&&mu_pt<30)
			  h_dphi12_1a->Fill(dPhi12*trk_charge);
			if(mu_pt>30&&mu_pt<60)
			  h_dphi12_2a->Fill(dPhi12*trk_charge);
			if(mu_pt>120&&mu_pt<250)
			  h_dphi12_3a->Fill(dPhi12*trk_charge);
		      }
		    
		   
		    h_eta_all_CSC->Fill(mu_eta);
		    if(trk_nHit==4)
		      h_eta_4hits_CSC->Fill(mu_eta);
		  }
		if(j==1)
		  {
		    h_eta_emtf_RPC->Fill(trk_eta);
                    h_phi_emtf_RPC->Fill(trk_phi);
                    h_pt_emtf_RPC->Fill(trk_pt);

		    h_eta_all_RPC->Fill(mu_eta);
		    if(trk_nHit==4)
		      h_eta_4hits_RPC->Fill(mu_eta);
		    h_dR_min_RPC->Fill(dRmin);
		  }
	      }

	  }
    
      }
    }   
  
  gStyle->SetOptStat(0);   
  TCanvas *c = new TCanvas("c","",500,500);
  TH1F *h_eff_CSC= new TH1F("h_eff_CSC","",50,1,2);
  h_eff_CSC->Divide(h_eta_4hits_CSC, h_eta_all_CSC);
  h_eff_CSC->SetLineColor(kBlue);
  h_eff_CSC->SetLineWidth(2);
  h_eff_CSC->GetXaxis()->SetTitle("GEN muon |#eta|");
  h_eff_CSC->GetYaxis()->SetTitle("Efficiency");
  h_eff_CSC->GetXaxis()->SetRangeUser(1.14,1.96);
  h_eff_CSC->GetYaxis()->SetRangeUser(0,1.1);
  h_eff_CSC->Draw();
  TH1F *h_eff_RPC= new TH1F("h_eff_RPC","",50,1,2);
  h_eff_RPC->Divide(h_eta_4hits_RPC, h_eta_all_RPC);
  h_eff_RPC->SetLineColor(kRed);
  h_eff_RPC->SetLineWidth(2);
  h_eff_RPC->GetXaxis()->SetTitle("GEN muon |#eta|");
  h_eff_RPC->GetYaxis()->SetTitle("Efficiency");
  h_eff_RPC->GetXaxis()->SetRangeUser(1.14,1.96);
  h_eff_RPC->Draw("SAMES");

  TCanvas *c0 = new TCanvas("c0","",500,500);
  h_nHits->Draw("hist");
   
  TCanvas *c1 = new TCanvas("c1","",500,500);
  h_dphi12_3->SetLineColor(kViolet+2);
  h_dphi12_3->GetYaxis()->SetRangeUser(0,1000);
  h_dphi12_3->SetLineWidth(2);
  h_dphi12_3->Draw();
  h_dphi12_2->SetLineWidth(2);
  h_dphi12_2->SetLineColor(kBlue);
  h_dphi12_2->Draw("SAMES");
  h_dphi12_1->SetLineWidth(2);
  h_dphi12_1->SetLineColor(kGreen);
  h_dphi12_1->Draw("SAMES");

  TCanvas *c2 = new TCanvas("c2","",500,500);
  h_dphi12_3a->SetLineColor(kViolet+2);
  h_dphi12_3a->GetYaxis()->SetRangeUser(0,1000);
  h_dphi12_3a->SetLineWidth(2);
  h_dphi12_3a->Draw();
  h_dphi12_2a->SetLineWidth(2);
  h_dphi12_2a->SetLineColor(kBlue);
  h_dphi12_2a->Draw("SAMES");
  h_dphi12_1a->SetLineWidth(2);
  h_dphi12_1a->SetLineColor(kGreen);
  h_dphi12_1a->Draw("SAMES");

  TCanvas *c3 = new TCanvas("c3","",1000,500);
  c3->Divide(2,1);
  TPad *p1 = (TPad*)c3->cd(1);
  p1->SetLogy();
  h_dR_CSC->Draw();
  TPad *p1a = (TPad*)c3->cd(2);
  p1a->SetLogy();
  h_dR_CSC_1->Draw();
  h_dR_min_CSC->SetFillColor(kRed);
  h_dR_min_CSC->SetFillStyle(3005);
  h_dR_min_CSC->Draw("sames");

  TCanvas *c3a = new TCanvas("c3a","",1000,500);
  c3a->Divide(2,1);
  TPad *p2 = (TPad*)c3a->cd(1);
  p2->SetLogy();
  h_dR_RPC->Draw();
  TPad *p2a = (TPad*)c3a->cd(2);
  p2a->SetLogy();
  h_dR_RPC_1->Draw();
  h_dR_min_RPC->SetFillColor(kRed);
  h_dR_min_RPC->SetFillStyle(3005);
  h_dR_min_RPC->Draw("sames");


  TCanvas *c4 = new TCanvas("c4","",1500,500);
  c4->Divide(3,1);
  c4->cd(1);
  h_dPhi->Draw();
  c4->cd(2);
  h_dEta->Draw();
  c4->cd(3);
  h_dPt->Draw();

  TCanvas *c5 = new TCanvas("c5","",1000,500);
  c5->Divide(2,1);
  c5->cd(1);
  h_eta_gen_CSC->Draw();
  h_eta_emtf_CSC->SetLineColor(kRed);
  h_eta_emtf_CSC->Draw("sames");
  c5->cd(2);
  h_eta_gen_RPC->Draw();
  h_eta_emtf_RPC->SetLineColor(kRed);
  h_eta_emtf_RPC->Draw("sames");

  TCanvas *c6 = new TCanvas("c6","",1000,500);
  c6->Divide(2,1);
  c6->cd(1);
  h_phi_gen_CSC->Draw();
  h_phi_emtf_CSC->SetLineColor(kRed);
  h_phi_emtf_CSC->Draw("sames");
  c6->cd(2);
  h_phi_gen_RPC->Draw();
  h_phi_emtf_RPC->SetLineColor(kRed);
  h_phi_emtf_RPC->Draw("sames");

  TCanvas *c7 = new TCanvas("c7","",1000,500);
  c7->Divide(2,1);
  c7->cd(1);
  h_pt_gen_CSC->Draw();
  h_pt_emtf_CSC->SetLineColor(kRed);
  h_pt_emtf_CSC->Draw("sames");
  c7->cd(2);
  h_pt_gen_RPC->Draw();
  h_pt_emtf_RPC->SetLineColor(kRed);
  h_pt_emtf_RPC->Draw("sames");


}
