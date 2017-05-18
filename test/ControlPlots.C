#include "TFile.h"
#include "TTree.h"

using namespace TMath;

const int MAX_EVT =   10000;
const int REPORT_EVT = 1000;
const double BIT = 0.000001; // Tiny value or offset
const double PI = 3.14159265359;

void ControlPlots()
{
  //HISTOS
  TH1F *h_eta_all = new TH1F("h_eta_all","",50,1,2);
  TH1F *h_eta_4hits = new TH1F("h_eta_4hits","",50,1,2);

  TH1F *h_nHits = new TH1F("h_nHits","h_nHits",5,0,5);
  TH1F *h_dphi12_1 = new TH1F("h_dphi12_1","",80,-1,0.2);
  TH1F *h_dphi12_2 = new TH1F("h_dphi12_2","",80,-1,0.2);
  TH1F *h_dphi12_3 = new TH1F("h_dphi12_3","",80,-1,0.2);

  TH1F *h_dphi12_1a = new TH1F("h_dphi12_1a","",80,-1,0.2);
  TH1F *h_dphi12_2a = new TH1F("h_dphi12_2a","",80,-1,0.2);
  TH1F *h_dphi12_3a = new TH1F("h_dphi12_3a","",80,-1,0.2);

  TH1F *h_dR = new TH1F("h_dR","",100,0,200);
  TH1F *h_dR_min = new TH1F("h_dR_min","",100,0,200);

  TH1F *h_dEta = new TH1F("h_dEta","",50,-0.2,0.2);
  TH1F *h_dPhi = new TH1F("h_dPhi","",50,-0.5,0.5);
  TH1F *h_dPt = new TH1F("h_dEta","",50,-1,1);

  TFile *ifile = new TFile("EMTF_MC_NTuple_SingleMu_noRPC_300k.root","READ");
  if (!ifile) {
    cout << "ERROR: could not open data file!" <<endl;
    exit(1);
  }
  if(ifile) cout<<"File found"<<endl;

   // Have to use TChain for both SetBranchAddress and GetEntry to work
   TChain *regChain = new TChain("ntuple/tree");
   regChain->Add("EMTF_MC_NTuple_SingleMu_noRPC_300k.root");

   TBranch *muon_br = regChain->GetBranch("muon");
   TBranch *hit_br   = regChain->GetBranch("hit");
   TBranch *track_br = regChain->GetBranch("track");


   Int_t nevtBin[40];
   Int_t nevtBinOK[40];
   Float_t eff[40],etaMin[40],etaMax[40];
   for(Int_t i=0;i<40;i++)
     {
       nevtBin[i]=0;
       nevtBinOK[i]=0;
       eff[i]=0;
       etaMin[0] = 1.2;
       etaMax[0] = 1.22;
       
       if(i>0){
	 etaMin[i]=etaMin[i-1]+0.02;
	 etaMax[i] =etaMax[i-1]+0.02;
       }
     }
   

   for (UInt_t iEvt = 0; iEvt < regChain->GetEntries(); iEvt++) {if (iEvt > MAX_EVT) break;
     if ( (iEvt % REPORT_EVT) == 0 )
       cout << "Looking at event " << iEvt <<endl;
     regChain->GetEntry(iEvt);

     UInt_t nMuons  = (muon_br->GetLeaf("nMuons"))->GetValue();
     UInt_t nTracks = (track_br->GetLeaf("nTracks"))->GetValue();
     cout << "There are " << nMuons << " GEN muons and " << nTracks << " EMTF tracks\n" << endl;

    for (UInt_t iMu = 0; iMu < nMuons; iMu++) 
      {
	Double_t mu_eta = (muon_br->GetLeaf("eta"))->GetValue(iMu);
	Double_t mu_phi = (muon_br->GetLeaf("phi"))->GetValue(iMu);
	Double_t mu_pt  = (muon_br->GetLeaf("pt"))->GetValue(iMu);
	Int_t mu_charge    = (muon_br->GetLeaf("charge"))->GetValue(iMu);
	if (mu_pt<10 || mu_pt>200) continue;
	if ( fabs( mu_eta ) < 1.10 ) continue;
	if ( fabs( mu_eta ) > 2.50 ) continue;
	
	Double_t dR[nTracks];
	
	for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) 
	  {
	    Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk);
	    Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk);
	    Double_t trk_theta = (track_br->GetLeaf("theta"))->GetValue(iTrk);
	    Double_t trk_phi   = (track_br->GetLeaf("phi"))->GetValue(iTrk);

	    Double_t dEta = trk_eta-mu_eta;
	    Double_t dPhi = trk_phi-mu_phi;
	    Double_t dPt = trk_pt-mu_pt;

	    dR[iTrk]= Sqrt(dEta*dEta+dPhi*dPhi);
	    h_dR->Fill(dR[iTrk]);
	    cout<<" dR: "<<dR[iTrk]<<endl;
	  }

	//I look for the track which has the lowest dR and I save its index
	Double_t dRmin = 1E10;
	Int_t iTrk_ok=0;
	for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) 
	  {
	    if(dR[iTrk]<dRmin)
	      {
		dRmin = dR[iTrk];
		iTrk_ok = iTrk;
	      }
	  }
	cout<<" dRmin: "<<dRmin<<endl;
	h_dR_min->Fill(dRmin);
       
	Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk_ok);
	Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk_ok);
	Double_t trk_theta = (track_br->GetLeaf("theta"))->GetValue(iTrk_ok);
	Double_t trk_phi   = (track_br->GetLeaf("phi"))->GetValue(iTrk_ok);
	Int_t trk_nHit    = (track_br->GetLeaf("nHits"))->GetValue(iTrk_ok);
	Int_t trk_charge    = (track_br->GetLeaf("charge"))->GetValue(iTrk_ok);
	Double_t phi1 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 0);
	Double_t phi2 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 1);
	Double_t phi3 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 2);
	Double_t phi4 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk_ok + 3);

	Double_t dEta = trk_eta-mu_eta;
	Double_t dPhi = trk_phi-mu_phi;
	Double_t dPt = trk_pt-mu_pt;

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
	
	h_eta_all->Fill(mu_eta);
	if(trk_nHit==4)
	  h_eta_4hits->Fill(mu_eta);
	for(Int_t i=0;i<40;i++)
	  {
	    if(mu_eta>etaMin[i]&&mu_eta<=etaMax[i])
	      {
		nevtBin[i]++;
		if(trk_nHit==4)
		  nevtBinOK[i]++;   
	      }
	  }
      }
   
   }
   
   Float_t eta[40];
   for(Int_t i=0;i<40;i++)
     {
       eff[i]=nevtBinOK[i]*1./nevtBin[i];
       eta[i] = (etaMin[i]+etaMax[i])/2.;
     }
   
   TCanvas *c = new TCanvas("c","",500,500);
   /*
   TGraph *g = new TGraph(40,eta,eff);
   g->SetFillColor(kBlue-1);
   g->Draw("AB");
   */
   gStyle->SetOptStat(0);

   TH1F *h_eff= new TH1F("h_eff","",50,1,2);
   h_eff->Divide(h_eta_4hits, h_eta_all);
   h_eff->SetLineColor(kBlue);
   h_eff->GetXaxis()->SetTitle("GEN muon |#eta|");
   h_eff->GetYaxis()->SetTitle("Efficiency");
   h_eff->GetXaxis()->SetRangeUser(1.14,1.96);
   h_eff->Draw();

   TCanvas *c0 = new TCanvas("c0","",500,500);
   h_nHits->Draw("hist");
   
   TCanvas *c1 = new TCanvas("c1","",500,500);
   h_dphi12_3->SetLineColor(kViolet+2);
   h_dphi12_3->GetYaxis()->SetRangeUser(0,100);
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
   h_dphi12_3a->GetYaxis()->SetRangeUser(0,50);
   h_dphi12_3a->SetLineWidth(2);
   h_dphi12_3a->Draw();
   h_dphi12_2a->SetLineWidth(2);
   h_dphi12_2a->SetLineColor(kBlue);
   h_dphi12_2a->Draw("SAMES");
   h_dphi12_1a->SetLineWidth(2);
   h_dphi12_1a->SetLineColor(kGreen);
   h_dphi12_1a->Draw("SAMES");

   TCanvas *c3 = new TCanvas("c3","",500,500);
   h_dR->Draw();
   h_dR_min->SetFillColor(kRed);
   h_dR_min->SetFillStyle(3005);
   h_dR_min->Draw("sames");

   TCanvas *c4 = new TCanvas("c4","",1500,500);
   c4->Divide(3,1);
   c4->cd(1);
   h_dPhi->Draw();
   c4->cd(2);
   h_dEta->Draw();
   c4->cd(3);
   h_dPt->Draw();
}
