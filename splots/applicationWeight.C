#define applicationWeight_cxx
#include "applicationWeight.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   
#include <iostream>  

// MC distributions w/wo application of weights computed comparing data and MC

void applicationWeight::Loop(const char* ratioFileName, int q2bin)
{
  // Define q^2 ranges
  float q2inf = 0;
  float q2sup = 1000;
  if (q2bin==0) { q2inf=1.05; q2sup=2.45; }
  if (q2bin==1) { q2inf=2.90; q2sup=3.20; }
  if (q2bin==2) { q2inf=3.55; q2sup=3.80; }
  if (q2bin==3) { q2inf=4.00; q2sup=4.80; }

  if (fChain == 0) return;

  // Histos: from file and from loop
  TFile *ratioFile = new TFile(ratioFileName);
  // 
  TH1F *h1mc_xgb_fromLoop_now = new TH1F("h1mc_xgb_fromLoop_now","h1mc_xgb_fromLoop_now",17,-2.,15.);
  TH1F *h1mc_xgb_fromLoop_ww  = new TH1F("h1mc_xgb_fromLoop_ww", "h1mc_xgb_fromLoop_ww", 17,-2.,15.);
  TH1F *h1mc_xgb_fromFile    = (TH1F*)ratioFile->Get("MC_bdt");
  TH1F *h1data_xgb_fromFile  = (TH1F*)ratioFile->Get("Data_bdt");
  TH1F *h1ratio_xgb_fromFile = (TH1F*)ratioFile->Get("ratio_bdt");
  h1mc_xgb_fromLoop_now -> Sumw2();
  h1mc_xgb_fromLoop_ww  -> Sumw2();

  // Load weights
  for (int i = 0; i<h1ratio_xgb_fromFile->GetNbinsX(); i++) {

    float myweight=h1ratio_xgb_fromFile->GetBinContent(i+1);
    anweight_.push_back(myweight);              

    float mylowedge=h1ratio_xgb_fromFile->GetBinLowEdge(i+1); 
    lowedge_.push_back(mylowedge);
  }


  // Counters
  int totEvents   = 0;
  int totEventsWW = 0;
  int totEventsBdtGt5   = 0;
  int totEventsBdtGt5WW = 0;

  // Loop over entries
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "nentries = " << nentries << endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    // Same selection as applied in data
    if (Bmass>5.7 || Bmass<4.7) continue;
    if (Mll>q2sup || Mll<q2inf) continue;
    if (xgb<0)        continue;
    if (KLmassD0<2.0) continue;

    // HLT and MC-match should be applied already at ntuple level

    // Fill BDT histos w/o weights
    h1mc_xgb_fromLoop_now -> Fill(xgb);

    // Fill BDT histos with weights    
    float thisEventWeight = GetAnBdtWeight(xgb);
    h1mc_xgb_fromLoop_ww -> Fill(xgb, thisEventWeight);

    // Statistics
    totEvents   = totEvents + 1;
    totEventsWW = totEventsWW + thisEventWeight;
    if (xgb>=5) {
      totEventsBdtGt5   = totEventsBdtGt5 + 1;
      totEventsBdtGt5WW = totEventsBdtGt5WW + thisEventWeight;
    }

  } // Loop over entries

  cout << "------------------------------------------" << endl;
  cout << "Bdt>5 / Total = " << (float)totEventsBdtGt5/(float)totEvents << endl;
  cout << "Weighted Bdt>5 / Total = " << (float)totEventsBdtGt5WW/(float)totEventsWW << endl;
  cout << "------------------------------------------" << endl;

  // Cosmetics
  h1mc_xgb_fromLoop_now -> SetLineColor(2);
  h1mc_xgb_fromLoop_ww  -> SetLineColor(3);
  h1mc_xgb_fromFile     -> SetLineColor(4);
  h1data_xgb_fromFile   -> SetMarkerColor(1);
  h1data_xgb_fromFile   -> SetMarkerSize(1);
  h1data_xgb_fromFile   -> SetMarkerStyle(20);

  // Plots
  gStyle->SetOptStat(0);

  TCanvas c1("c1","Sanity check",1);
  h1mc_xgb_fromLoop_now -> DrawNormalized();
  h1mc_xgb_fromFile     -> DrawNormalized("same");
  c1.SaveAs("CheckMC.png");

  TCanvas c2("c2","Data and MC w/o weights",1);
  h1mc_xgb_fromLoop_now -> SetTitle("");
  h1data_xgb_fromFile   -> SetTitle("");
  h1mc_xgb_fromLoop_now -> DrawNormalized();
  h1data_xgb_fromFile   -> DrawNormalized("sameP");
  c2.SaveAs("DataMcNoWeight.png");

  TCanvas c3("c3","Data and MC with weights",1);
  h1mc_xgb_fromLoop_ww -> SetTitle("");
  h1data_xgb_fromFile  -> SetTitle("");
  c3.Divide(1,2);
  c3.cd(1); h1mc_xgb_fromLoop_ww -> DrawNormalized();
  h1data_xgb_fromFile  -> Scale(1./h1data_xgb_fromFile->Integral());
  c3.cd(2); h1data_xgb_fromFile -> Draw();
  c3.SaveAs("DataMcWithWeight.png");

}

float applicationWeight::GetAnBdtWeight(float bdt) {

  int thesizem1 = lowedge_.size()-1;
  
  float weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedge_[i]<=bdt && lowedge_[i+1]>bdt) weight = anweight_[i];  
  }
  
  return weight;
}
