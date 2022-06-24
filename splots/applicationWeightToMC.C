#define applicationWeightToMC_cxx
#include "applicationWeightToMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   
#include <iostream>  

// Macro to be ran on jpsi events to get the same pT(ele1), pT(ele2) or 
// dR(ele1, ele2) distribution as in lowq2 events
// To be ran on MC, to be sure it works

using namespace std;

void applicationWeightToMC::Loop()      
{
  // Define q^2 ranges
  float q2infLowQ2 = 1.05;
  float q2supLowQ2 = 2.45;
  float q2infJPsi  = 2.90;
  float q2supJPsi  = 3.20;

  if (fChain == 0) return;

  // Counters for efficiencies
  float numerator_now = 0.;
  float numerator_wwele1pt = 0.;
  float numerator_wwele2pt = 0.;
  float numerator_wwelesdr = 0.;
  float denominator_now = 0.;
  float denominator_wwele1pt = 0.;
  float denominator_wwele2pt = 0.;
  float denominator_wwelesdr = 0.;

  // File with Ratios
  TFile *ratioFileMCBins = new TFile("outputFiles_PFLP/JPsiToLowQ2BinRatio__PFLP__noAntiD0cut__wp0.0.root");      // ratio = JPsi / LowQ2
  TH1F *ratioEle1pt_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_Ele1Pt-JPsiToLowQ2");
  TH1F *ratioEle2pt_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_Ele2Pt-JPsiToLowQ2");
  TH1F *ratioElesDr_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_ElesDr-JPsiToLowQ2");
  
  // Files with MC distributions
  TFile *mcJPsiBinFile  = new TFile("outputFiles_PFLP/mcForBinRatio_PFLP_JPsiBin__noAntiD0cut__wp0.0.root");
  TFile *mcLowQ2BinFile = new TFile("outputFiles_PFLP/mcForBinRatio_PFLP_LowQ2Bin__noAntiD0cut__wp0.0.root");

  // Histos with MC distributions 
  TH1F *mcJPsiEle1pt_fromFile  = (TH1F*)mcJPsiBinFile ->Get("h1mc_Ele1pt");
  TH1F *mcLowQ2Ele1pt_fromFile = (TH1F*)mcLowQ2BinFile->Get("h1mc_Ele1pt");
  TH1F *mcJPsiEle2pt_fromFile  = (TH1F*)mcJPsiBinFile ->Get("h1mc_Ele2pt");
  TH1F *mcLowQ2Ele2pt_fromFile = (TH1F*)mcLowQ2BinFile->Get("h1mc_Ele2pt");
  TH1F *mcJPsiElesDr_fromFile  = (TH1F*)mcJPsiBinFile ->Get("h1mc_ElesDr");
  TH1F *mcLowQ2ElesDr_fromFile = (TH1F*)mcLowQ2BinFile->Get("h1mc_ElesDr");

  // Histos filled during the loop on jpsi events
  TH1F *mcJPsiEle1pt_fromLoop_now = new TH1F("mcJPsiEle1pt_fromLoop_now","mcJPsiEle1pt_fromLoop_now", 80,  0.0, 40.);
  TH1F *mcJPsiEle1pt_fromLoop_ww  = new TH1F("mcJPsiEle1pt_fromLoop_ww", "mcJPsiEle1pt_fromLoop_ww",  80,  0.0, 40.);
  TH1F *mcJPsiEle2pt_fromLoop_now = new TH1F("mcJPsiEle2pt_fromLoop_now","mcJPsiEle2pt_fromLoop_now", 40,  0.0, 20.);
  TH1F *mcJPsiEle2pt_fromLoop_ww  = new TH1F("mcJPsiEle2pt_fromLoop_ww", "mcJPsiEle2pt_fromLoop_ww",  40,  0.0, 20.);
  TH1F *mcJPsiElesDr_fromLoop_now = new TH1F("mcJPsiElesDr_fromLoop_now","mcJPsiElesDr_fromLoop_now", 20,  0.0, 2.);
  TH1F *mcJPsiElesDr_fromLoop_ww  = new TH1F("mcJPsiElesDr_fromLoop_ww", "mcJPsiElesDr_fromLoop_ww",  20,  0.0, 2.);
  mcJPsiEle1pt_fromLoop_now -> Sumw2();
  mcJPsiEle1pt_fromLoop_ww  -> Sumw2();
  mcJPsiEle2pt_fromLoop_now -> Sumw2();
  mcJPsiEle2pt_fromLoop_ww  -> Sumw2();
  mcJPsiElesDr_fromLoop_now -> Sumw2();
  mcJPsiElesDr_fromLoop_ww  -> Sumw2();

  // Load weights: ele1 pt
  for (int i = 0; i<ratioEle1pt_fromFile->GetNbinsX(); i++) {

    float myweightEle1pt=ratioEle1pt_fromFile->GetBinContent(i+1);
    anweightEle1pt_.push_back(myweightEle1pt);              

    float mylowedgeEle1pt=ratioEle1pt_fromFile->GetBinLowEdge(i+1); 
    lowedgeEle1pt_.push_back(mylowedgeEle1pt);
  }

  // Load weights: ele2 pt
  for (int i = 0; i<ratioEle2pt_fromFile->GetNbinsX(); i++) {

    float myweightEle2pt=ratioEle2pt_fromFile->GetBinContent(i+1);
    anweightEle2pt_.push_back(myweightEle2pt);              

    float mylowedgeEle2pt=ratioEle2pt_fromFile->GetBinLowEdge(i+1); 
    lowedgeEle2pt_.push_back(mylowedgeEle2pt);
  }

  // Load weights: dR(ele1, ele2)
  for (int i = 0; i<ratioElesDr_fromFile->GetNbinsX(); i++) {

    float myweightElesDr=ratioElesDr_fromFile->GetBinContent(i+1);
    anweightElesDr_.push_back(myweightElesDr);              

    float mylowedgeElesDr=ratioElesDr_fromFile->GetBinLowEdge(i+1); 
    lowedgeElesDr_.push_back(mylowedgeElesDr);
  }


  // Loop over entries: this is the jpsi file
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "nentries = " << nentries << endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    // JPsi Q2 range
    if (Mll>q2supJPsi || Mll<q2infJPsi) continue;

    // Same selection as applied in data
    if (Bmass>5.7 || Bmass<4.7) continue;
    if (xgb<0)        continue;
    // if (KLmassD0<2.0) continue;

   
    // Fill histos w/o weights
    mcJPsiEle1pt_fromLoop_now -> Fill(L1pt);
    mcJPsiEle2pt_fromLoop_now -> Fill(L2pt);
    mcJPsiElesDr_fromLoop_now -> Fill(L1L2dr);

    // Fill histos with weights    
    float thisEventWeightEle1pt = 1./GetEle1ptWeight(L1pt);
    mcJPsiEle1pt_fromLoop_ww -> Fill(L1pt, thisEventWeightEle1pt);

    float thisEventWeightEle2pt = 1./GetEle2ptWeight(L2pt);
    mcJPsiEle2pt_fromLoop_ww -> Fill(L2pt, thisEventWeightEle2pt);

    float thisEventWeightElesDr = 1./GetElesDrWeight(L1L2dr);
    if (GetElesDrWeight(L1L2dr)!=0) {
      mcJPsiElesDr_fromLoop_ww -> Fill(L1L2dr, thisEventWeightElesDr);
    }

    // To compute efficiencies: denominator
    denominator_now = denominator_now+1;
    denominator_wwele1pt  = denominator_wwele1pt+thisEventWeightEle1pt;
    denominator_wwele2pt  = denominator_wwele2pt+thisEventWeightEle2pt;
    if (GetElesDrWeight(L1L2dr)!=0) denominator_wwelesdr = denominator_wwelesdr+thisEventWeightElesDr;

    // To compute efficiencies: numerator
    if (xgb<8.3) continue;
    if (KLmassD0<2.0) continue;
    numerator_now = numerator_now+1;
    numerator_wwele1pt  = numerator_wwele1pt+thisEventWeightEle1pt;
    numerator_wwele2pt  = numerator_wwele2pt+thisEventWeightEle2pt;
    if (GetElesDrWeight(L1L2dr)!=0) numerator_wwelesdr = numerator_wwelesdr+thisEventWeightElesDr;


  } // Loop over entries


  // Statistics / efficiencies
  cout << endl;
  cout << "---------------------------------------" << endl;
  cout << "No weight: Numerator = "        << numerator_now      << ", denominator = " << denominator_now      << " => eff = " << numerator_now/denominator_now << endl;
  cout << "Ele1 pt weight: Numerator = "   << numerator_wwele1pt << ", denominator = " << denominator_wwele1pt << " => eff = " << numerator_wwele1pt/denominator_wwele1pt << endl;
  cout << "Ele2 pt weight: Numerator = "   << numerator_wwele2pt << ", denominator = " << denominator_wwele2pt << " => eff = " << numerator_wwele2pt/denominator_wwele2pt << endl;
  cout << "Dr(e1,e2) weight: Numerator = " << numerator_wwelesdr << ", denominator = " << denominator_wwelesdr << " => eff = " << numerator_wwelesdr/denominator_wwelesdr << endl;
  cout << "---------------------------------------" << endl;
  cout << endl;


  // Cosmetics
  mcJPsiEle1pt_fromLoop_now -> SetLineColor(2);  
  mcJPsiEle1pt_fromLoop_ww  -> SetLineColor(1);  
  mcJPsiEle2pt_fromLoop_now -> SetLineColor(2);  
  mcJPsiEle2pt_fromLoop_ww  -> SetLineColor(1);  
  mcJPsiElesDr_fromLoop_now -> SetLineColor(2);  
  mcJPsiElesDr_fromLoop_ww  -> SetLineColor(1);  
  mcJPsiEle1pt_fromFile     -> SetLineColor(6);  
  mcLowQ2Ele1pt_fromFile    -> SetLineColor(4);  
  mcJPsiEle2pt_fromFile     -> SetLineColor(6);  
  mcLowQ2Ele2pt_fromFile    -> SetLineColor(4);  
  mcJPsiElesDr_fromFile     -> SetLineColor(6);  
  mcLowQ2ElesDr_fromFile    -> SetLineColor(4);  
  mcJPsiEle1pt_fromLoop_now -> SetLineWidth(2);  
  mcJPsiEle1pt_fromLoop_ww  -> SetLineWidth(2);  
  mcJPsiEle2pt_fromLoop_now -> SetLineWidth(2);  
  mcJPsiEle2pt_fromLoop_ww  -> SetLineWidth(2);  
  mcJPsiElesDr_fromLoop_now -> SetLineWidth(2);  
  mcJPsiElesDr_fromLoop_ww  -> SetLineWidth(2);  
  mcJPsiEle1pt_fromFile     -> SetLineWidth(2);  
  mcLowQ2Ele1pt_fromFile    -> SetLineWidth(2);  
  mcJPsiEle2pt_fromFile     -> SetLineWidth(2);  
  mcLowQ2Ele2pt_fromFile    -> SetLineWidth(2);  
  mcJPsiElesDr_fromFile     -> SetLineWidth(2);  
  mcLowQ2ElesDr_fromFile    -> SetLineWidth(2);  

  TLegend leg2 (0.55,0.55,0.90,0.85);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(mcJPsiEle1pt_fromLoop_now, "MC, JPsi bin","PLE");
  leg2.AddEntry(mcJPsiEle1pt_fromLoop_ww,  "MC, JPsi bin. Weighted","L");
  leg2.AddEntry(mcLowQ2Ele1pt_fromFile,    "MC, Low-q2 bin","L");

  // Plots
  gStyle->SetOptStat(0);

  TCanvas c1a("c1a","Sanity check",1);
  mcJPsiEle1pt_fromFile     -> DrawNormalized("hist");
  mcJPsiEle1pt_fromLoop_now -> DrawNormalized("samehist");
  c1a.SaveAs("CheckMcEle1pt.png");

  TCanvas c1b("c1b","Sanity check",1);
  mcJPsiEle2pt_fromFile     -> DrawNormalized("hist");
  mcJPsiEle2pt_fromLoop_now -> DrawNormalized("samehist");
  c1b.SaveAs("CheckMcEle2pt.png");

  TCanvas c1c("c1c","Sanity check",1);
  mcJPsiElesDr_fromFile     -> DrawNormalized("hist");
  mcJPsiElesDr_fromLoop_now -> DrawNormalized("samehist");
  c1c.SaveAs("CheckMcElesDr.png");

  TCanvas c2a("c2a","MC with and wo weights",1);
  mcLowQ2Ele1pt_fromFile    -> Rebin(4);
  mcJPsiEle1pt_fromLoop_now -> Rebin(4);
  mcJPsiEle1pt_fromLoop_ww  -> Rebin(4);
  mcLowQ2Ele1pt_fromFile    -> DrawNormalized("hist");
  mcJPsiEle1pt_fromLoop_now -> DrawNormalized("samehist");
  mcJPsiEle1pt_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c2a.SaveAs("Ele1pt_withAndWoWeight.png");

  TCanvas c2b("c2b","MC with and wo weights",1);
  mcLowQ2Ele2pt_fromFile    -> Rebin();
  mcJPsiEle2pt_fromLoop_now -> Rebin();
  mcJPsiEle2pt_fromLoop_ww  -> Rebin();
  mcLowQ2Ele2pt_fromFile    -> DrawNormalized("hist");
  mcJPsiEle2pt_fromLoop_now -> DrawNormalized("samehist");
  mcJPsiEle2pt_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c2b.SaveAs("Ele2pt_withAndWoWeight.png");

  TCanvas c2c("c2c","MC with and wo weights",1);
  mcLowQ2ElesDr_fromFile    -> DrawNormalized("hist");
  mcJPsiElesDr_fromLoop_now -> DrawNormalized("samehist");
  mcJPsiElesDr_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c2c.SaveAs("ElesDr_withAndWoWeight.png");
}


float applicationWeightToMC::GetEle1ptWeight(float ele1pt) {

  int thesizem1 = lowedgeEle1pt_.size()-1;
  
  float weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedgeEle1pt_[i]<=ele1pt && lowedgeEle1pt_[i+1]>ele1pt) weight = anweightEle1pt_[i];  
  }
  
  return weight;
}

float applicationWeightToMC::GetEle2ptWeight(float ele2pt) {

  int thesizem1 = lowedgeEle2pt_.size()-1;
  
  float weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedgeEle2pt_[i]<=ele2pt && lowedgeEle2pt_[i+1]>ele2pt) weight = anweightEle2pt_[i];  
  }
  
  return weight;
}

float applicationWeightToMC::GetElesDrWeight(float elesdr) {

  int thesizem1 = lowedgeElesDr_.size()-1;
  
  float weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedgeElesDr_[i]<=elesdr && lowedgeElesDr_[i+1]>elesdr) weight = anweightElesDr_[i];  
  }
  
  return weight;
}



