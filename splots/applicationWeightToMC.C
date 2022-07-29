#define applicationWeightToMC_cxx
#include "applicationWeightToMC.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   
#include <iostream>  

// Macro to be ran on jpsi events to get the same pT(ele1), pT(ele2) or 
// dR(ele1, ele2) distribution as in lowq2 events
// To be ran on MC, to be sure it works

using namespace std;

void drawTH1pair(TH1* h1, TH1* h2, 
		 const string& xAxisNameTmp = "", const string& yAxisName = "Events", 
		 float lumi=-1, const string& canvasName = "default", 
		 const string& outputDIR = "./", int mycolor=2, int logy=0,
		 const string& legEntry1 = "data", const string& legEntry2 = "MC", const string& ratioPadYaxisName = "data/MC", const string& outputFILE = "outFile.root") 
{
  string xAxisName = "";
  string separator = "::";
  Bool_t setXAxisRangeFromUser = false;
  Double_t xmin = 0;
  Double_t xmax = 0;

  size_t pos = xAxisNameTmp.find(separator);
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    xAxisName.assign(xAxisNameTmp, 0, pos); 
    xrange.assign(xAxisNameTmp, pos + separator.size(), string::npos);
    separator = ",";
    pos = xrange.find(separator);
    string numString = ""; 
    numString.assign(xrange,0,pos);
    xmin = std::stod(numString);
    numString.assign(xrange,pos + separator.size(), string::npos);
    xmax = std::stod(numString);
  } else {
    xAxisName = xAxisNameTmp;
  }

  if (yAxisName == "a.u.") {
    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
  }
  else if (lumi>-1) {
    h1->Scale(lumi/h1->Integral());
    h2->Scale(lumi/h2->Integral());
  }

  h1->SetStats(0);
  h2->SetStats(0);

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h1->SetTitle("");
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitle(yAxisName.c_str());
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetMaximum(),h2->GetMaximum()) * 1.2);
  if (setXAxisRangeFromUser) h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->Draw("hist E");

  h2->SetTitle("");
  h2->SetLineColor(mycolor);
  h2->SetLineWidth(2);
  h2->Draw("hist E same");

  //TLegend leg2 (0.15,0.65,0.45,0.85);
  TLegend leg2 (0.65,0.65,0.95,0.85);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(h1,legEntry1.c_str(),"L");
  leg2.AddEntry(h2,legEntry2.c_str(),"L");

  leg2.Draw("same");
  canvas->RedrawAxis("sameaxis");

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.5,1.5);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
  if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
  frame->GetXaxis()->SetTitleSize(0.05);
  
  TH1D* ratio = (TH1D*) h1->Clone("ratio");
  TH1D* den_noerr = (TH1D*) h2->Clone("den_noerr");
  TH1D* den = (TH1D*) h2->Clone("den");
  for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++) den_noerr->SetBinError(iBin,0.);

  ratio->Divide(den_noerr);
  den->Divide(den_noerr);
  den->SetFillColor(kGray);
  frame->Draw();
  ratio->SetMarkerSize(0.85);
  ratio->Draw("EPsame");
  den->Draw("E2same");

  TF1* line = new TF1("horiz_line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(mycolor);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  canvas->SetLogy(logy);
  canvas->SaveAs((outputDIR + canvasName + ".png").c_str());

  if (yAxisName == "a.u.") h1->GetYaxis()->SetRangeUser(max(0.0001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
  else h1->GetYaxis()->SetRangeUser(max(0.001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);

  delete canvas;
  frame->Reset("ICES");

  TFile fileOut(outputFILE.c_str(), "UPDATE");
  fileOut.cd();
  ratio->Write( ("ratio_" + canvasName).c_str());
  h1->Write( ("MC1_" + canvasName).c_str());
  h2->Write( ("MC2_" + canvasName).c_str());
  fileOut.Close();
}

void applicationWeightToMC::Loop(int theVariable)       
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
  float numerator_wwkpt = 0.;
  float numerator_wwelesdr = 0.;
  float denominator_now = 0.;
  float denominator_wwele1pt = 0.;
  float denominator_wwele2pt = 0.;
  float denominator_wwkpt = 0.;
  float denominator_wwelesdr = 0.;

  // File with Ratios
  TFile *ratioFileMCBins = new TFile("JPsiToLowQ2BinRatio__PFPF_withAntiD0cut__wp0.0.root");      // ratio = JPsi / LowQ2
  TH1F *ratioEle1pt_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_Ele1Pt-JPsiToLowQ2");
  TH1F *ratioEle2pt_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_Ele2Pt-JPsiToLowQ2");
  TH1F *ratioKpt_fromFile    = (TH1F*)ratioFileMCBins->Get("ratio_KPt-JPsiToLowQ2");
  TH1F *ratioElesDr_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_ElesDr-JPsiToLowQ2");
  
  // Files with MC distributions
  TFile *mcJPsiBinFile  = new TFile("mcForBinRatio_PFPF_JPsiBin__withAntiD0cut__wp0.0.root");
  TFile *mcLowQ2BinFile = new TFile("mcForBinRatio_PFPF_LowQ2Bin__withAntiD0cut__wp0.0.root");

  // Histos with MC distributions 
  TH1F *mcJPsiEle1pt_fromFile    = (TH1F*)mcJPsiBinFile ->Get("h1mc_Ele1pt");
  TH1F *mcLowQ2Ele1pt_fromFile   = (TH1F*)mcLowQ2BinFile->Get("h1mc_Ele1pt");
  TH1F *mcJPsiEle2pt_fromFile    = (TH1F*)mcJPsiBinFile ->Get("h1mc_Ele2pt");
  TH1F *mcLowQ2Ele2pt_fromFile   = (TH1F*)mcLowQ2BinFile->Get("h1mc_Ele2pt");
  TH1F *mcJPsiKpt_fromFile       = (TH1F*)mcJPsiBinFile ->Get("h1mc_Kpt");
  TH1F *mcLowQ2Kpt_fromFile      = (TH1F*)mcLowQ2BinFile->Get("h1mc_Kpt");
  TH1F *mcJPsiElesDr_fromFile    = (TH1F*)mcJPsiBinFile ->Get("h1mc_ElesDr");
  TH1F *mcLowQ2ElesDr_fromFile   = (TH1F*)mcLowQ2BinFile->Get("h1mc_ElesDr");
  TH1F *mcJPsiXgb_fromFile       = (TH1F*)mcJPsiBinFile ->Get("h1mc_xgb");
  TH1F *mcLowQ2Xgb_fromFile      = (TH1F*)mcLowQ2BinFile->Get("h1mc_xgb");
  TH1F *mcJPsiKLmassD0_fromFile  = (TH1F*)mcJPsiBinFile ->Get("h1mc_KLmassD0");
  TH1F *mcLowQ2KLmassD0_fromFile = (TH1F*)mcLowQ2BinFile->Get("h1mc_KLmassD0");
  
  // Histos filled during the loop on jpsi events
  TH1F *mcJPsiEle1pt_fromLoop_now = new TH1F("mcJPsiEle1pt_fromLoop_now","mcJPsiEle1pt_fromLoop_now", 80,  0.0, 40.);
  TH1F *mcJPsiEle1pt_fromLoop_ww  = new TH1F("mcJPsiEle1pt_fromLoop_ww", "mcJPsiEle1pt_fromLoop_ww",  80,  0.0, 40.);
  TH1F *mcJPsiEle2pt_fromLoop_now = new TH1F("mcJPsiEle2pt_fromLoop_now","mcJPsiEle2pt_fromLoop_now", 40,  0.0, 20.);
  TH1F *mcJPsiEle2pt_fromLoop_ww  = new TH1F("mcJPsiEle2pt_fromLoop_ww", "mcJPsiEle2pt_fromLoop_ww",  40,  0.0, 20.);
  TH1F *mcJPsiKpt_fromLoop_now    = new TH1F("mcJPsiKpt_fromLoop_now",   "mcJPsiKpt_fromLoop_now",  40,  0.0, 20.);
  TH1F *mcJPsiKpt_fromLoop_ww     = new TH1F("mcJPsiKpt_fromLoop_ww",    "mcJPsiKpt_fromLoop_ww",   40,  0.0, 20.);
  TH1F *mcJPsiElesDr_fromLoop_now = new TH1F("mcJPsiElesDr_fromLoop_now","mcJPsiElesDr_fromLoop_now", 20,  0.0, 2.);
  TH1F *mcJPsiElesDr_fromLoop_ww  = new TH1F("mcJPsiElesDr_fromLoop_ww", "mcJPsiElesDr_fromLoop_ww",  20,  0.0, 2.);
  TH1F *mcJPsiXgb_fromLoop_now      = new TH1F("mcJPsiXgb_fromLoop_now", "mcJPsiXgb_fromLoop_now", 45,  0.0, 15.);
  TH1F *mcJPsiXgb_fromLoop_ww       = new TH1F("mcJPsiXgb_fromLoop_ww",  "mcJPsiXgb_fromLoop_ww",  45,  0.0, 15.);
  TH1F *mcJPsiKLmassD0_fromLoop_now = new TH1F("mcJPsiKLmassD0_fromLoop_now", "mcJPsiKLmassD0_fromLoop_now", 30,  0.0, 6.);
  TH1F *mcJPsiKLmassD0_fromLoop_ww  = new TH1F("mcJPsiKLmassD0_fromLoop_ww",  "mcJPsiKLmassD0_fromLoop_ww",  30,  0.0, 6.);
  //
  mcJPsiEle1pt_fromLoop_now -> Sumw2();
  mcJPsiEle1pt_fromLoop_ww  -> Sumw2();
  mcJPsiEle2pt_fromLoop_now -> Sumw2();
  mcJPsiEle2pt_fromLoop_ww  -> Sumw2();
  mcJPsiKpt_fromLoop_now    -> Sumw2();
  mcJPsiKpt_fromLoop_ww     -> Sumw2();
  mcJPsiElesDr_fromLoop_now -> Sumw2();
  mcJPsiElesDr_fromLoop_ww  -> Sumw2();
  mcJPsiXgb_fromLoop_now    -> Sumw2();
  mcJPsiXgb_fromLoop_ww     -> Sumw2();
  mcJPsiKLmassD0_fromLoop_now -> Sumw2();
  mcJPsiKLmassD0_fromLoop_ww  -> Sumw2();

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
    if (KLmassD0<2.0) continue;

    double theWeight;
    if (theVariable==0) theWeight = weightElePt1;
    if (theVariable==1) theWeight = weightElePt2;
    if (theVariable==2) theWeight = weightKPt;
    if (theVariable==3) theWeight = weightElesDr;
   
    // Fill histos w/o weights
    mcJPsiEle1pt_fromLoop_now   -> Fill(L1pt);
    mcJPsiEle2pt_fromLoop_now   -> Fill(L2pt);
    mcJPsiKpt_fromLoop_now      -> Fill(Kpt);
    mcJPsiElesDr_fromLoop_now   -> Fill(L1L2dr);
    mcJPsiXgb_fromLoop_now      -> Fill(xgb);
    mcJPsiKLmassD0_fromLoop_now -> Fill(KLmassD0);
    mcJPsiEle1pt_fromLoop_ww    -> Fill(L1pt, weightElePt1);
    mcJPsiEle2pt_fromLoop_ww    -> Fill(L2pt, weightElePt2);
    mcJPsiKpt_fromLoop_ww       -> Fill(Kpt, weightKPt);
    mcJPsiElesDr_fromLoop_ww    -> Fill(L1L2dr, weightElesDr);
    mcJPsiXgb_fromLoop_ww       -> Fill(xgb, theWeight);
    mcJPsiKLmassD0_fromLoop_ww  -> Fill(KLmassD0, theWeight);

     // To compute efficiencies: denominator
    denominator_now = denominator_now+1;
    denominator_wwele1pt  = denominator_wwele1pt+weightElePt1;
    denominator_wwele2pt  = denominator_wwele2pt+weightElePt2;
    denominator_wwkpt     = denominator_wwkpt+weightKPt;
    denominator_wwelesdr  = denominator_wwelesdr+weightElesDr;

    // To compute efficiencies: numerator
    if (xgb<8.3) continue;
    if (KLmassD0<2.0) continue;
    numerator_now = numerator_now+1;
    numerator_wwele1pt  = numerator_wwele1pt+weightElePt1;
    numerator_wwele2pt  = numerator_wwele2pt+weightElePt2;
    numerator_wwkpt     = numerator_wwkpt+weightKPt;
    numerator_wwelesdr  = numerator_wwelesdr+weightElesDr;


  } // Loop over entries


  // Statistics / efficiencies
  cout << endl;
  cout << "---------------------------------------" << endl;
  cout << "No weight: Numerator = "        << numerator_now      << ", denominator = " << denominator_now      << " => eff = " << numerator_now/denominator_now << endl;
  cout << "Ele1 pt weight: Numerator = "   << numerator_wwele1pt << ", denominator = " << denominator_wwele1pt << " => eff = " << numerator_wwele1pt/denominator_wwele1pt << endl;
  cout << "Ele2 pt weight: Numerator = "   << numerator_wwele2pt << ", denominator = " << denominator_wwele2pt << " => eff = " << numerator_wwele2pt/denominator_wwele2pt << endl;
  cout << "K pt weight: Numerator = "      << numerator_wwkpt    << ", denominator = " << denominator_wwkpt    << " => eff = " << numerator_wwkpt/denominator_wwkpt << endl;
  cout << "Dr(e1,e2) weight: Numerator = " << numerator_wwelesdr << ", denominator = " << denominator_wwelesdr << " => eff = " << numerator_wwelesdr/denominator_wwelesdr << endl;
  cout << "---------------------------------------" << endl;
  cout << endl;


  // Cosmetics
  mcJPsiEle1pt_fromLoop_now -> SetLineColor(2);  
  mcJPsiEle1pt_fromLoop_ww  -> SetLineColor(1);  
  mcJPsiEle2pt_fromLoop_now -> SetLineColor(2);  
  mcJPsiEle2pt_fromLoop_ww  -> SetLineColor(1);  
  mcJPsiKpt_fromLoop_now    -> SetLineColor(2);  
  mcJPsiKpt_fromLoop_ww     -> SetLineColor(1);  
  mcJPsiElesDr_fromLoop_now -> SetLineColor(2);  
  mcJPsiElesDr_fromLoop_ww  -> SetLineColor(1);  
  mcJPsiXgb_fromLoop_now      -> SetLineColor(2);
  mcJPsiXgb_fromLoop_ww       -> SetLineColor(1);
  mcJPsiKLmassD0_fromLoop_now -> SetLineColor(2);
  mcJPsiKLmassD0_fromLoop_ww  -> SetLineColor(1);
  mcJPsiEle1pt_fromFile     -> SetLineColor(6);  
  mcLowQ2Ele1pt_fromFile    -> SetLineColor(4);  
  mcJPsiEle2pt_fromFile     -> SetLineColor(6);  
  mcLowQ2Ele2pt_fromFile    -> SetLineColor(4);  
  mcJPsiKpt_fromFile        -> SetLineColor(6);  
  mcLowQ2Kpt_fromFile       -> SetLineColor(4);  
  mcJPsiElesDr_fromFile     -> SetLineColor(6);  
  mcLowQ2ElesDr_fromFile    -> SetLineColor(4);  
  mcJPsiEle1pt_fromLoop_now -> SetLineWidth(2);  
  mcJPsiEle1pt_fromLoop_ww  -> SetLineWidth(2);  
  mcJPsiEle2pt_fromLoop_now -> SetLineWidth(2);  
  mcJPsiEle2pt_fromLoop_ww  -> SetLineWidth(2);  
  mcJPsiKpt_fromLoop_now    -> SetLineWidth(2);  
  mcJPsiKpt_fromLoop_ww     -> SetLineWidth(2);  
  mcJPsiElesDr_fromLoop_now -> SetLineWidth(2);  
  mcJPsiElesDr_fromLoop_ww  -> SetLineWidth(2);  
  mcJPsiXgb_fromLoop_now      -> SetLineWidth(2);
  mcJPsiXgb_fromLoop_ww       -> SetLineWidth(2);
  mcJPsiKLmassD0_fromLoop_now -> SetLineWidth(2);
  mcJPsiKLmassD0_fromLoop_ww  -> SetLineWidth(2);
  mcJPsiEle1pt_fromFile     -> SetLineWidth(2);  
  mcLowQ2Ele1pt_fromFile    -> SetLineWidth(2);  
  mcJPsiEle2pt_fromFile     -> SetLineWidth(2);  
  mcLowQ2Ele2pt_fromFile    -> SetLineWidth(2);  
  mcJPsiKpt_fromFile        -> SetLineWidth(2);  
  mcLowQ2Kpt_fromFile       -> SetLineWidth(2);  
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

  TCanvas c1d("c1d","Sanity check",1);
  mcJPsiKpt_fromFile     -> DrawNormalized("hist");
  mcJPsiKpt_fromLoop_now -> DrawNormalized("samehist");
  c1d.SaveAs("CheckMcKpt.png");

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

  TCanvas c2d("c2d","MC with and wo weights",1);
  mcLowQ2Kpt_fromFile    -> Rebin();
  mcJPsiKpt_fromLoop_now -> Rebin();
  mcJPsiKpt_fromLoop_ww  -> Rebin();
  mcJPsiKpt_fromLoop_now -> DrawNormalized("hist");
  mcLowQ2Kpt_fromFile    -> DrawNormalized("samehist");
  mcJPsiKpt_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c2d.SaveAs("Kpt_withAndWoWeight.png");

  TCanvas c2c("c2c","MC with and wo weights",1);
  mcLowQ2ElesDr_fromFile    -> DrawNormalized("hist");
  mcJPsiElesDr_fromLoop_now -> DrawNormalized("samehist");
  mcJPsiElesDr_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c2c.SaveAs("ElesDr_withAndWoWeight.png");

  TCanvas c3a("c3a","MC with and wo weights",1);
  mcLowQ2Xgb_fromFile    -> DrawNormalized("hist");
  mcJPsiXgb_fromLoop_now -> DrawNormalized("samehist");
  mcJPsiXgb_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c3a.SaveAs("Xgb_withAndWoWeight.png");
  c3a.SetLogy();
  c3a.SaveAs("Xgb_withAndWoWeightLog.png");

  drawTH1pair(mcJPsiXgb_fromLoop_now, mcLowQ2Xgb_fromFile,  "BDT-JPsiToLowQ2 NoWeight",  "a.u.",  1.,"bdt-JPsiToLowQ2-noWeight",  "./",2, 0, "JPsi","LowQ2");
  drawTH1pair(mcJPsiXgb_fromLoop_ww, mcLowQ2Xgb_fromFile,  "BDT-JPsiToLowQ2 WithWeight",  "a.u.",  1.,"bdt-JPsiToLowQ2-withWeight",  "./",2, 0, "JPsi","LowQ2");

  TCanvas c3b("c3b","MC with and wo weights",1);
  mcLowQ2KLmassD0_fromFile    -> DrawNormalized("hist");
  mcJPsiKLmassD0_fromLoop_now -> DrawNormalized("samehist");
  mcJPsiKLmassD0_fromLoop_ww  -> DrawNormalized("samehist");
  leg2.Draw();
  c3b.SaveAs("KLmassD0_withAndWoWeight.png");
  c3b.SetLogy();
  c3b.SaveAs("KLmassD0_withAndWoWeightLog.png");
}



