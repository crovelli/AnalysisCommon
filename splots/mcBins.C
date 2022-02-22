#define mcBins_cxx

// Compare MC signal distributions for different q2 bins

// ROOT includes 
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <iostream>

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

  TLegend leg2 (0.15,0.65,0.45,0.85);
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
  frame->GetYaxis()->SetRangeUser(0.,4.);
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
  h2->Write( ("MC_" + canvasName).c_str());
  fileOut.Close();
}

void mcBins()
{
  // Input files 
  TFile *fileJpsi  = new TFile("files_nonRegr/myFileMC_PFLP_JPsi.root");
  TFile *fileLowQ2 = new TFile("files_nonRegr/myFileMC_PFLP_LowQ2.root");
  TFile *filePsi2s = new TFile("files_nonRegr/myFileMC_PFLP_Psi2s.root");
  
  // Histos 
  TH1F *h1jpsi_xgb  = (TH1F*)fileJpsi ->Get("h1mc_xgb");
  TH1F *h1lowq2_xgb = (TH1F*)fileLowQ2->Get("h1mc_xgb");
  TH1F *h1psi2s_xgb = (TH1F*)filePsi2s->Get("h1mc_xgb");
  //
  TH1F *h1jpsi_mll  = (TH1F*)fileJpsi ->Get("h1mc_mll");
  TH1F *h1lowq2_mll = (TH1F*)fileLowQ2->Get("h1mc_mll");
  TH1F *h1psi2s_mll = (TH1F*)filePsi2s->Get("h1mc_mll");


  // Rebin
  h1jpsi_xgb->Rebin(3);
  h1lowq2_xgb->Rebin(3);
  h1psi2s_xgb->Rebin(3);

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  // 
  drawTH1pair(h1jpsi_xgb,   h1lowq2_xgb,  "BDT-JPsiToLowQ2",  "a.u.",  1.,"bdt-JPsiToLowQ2",  "./",2, 0, "JPsi","LowQ2");
  drawTH1pair(h1jpsi_xgb,   h1psi2s_xgb,  "BDT-JPsiToPsi2s",  "a.u.",  1.,"bdt-JPsiToPsi2s",  "./",2, 0, "JPsi","Psi2s");
  drawTH1pair(h1jpsi_mll,   h1lowq2_mll,  "Mll-JPsiToLowQ2",  "a.u.",  1.,"mll-JPsiToLowQ2",  "./",2, 0, "JPsi","LowQ2");
  drawTH1pair(h1jpsi_mll,   h1psi2s_mll,  "Mll-JPsiToPsi2s",  "a.u.",  1.,"mll-JPsiToPsi2s",  "./",2, 0, "JPsi","Psi2s");
}
