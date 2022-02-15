#define mcVsSplots_cxx

// Data vs MC signal distributions

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

  // To deal with splots
  for (int ii=0; ii<h1->GetNbinsX(); ii++){
    int iip1 = ii+1;
    if (h1->GetBinContent(iip1)<=0) { 
      h1->SetBinContent(iip1,0);
      h1->SetBinError(iip1,0);
    } 
    if (h2->GetBinContent(iip1)<=0) {
      h2->SetBinContent(iip1,0);
      h2->SetBinError(iip1,0);
    } 
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
  h1->Draw("EP");

  h2->SetTitle("");
  h2->SetLineColor(mycolor);
  h2->SetLineWidth(2);
  h2->Draw("hist E same");

  TLegend leg2 (0.65,0.65,0.95,0.85);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(h1,legEntry1.c_str(),"PLE");
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

void mcVsSplots(int q2bin)
{
  // Input files 
  TFile *fileMC = new TFile("./files_nonRegr/mc_PFPF_JPsiBin__wp5.0.root");
  TFile *fileSPlots = new TFile("./files_nonRegr/splots_PFPF_JPsiBin___testA__wp5.0.root");
  
  // MC histos 
  TH1F *h1mc_xgb   = (TH1F*)fileMC->Get("h1mc_xgb");
  TH1F *h1mc_L1id  = (TH1F*)fileMC->Get("h1mc_L1id");
  TH1F *h1mc_L2id  = (TH1F*)fileMC->Get("h1mc_L2id");
  TH1F *h1mc_Bprob = (TH1F*)fileMC->Get("h1mc_Bprob");
  TH1F *h1mc_BsLxy = (TH1F*)fileMC->Get("h1mc_BsLxy");
  TH1F *h1mc_Bcos  = (TH1F*)fileMC->Get("h1mc_Bcos");
  TH1F *h1mc_L1pt  = (TH1F*)fileMC->Get("h1mc_L1pt");
  TH1F *h1mc_L2pt  = (TH1F*)fileMC->Get("h1mc_L2pt");
  TH1F *h1mc_Bpt   = (TH1F*)fileMC->Get("h1mc_Bpt");
  TH1F *h1mc_Kpt    = (TH1F*)fileMC->Get("h1mc_Kpt");
  TH1F *h1mc_LKdz   = (TH1F*)fileMC->Get("h1mc_LKdz");
  TH1F *h1mc_L1L2dr = (TH1F*)fileMC->Get("h1mc_L1L2dr");
  TH1F *h1mc_LKdr   = (TH1F*)fileMC->Get("h1mc_LKdr");
  TH1F *h1mc_L1iso  = (TH1F*)fileMC->Get("h1mc_L1iso");
  TH1F *h1mc_Kiso   = (TH1F*)fileMC->Get("h1mc_Kiso");
  TH1F *h1mc_Passymetry = (TH1F*)fileMC->Get("h1mc_Passymetry");
  TH1F *h1mc_Kip3d = (TH1F*)fileMC->Get("h1mc_Kip3d");
  TH1F *h1mc_KLmassD0 = (TH1F*)fileMC->Get("h1mc_KLmassD0");

  // s-Plots output 
  TH1F *h1data_xgb   = (TH1F*)fileSPlots->Get("h1_xgb__xgb");
  TH1F *h1data_L1id  = (TH1F*)fileSPlots->Get("h1_L1id__L1id");
  TH1F *h1data_L2id  = (TH1F*)fileSPlots->Get("h1_L2id__L2id");
  TH1F *h1data_Bprob = (TH1F*)fileSPlots->Get("h1_Bprob__Bprob");
  TH1F *h1data_BsLxy = (TH1F*)fileSPlots->Get("h1_BsLxy__BsLxy");
  TH1F *h1data_Bcos  = (TH1F*)fileSPlots->Get("h1_Bcos__Bcos");
  TH1F *h1data_L1pt  = (TH1F*)fileSPlots->Get("h1_L1pt__L1pt");
  TH1F *h1data_L2pt  = (TH1F*)fileSPlots->Get("h1_L2pt__L2pt");
  TH1F *h1data_Bpt   = (TH1F*)fileSPlots->Get("h1_Bpt__Bpt");
  TH1F *h1data_Kpt    = (TH1F*)fileSPlots->Get("h1_Kpt__Kpt");
  TH1F *h1data_LKdz   = (TH1F*)fileSPlots->Get("h1_LKdz__LKdz");
  TH1F *h1data_L1L2dr = (TH1F*)fileSPlots->Get("h1_L1L2dr__L1L2dr");
  TH1F *h1data_LKdr   = (TH1F*)fileSPlots->Get("h1_LKdr__LKdr");
  TH1F *h1data_L1iso  = (TH1F*)fileSPlots->Get("h1_L1iso__L1iso");
  TH1F *h1data_Kiso   = (TH1F*)fileSPlots->Get("h1_Kiso__Kiso");
  TH1F *h1data_Passymetry = (TH1F*)fileSPlots->Get("h1_Passymetry__Passymetry");
  TH1F *h1data_Kip3d = (TH1F*)fileSPlots->Get("h1_Kip3d__Kip3d");
  TH1F *h1data_KLmassD0 = (TH1F*)fileSPlots->Get("h1_KLmassD0__KLmassD0");

  // Rebin
  h1mc_L1pt->Rebin();
  h1mc_L2pt->Rebin();
  h1mc_Bpt->Rebin();
  h1mc_L1iso->Rebin();
  h1mc_Kiso->Rebin();
  h1data_L1pt->Rebin();
  h1data_L2pt->Rebin();
  h1data_Bpt->Rebin();
  h1data_L1iso->Rebin();
  h1data_Kiso->Rebin();
  h1mc_KLmassD0->Rebin();
  h1data_KLmassD0->Rebin();

  if (q2bin==2) {
    h1mc_Kiso->Rebin();
    h1data_Kiso->Rebin();
    h1data_Kpt->Rebin();
    h1mc_Kpt->Rebin();
    h1data_L1pt->Rebin();
    h1mc_L1pt->Rebin();
    h1mc_L1id->Rebin();
    h1data_L1id->Rebin();
    h1data_L2pt->Rebin();
    h1mc_L2pt->Rebin();
    h1mc_L2id->Rebin();
    h1data_L2id->Rebin();
    h1data_LKdr->Rebin();
    h1mc_LKdr->Rebin();
  }

  if (q2bin==0 || q2bin==2) {
    h1data_xgb->Rebin();
    h1mc_xgb->Rebin();
    h1data_Kpt->Rebin();
    h1mc_Kpt->Rebin();
    h1data_LKdr->Rebin();
    h1mc_LKdr->Rebin();
  }

  // In sPlots, put bins with <0 weight to zero
  int nBinsData_xgb   = h1data_xgb->GetNbinsX();
  int nBinsData_L1id  = h1data_L1id->GetNbinsX();
  int nBinsData_L2id  = h1data_L2id->GetNbinsX();
  int nBinsData_Bprob = h1data_Bprob->GetNbinsX();
  int nBinsData_BsLxy = h1data_BsLxy->GetNbinsX();
  int nBinsData_Bcos  = h1data_Bcos->GetNbinsX();
  int nBinsData_L1pt  = h1data_L1pt->GetNbinsX();
  int nBinsData_L2pt  = h1data_L2pt->GetNbinsX();
  int nBinsData_Bpt   = h1data_Bpt->GetNbinsX();
  int nBinsData_Kpt   = h1data_Kpt->GetNbinsX();
  int nBinsData_LKdz  = h1data_LKdz ->GetNbinsX();
  int nBinsData_L1L2dr = h1data_L1L2dr->GetNbinsX();
  int nBinsData_LKdr   = h1data_LKdr->GetNbinsX();
  int nBinsData_L1iso  = h1data_L1iso->GetNbinsX();
  int nBinsData_Kiso   = h1data_Kiso->GetNbinsX();
  int nBinsData_Passymetry = h1data_Passymetry->GetNbinsX();
  int nBinsData_Kip3d  = h1data_Kip3d->GetNbinsX();
  int nBinsData_KLmassD0 = h1data_KLmassD0->GetNbinsX(); 

  for (int ii=1; ii<=nBinsData_xgb; ii++)  { 
    if (h1data_xgb->GetBinContent(ii)<0) h1data_xgb->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_L1id; ii++) { 
    if (h1data_L1id->GetBinContent(ii)<0) h1data_L1id->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_L2id; ii++) { 
    if (h1data_L2id->GetBinContent(ii)<0) h1data_L2id->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Bprob; ii++) { 
    if (h1data_Bprob->GetBinContent(ii)<0) h1data_Bprob->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_BsLxy; ii++) { 
    if (h1data_BsLxy->GetBinContent(ii)<0) h1data_BsLxy->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Bcos; ii++) { 
    if (h1data_Bcos->GetBinContent(ii)<0) h1data_Bcos->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_L1pt; ii++) { 
    if (h1data_L1pt->GetBinContent(ii)<0) h1data_L1pt->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_L2pt; ii++) { 
    if (h1data_L2pt->GetBinContent(ii)<0) h1data_L2pt->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Bpt; ii++) { 
    if (h1data_Bpt->GetBinContent(ii)<0) h1data_Bpt->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Kpt; ii++) { 
    if (h1data_Kpt->GetBinContent(ii)<0) h1data_Kpt->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_LKdz; ii++) { 
    if (h1data_LKdz->GetBinContent(ii)<0) h1data_LKdz->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_L1L2dr; ii++) { 
    if (h1data_L1L2dr->GetBinContent(ii)<0) h1data_L1L2dr->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_LKdr; ii++) { 
    if (h1data_LKdr->GetBinContent(ii)<0) h1data_LKdr->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_L1iso; ii++) { 
    if (h1data_L1iso->GetBinContent(ii)<0) h1data_L1iso->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Kiso; ii++) { 
    if (h1data_Kiso->GetBinContent(ii)<0) h1data_Kiso->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Passymetry; ii++) { 
    if (h1data_Passymetry->GetBinContent(ii)<0) h1data_Passymetry->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_Kip3d; ii++) { 
    if (h1data_Kip3d->GetBinContent(ii)<0) h1data_Kip3d->SetBinContent(ii,0); }
  for (int ii=1; ii<=nBinsData_KLmassD0; ii++) { 
    if (h1data_KLmassD0->GetBinContent(ii)<0) h1data_KLmassD0->SetBinContent(ii,0); }

  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  // 
  drawTH1pair(h1data_xgb,   h1mc_xgb,    "BDT",    "a.u.",1.,"bdt",  "./",2, 0, "Data","MC");
  drawTH1pair(h1data_L1id,  h1mc_L1id,   "L1id",   "a.u.",1.,"L1id", "./",2, 0, "Data","MC");
  drawTH1pair(h1data_L2id,  h1mc_L2id,   "L2id",   "a.u.",1.,"L2id", "./",2, 0, "Data","MC");
  drawTH1pair(h1data_Bprob, h1mc_Bprob,  "Bprob",  "a.u.",1.,"Bprob","./",2, 0, "Data","MC");
  drawTH1pair(h1data_BsLxy, h1mc_BsLxy,  "BsLxy",  "a.u.",1.,"BsLxy","./",2, 0, "Data","MC");
  drawTH1pair(h1data_Bcos,  h1mc_Bcos,   "Bcos",   "a.u.",1.,"Bcos", "./",2, 0, "Data","MC");       
  drawTH1pair(h1data_L1pt,  h1mc_L1pt,   "L1pt",   "a.u.",1.,"L1pt", "./",2, 0, "Data","MC");       
  drawTH1pair(h1data_L2pt,  h1mc_L2pt,   "L2pt",   "a.u.",1.,"L2pt", "./",2, 0, "Data","MC");        
  drawTH1pair(h1data_Bpt,   h1mc_Bpt,    "Bpt",    "a.u.",1.,"Bpt", "./", 2, 0, "Data","MC");        
  drawTH1pair(h1data_Kpt,   h1mc_Kpt,    "Kpt",    "a.u.",1.,"Kpt",  "./",2, 0, "Data","MC");
  drawTH1pair(h1data_LKdz,  h1mc_LKdz,   "LKdz",   "a.u.",1.,"LKdz", "./",2, 0, "Data","MC");
  drawTH1pair(h1data_L1L2dr, h1mc_L1L2dr,  "L1L2dr",  "a.u.",1.,"L1L2dr","./",2, 0, "Data","MC");       
  drawTH1pair(h1data_LKdr,  h1mc_LKdr,   "LKdr",   "a.u.",1.,"LKdr", "./",2, 0, "Data","MC");
  drawTH1pair(h1data_L1iso, h1mc_L1iso,  "L1iso",  "a.u.",1.,"L1iso","./",2, 0, "Data","MC");       
  drawTH1pair(h1data_Kiso,  h1mc_Kiso,   "Kiso",   "a.u.",1.,"Kiso", "./",2, 0, "Data","MC");
  drawTH1pair(h1data_Passymetry, h1mc_Passymetry,  "Passymetry",  "a.u.",1.,"Passymetry","./",2, 0, "Data","MC");      
  drawTH1pair(h1data_Kip3d, h1mc_Kip3d,  "Kip3d",  "a.u.",1.,"Kip3d","./",2, 0, "Data","MC");        
  drawTH1pair(h1data_KLmassD0, h1mc_KLmassD0,  "KLmassD0",  "a.u.",1.,"KLmassD0","./",2, 0, "Data","MC");  
}
