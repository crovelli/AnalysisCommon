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
  frame->GetYaxis()->SetRangeUser(-1.,4.);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
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

  delete canvas;
  frame->Reset("ICES");

  TFile fileOut(outputFILE.c_str(), "UPDATE");
  fileOut.cd();
  ratio->Write( ("ratio_" + canvasName).c_str());
  h2->Write( ("MC_" + canvasName).c_str());
  h1->Write( ("Data_" + canvasName).c_str());
  fileOut.Close();
}

void mcVsSplots(int q2bin)
{
  // Input files 
  TFile *fileMC     = new TFile("./myFileMC.root");
  TFile *fileSPlots = new TFile("./myFileSPlots.root");
  
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

  // all bins
  h1mc_L1pt->Rebin(4);
  h1mc_L1id->Rebin(2);
  h1mc_L2id->Rebin(2);
  h1mc_L2pt->Rebin();
  h1mc_LKdr->Rebin();
  h1mc_Bpt->Rebin(4);
  h1mc_Kpt->Rebin(4);
  h1mc_L1iso->Rebin();
  h1mc_Kiso->Rebin();
  h1data_L1pt->Rebin(4);
  h1data_L1id->Rebin(2);
  h1data_L2id->Rebin(2);
  h1data_L2pt->Rebin();
  h1data_LKdr->Rebin();
  h1data_Bpt->Rebin(4);
  h1data_Kpt->Rebin(4);
  h1data_L1iso->Rebin();
  h1data_Kiso->Rebin();
  h1mc_KLmassD0->Rebin();
  h1data_KLmassD0->Rebin();

  if (q2bin==2) {
    h1mc_L1L2dr->Rebin();
    h1data_L1L2dr->Rebin();
    h1mc_L2id->Rebin();
    h1data_L2id->Rebin();
    h1mc_L2pt->Rebin();
    h1data_L2pt->Rebin();
    h1mc_LKdr->Rebin();
    h1data_LKdr->Rebin();
    h1mc_xgb->Rebin();
    h1data_xgb->Rebin();
  }

  /*
  if (q2bin==0 || q2bin==2) {
    h1data_xgb->Rebin();
    h1mc_xgb->Rebin();
    h1data_Kpt->Rebin();
    h1mc_Kpt->Rebin();
    h1data_LKdr->Rebin();
    h1mc_LKdr->Rebin();
  }
  */


  // Plots
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  // 
  cout << h1data_xgb->GetMinimum() << " " << h1mc_xgb->GetMinimum() << endl;
  drawTH1pair(h1data_xgb,   h1mc_xgb,    "BDT",    "a.u.",1.,"bdt",  "./",2, 0, "Data","MC");
  /*
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
  */
  drawTH1pair(h1data_KLmassD0, h1mc_KLmassD0,  "KLmassD0",  "a.u.",1.,"KLmassD0","./",2, 0, "Data","MC");  
}
