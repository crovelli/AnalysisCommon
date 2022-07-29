#define correlationsBdtAntiD0_cxx
#include "correlationsBdtAntiD0.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLegend.h"
#include "TPaveText.h"
#include <iostream> 

using namespace std;

void correlationsBdtAntiD0::Loop()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // chiara
  bool isSonly = 1;
  bool isBonly = 0;
  bool isGenericB = 0;
  bool bdtCut0 = 0;
  float nominalBDT = 8.3;
  // chiara

  TH2F *corrS;
  TH2F *corrB;
  if (bdtCut0==0) {
    corrS = new TH2F("corrS", "corrS", 100, -15, 12, 50, 0., 6.);
    corrB = new TH2F("corrB", "corrB", 100, -15, 12, 50, 0., 6.); 
  }
  else {
    corrS = new TH2F("corrS", "corrS", 100, 0, 12, 50, 0., 6.);
    corrB = new TH2F("corrB", "corrB", 100, 0., 12, 50, 0., 6.); 
  }
  TH1F *antiD0passingBdtS = new TH1F("antiD0passingBdtS", "antiD0passingBdtS", 50, 0., 6.);  
  TH1F *antiD0passingBdtB = new TH1F("antiD0passingBdtB", "antiD0passingBdtB", 50, 0., 6.);  
  TH1F *antiD0failingBdtS = new TH1F("antiD0failingBdtS", "antiD0failingBdtS", 50, 0., 6.);  
  TH1F *antiD0failingBdtB = new TH1F("antiD0failingBdtB", "antiD0failingBdtB", 50, 0., 6.);  
  TH1F *antiD0BdtLt0S = new TH1F("antiD0Lt0BdtS", "antiD0Lt0BdtS", 50, 0., 6.);  
  TH1F *antiD0BdtLt0B = new TH1F("antiD0Lt0BdtB", "antiD0Lt0BdtB", 50, 0., 6.);  
  //
  antiD0passingBdtS->GetXaxis()->SetTitle("antiD0");
  antiD0passingBdtB->GetXaxis()->SetTitle("antiD0");
  antiD0failingBdtS->GetXaxis()->SetTitle("antiD0");
  antiD0failingBdtB->GetXaxis()->SetTitle("antiD0");
  antiD0BdtLt0B->GetXaxis()->SetTitle("antiD0");
  antiD0BdtLt0S->GetXaxis()->SetTitle("antiD0");
  //
  antiD0passingBdtS->SetTitle("");
  antiD0passingBdtB->SetTitle("");
  antiD0failingBdtS->SetTitle("");
  antiD0failingBdtB->SetTitle("");
  antiD0BdtLt0B ->SetTitle("");
  antiD0BdtLt0S ->SetTitle("");
  //
  corrS -> GetXaxis()->SetTitle("BDT output");
  corrS -> GetYaxis()->SetTitle("antiD0");
  corrS -> SetTitle("");
  corrB -> GetXaxis()->SetTitle("BDT output");
  corrB -> GetYaxis()->SetTitle("antiD0");
  corrB -> SetTitle("");
    

  // Loop over events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
  
    // Basic selection
    if ( Bmass>5.7 || Bmass<4.7 ) continue;
    if ( Mll>3.2 || Mll<2.9 )  continue;
    if ( bdtCut0==1 && xgb<0 ) continue;

    if (isSonly) { 
      if (xgb<=0) antiD0BdtLt0S -> Fill(KLmassD0);
      if (xgb>0 && xgb<=nominalBDT) antiD0failingBdtS -> Fill(KLmassD0);
      if (xgb>nominalBDT) antiD0passingBdtS -> Fill(KLmassD0);
      corrS -> Fill(xgb, KLmassD0);
    }

    if (isBonly) {
      if (xgb<=0) antiD0BdtLt0B -> Fill(KLmassD0);
      if (xgb>0 && xgb<=nominalBDT) antiD0failingBdtB -> Fill(KLmassD0);
      if (xgb>nominalBDT) antiD0passingBdtB -> Fill(KLmassD0);
      corrB -> Fill(xgb, KLmassD0);
    }

    if (isGenericB) {
      //if (isOtherB==1) {
      if (xgb<=0) antiD0BdtLt0B -> Fill(KLmassD0);
      if (xgb>0 && xgb<=nominalBDT) antiD0failingBdtB -> Fill(KLmassD0);
      if (xgb>nominalBDT) antiD0passingBdtB -> Fill(KLmassD0);
      corrB -> Fill(xgb, KLmassD0);
    }//}

  } // Loop over events


  // cosmetics
  antiD0passingBdtS -> SetLineWidth(2);
  antiD0failingBdtS -> SetLineWidth(2);
  antiD0BdtLt0S     -> SetLineWidth(2);
  antiD0passingBdtS -> SetLineColor(2);
  antiD0failingBdtS -> SetLineColor(3);
  antiD0BdtLt0S     -> SetLineColor(1);
  //
  antiD0passingBdtB -> SetLineWidth(2);
  antiD0failingBdtB -> SetLineWidth(2);
  antiD0BdtLt0B     -> SetLineWidth(2);
  antiD0passingBdtB -> SetLineColor(2);
  antiD0failingBdtB -> SetLineColor(3);
  antiD0BdtLt0B     -> SetLineColor(1);


  // Legend
  TLegend *leg;
  leg = new TLegend(0.15,0.70,0.50,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(antiD0BdtLt0S, "BDT<0", "l");
  leg->AddEntry(antiD0failingBdtS, "0<BDT<8.3", "l");
  leg->AddEntry(antiD0passingBdtS, "BDT>8.3", "l");

  // Plots
  gStyle->SetOptStat(0);

  if (isSonly) {   
    TCanvas c1("c1","",1);
    antiD0BdtLt0S->DrawNormalized();
    antiD0passingBdtS->DrawNormalized("same");
    antiD0failingBdtS->DrawNormalized("same");
    leg->Draw();
    c1.SaveAs("S_antiD0andBdt.png");
  }

  if (isBonly || isGenericB) { 
    TCanvas c2("c2","",1);
    antiD0passingBdtB->DrawNormalized();
    antiD0failingBdtB->DrawNormalized("same");
    antiD0BdtLt0B->DrawNormalized("same");
    leg->Draw();
    c2.SaveAs("B_antiD0andBdt.png");
  }

  if (isSonly) { 
    TCanvas c3("c3","",1);
    corrS->Draw("colz");
    cout << "CorrS: corr = " << corrS->GetCorrelationFactor() << endl;
    c3.SaveAs("corrS.png");
  }

  if (isBonly || isGenericB) { 
    TCanvas c4("c4","",1);
    corrB->Draw("colz");
    cout << "CorrB: corr = " << corrB->GetCorrelationFactor() << endl;
    c4.SaveAs("corrB.png");
  }


}
