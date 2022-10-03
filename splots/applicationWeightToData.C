#define applicationWeightToData_cxx
#include "applicationWeightToData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   

void applicationWeightToData::Loop(int theVariable)       // theVariable must be coherent with .h : 0=ele1pt, 1=ele2pt, 2=kpt, 3=elesDr
{
  TH1F *jPsiEle1pt_fromLoop_now = new TH1F("jPsiEle1pt_fromLoop_now","jPsiEle1pt_fromLoop_now", 80,  0.0, 40.);
  TH1F *jPsiEle1pt_fromLoop_ww  = new TH1F("jPsiEle1pt_fromLoop_ww", "jPsiEle1pt_fromLoop_ww",  80,  0.0, 40.);
  TH1F *jPsiEle2pt_fromLoop_now = new TH1F("jPsiEle2pt_fromLoop_now","jPsiEle2pt_fromLoop_now", 40,  0.0, 20.);
  TH1F *jPsiEle2pt_fromLoop_ww  = new TH1F("jPsiEle2pt_fromLoop_ww", "jPsiEle2pt_fromLoop_ww",  40,  0.0, 20.);
  TH1F *jPsiKpt_fromLoop_now    = new TH1F("jPsiKpt_fromLoop_now","jPsiKpt_fromLoop_now", 40,  0.0, 20.);
  TH1F *jPsiKpt_fromLoop_ww     = new TH1F("jPsiKpt_fromLoop_ww", "jPsiKpt_fromLoop_ww",  40,  0.0, 20.);
  TH1F *jPsiElesDr_fromLoop_now = new TH1F("jPsiElesDr_fromLoop_now","jPsiElesDr_fromLoop_now", 20,  0.0, 2.);
  TH1F *jPsiElesDr_fromLoop_ww  = new TH1F("jPsiElesDr_fromLoop_ww", "jPsiElesDr_fromLoop_ww",  20,  0.0, 2.);
  TH1F *jPsiBmass_fromLoop_now  = new TH1F("jPsiBmass_fromLoop_now", "jPsiBmass_fromLoop_now",  50,  4.7, 5.7);
  TH1F *jPsiBmass_fromLoop_ww   = new TH1F("jPsiBmass_fromLoop_ww",  "jPsiBmass_fromLoop_ww",   50,  4.7, 5.7);
  TH1F *jPsiXgb_fromLoop_now    = new TH1F("jPsiXgb_fromLoop_now", "jPsiXgb_fromLoop_now",  45,  0.0, 15.);
  TH1F *jPsiXgb_fromLoop_ww     = new TH1F("jPsiXgb_fromLoop_ww",  "jPsiXgb_fromLoop_ww",   45,  0.0, 15.);
  //
  jPsiEle1pt_fromLoop_now -> Sumw2();
  jPsiEle1pt_fromLoop_ww  -> Sumw2();
  jPsiEle2pt_fromLoop_now -> Sumw2();
  jPsiEle2pt_fromLoop_ww  -> Sumw2();
  jPsiKpt_fromLoop_now    -> Sumw2();
  jPsiKpt_fromLoop_ww     -> Sumw2();
  jPsiElesDr_fromLoop_now -> Sumw2();
  jPsiElesDr_fromLoop_ww  -> Sumw2();
  jPsiBmass_fromLoop_now  -> Sumw2();
  jPsiBmass_fromLoop_ww   -> Sumw2();
  jPsiXgb_fromLoop_now    -> Sumw2();
  jPsiXgb_fromLoop_ww     -> Sumw2();
  //
  jPsiEle1pt_fromLoop_now -> Rebin(4);
  jPsiEle1pt_fromLoop_ww  -> Rebin(4);
  jPsiEle2pt_fromLoop_now -> Rebin(2);
  jPsiEle2pt_fromLoop_ww  -> Rebin(2);
  //
  //
  // cosmetics
  jPsiEle1pt_fromLoop_now -> SetTitle("");
  jPsiEle1pt_fromLoop_ww  -> SetTitle("");
  jPsiEle2pt_fromLoop_now -> SetTitle("");
  jPsiEle2pt_fromLoop_ww  -> SetTitle("");
  jPsiKpt_fromLoop_now    -> SetTitle("");
  jPsiKpt_fromLoop_ww     -> SetTitle("");
  jPsiElesDr_fromLoop_now -> SetTitle("");
  jPsiElesDr_fromLoop_ww  -> SetTitle("");
  jPsiBmass_fromLoop_now  -> SetTitle("");
  jPsiBmass_fromLoop_ww   -> SetTitle("");
  jPsiXgb_fromLoop_now    -> SetTitle("");
  jPsiXgb_fromLoop_ww     -> SetTitle("");
  jPsiEle1pt_fromLoop_now -> SetLineColor(2);   
  jPsiEle1pt_fromLoop_ww  -> SetLineColor(1);   
  jPsiEle2pt_fromLoop_now -> SetLineColor(2);   
  jPsiEle2pt_fromLoop_ww  -> SetLineColor(1);   
  jPsiKpt_fromLoop_now    -> SetLineColor(2);   
  jPsiKpt_fromLoop_ww     -> SetLineColor(1);   
  jPsiElesDr_fromLoop_now -> SetLineColor(2);
  jPsiElesDr_fromLoop_ww  -> SetLineColor(1);
  jPsiBmass_fromLoop_now  -> SetLineColor(2);
  jPsiBmass_fromLoop_ww   -> SetLineColor(1);
  jPsiXgb_fromLoop_now    -> SetLineColor(2);
  jPsiXgb_fromLoop_ww     -> SetLineColor(1);
  jPsiEle1pt_fromLoop_now -> SetLineWidth(2);   
  jPsiEle1pt_fromLoop_ww  -> SetLineWidth(2);   
  jPsiEle2pt_fromLoop_now -> SetLineWidth(2);   
  jPsiEle2pt_fromLoop_ww  -> SetLineWidth(2);   
  jPsiKpt_fromLoop_now    -> SetLineWidth(2);   
  jPsiKpt_fromLoop_ww     -> SetLineWidth(2);   
  jPsiElesDr_fromLoop_now -> SetLineWidth(2);
  jPsiElesDr_fromLoop_ww  -> SetLineWidth(2);
  jPsiBmass_fromLoop_now  -> SetLineWidth(2);
  jPsiBmass_fromLoop_ww   -> SetLineWidth(2);
  jPsiXgb_fromLoop_now    -> SetLineWidth(2);
  jPsiXgb_fromLoop_ww     -> SetLineWidth(2);
  

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Selection (denominator)
    if (Mll>3.2 || Mll<2.9) continue;
    if (Bmass>5.7 || Bmass<4.7) continue;
    if (xgb<0)        continue;
    if (KLmassD0<2.0) continue;      

    float weight;
    if (theVariable==0) weight = weightElePt1;
    if (theVariable==1) weight = weightElePt2;
    if (theVariable==2) weight = weightKPt;
    if (theVariable==3) weight = weightElesDr;

    jPsiEle1pt_fromLoop_now->Fill(L1pt);
    jPsiEle1pt_fromLoop_ww ->Fill(L1pt, weight);       
    jPsiEle2pt_fromLoop_now->Fill(L2pt);
    jPsiEle2pt_fromLoop_ww ->Fill(L2pt, weight);       
    jPsiKpt_fromLoop_now   ->Fill(Kpt);
    jPsiKpt_fromLoop_ww    ->Fill(Kpt, weight);       
    jPsiElesDr_fromLoop_now->Fill(L1L2dr);
    jPsiElesDr_fromLoop_ww ->Fill(L1L2dr, weight);       
    jPsiBmass_fromLoop_now ->Fill(Bmass);
    jPsiBmass_fromLoop_ww  ->Fill(Bmass, weight);       
    jPsiXgb_fromLoop_now   ->Fill(xgb);
    jPsiXgb_fromLoop_ww    ->Fill(xgb, weight);       
  }


  // plots
  gStyle->SetOptStat(0);

  TLegend leg2 (0.55,0.55,0.90,0.85);
  leg2.SetFillColor(0);
  leg2.SetFillStyle(0);
  leg2.SetBorderSize(0);
  leg2.AddEntry(jPsiElesDr_fromLoop_now, "Data, JPsi bin","PLE");
  leg2.AddEntry(jPsiElesDr_fromLoop_ww,  "Data, JPsi bin. Weighted","L");
 
  if (theVariable==0) {
    TCanvas c2a("c2a","With and wo weights",1);
    jPsiEle1pt_fromLoop_ww  -> DrawNormalized("hist");
    jPsiEle1pt_fromLoop_now -> DrawNormalized("samehist");
    leg2.Draw();
    c2a.SaveAs("Ele1pt_withAndWoWeight.png");
  }

  if (theVariable==1) {
    TCanvas c2b("c2b","With and wo weights",1);
    jPsiEle2pt_fromLoop_ww  -> DrawNormalized("hist");
    jPsiEle2pt_fromLoop_now -> DrawNormalized("samehist");
    leg2.Draw();
    c2b.SaveAs("Ele2pt_withAndWoWeight.png");
  }

  if (theVariable==2) {
    TCanvas c2c("c2c","With and wo weights",1);
    jPsiKpt_fromLoop_ww  -> DrawNormalized("hist");
    jPsiKpt_fromLoop_now -> DrawNormalized("samehist");
    leg2.Draw();
    c2c.SaveAs("Kpt_withAndWoWeight.png");
  }

  if (theVariable==3) {
    TCanvas c2d("c2d","With and wo weights",1);
    jPsiElesDr_fromLoop_ww  -> DrawNormalized("hist");
    jPsiElesDr_fromLoop_now -> DrawNormalized("samehist");
    leg2.Draw();
    c2d.SaveAs("ElesDr_withAndWoWeight.png");
  }

  TCanvas c2e("c2e","With and wo weights",1);
  jPsiBmass_fromLoop_ww  -> DrawNormalized("hist");
  jPsiBmass_fromLoop_now -> DrawNormalized("samehist");
  leg2.Draw();
  if (theVariable==0) c2e.SaveAs("Bmass_withAndWoEle1PtWeight.png");
  if (theVariable==1) c2e.SaveAs("Bmass_withAndWoEle2PtWeight.png");
  if (theVariable==2) c2e.SaveAs("Bmass_withAndWoKPtWeight.png");
  if (theVariable==3) c2e.SaveAs("Bmass_withAndWoElesDrWeight.png");

  TCanvas c2f("c2f","With and wo weights",1);
  jPsiXgb_fromLoop_ww  -> DrawNormalized("hist");
  jPsiXgb_fromLoop_now -> DrawNormalized("samehist");
  leg2.Draw();
  if (theVariable==0) c2f.SaveAs("Xgb_withAndWoEle1PtWeight.png");
  if (theVariable==1) c2f.SaveAs("Xgb_withAndWoEle2PtWeight.png");
  if (theVariable==2) c2f.SaveAs("Xgb_withAndWoKPtWeight.png");
  if (theVariable==3) c2f.SaveAs("Xgb_withAndWoElesDrWeight.png");
}
