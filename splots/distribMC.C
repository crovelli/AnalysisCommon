#define distribMC_cxx
#include "distribMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   
#include <iostream>  

// To be run on MC tnp formatted ntuples to get distributions for 
// comparison with sPlots 
// Select signal and background based on match with MC-truth 

void distribMC::Loop(int q2bin, float bdtCut, int isPFPF)
{
  // Define q^2 ranges
  float q2inf = 0;
  float q2sup = 1000;
  if (q2bin==0) { q2inf=1.05; q2sup=2.45; }
  if (q2bin==1) { q2inf=2.90; q2sup=3.20; }
  if (q2bin==2) { q2inf=3.55; q2sup=3.80; }
  if (q2bin==3) { q2inf=4.00; q2sup=4.80; }

  if (fChain == 0) return;
  
  float theinf = bdtCut-2;
  float thedelta = 15.-theinf;
  int thebin = thedelta/0.5;
  if (q2bin==0) {
    theinf = bdtCut-2;
    thedelta = 12.5-theinf;
    thebin = thedelta/0.5;
  }

  TH1F *h1mc_mll    = new TH1F("h1mc_mll","h1mc_mll",60,  1., 5.);
  TH1F *h1mc_bmass  = new TH1F("h1mc_bmass","h1mc_bmass",60,  4.5, 5.9);

  TH1F *h1mc_xgb;
  //if (q2bin!=0) h1mc_xgb = new TH1F("h1mc_xgb","h1mc_xgb",thebin,theinf,15);
  //if (q2bin==0) h1mc_xgb = new TH1F("h1mc_xgb","h1mc_xgb",thebin,theinf,12.5);
  h1mc_xgb = new TH1F("h1mc_xgb","h1mc_xgb",24,3.0,15.);

  TH1F *h1mc_L1id;
  TH1F *h1mc_L2id;
  if (isPFPF==1) {
    if (q2bin!=0) { 
      h1mc_L1id   = new TH1F("h1mc_L1id","h1mc_L1id",32,-4.,7.);
      h1mc_L2id   = new TH1F("h1mc_L2id","h1mc_L2id",32,-4.,7.);
    } 
    if (q2bin==0) { 
      h1mc_L1id   = new TH1F("h1mc_L1id","h1mc_L1id",22,-4.,7.);
      h1mc_L2id   = new TH1F("h1mc_L2id","h1mc_L2id",22,-4.,7.);
    } 
  } else {
    if (q2bin!=0) { 
      h1mc_L1id   = new TH1F("h1mc_L1id","h1mc_L1id",40,-4.,16.);
      h1mc_L2id   = new TH1F("h1mc_L2id","h1mc_L2id",40,-4.,16.);
    } 
    if (q2bin==0) { 
      h1mc_L1id   = new TH1F("h1mc_L1id","h1mc_L1id",20,-4.,16.);
      h1mc_L2id   = new TH1F("h1mc_L2id","h1mc_L2id",20,-4.,16.);
    } 
  }
  TH1F *h1mc_Bprob  = new TH1F("h1mc_Bprob","h1mc_Bprob",10,0.,1.);
  TH1F *h1mc_BsLxy  = new TH1F("h1mc_BsLxy","h1mc_BsLxy",10,0.,100.);
  TH1F *h1mc_Bcos   = new TH1F("h1mc_Bcos","h1mc_Bcos",10,0.99,1.);
  TH1F *h1mc_L1pt   = new TH1F("h1mc_L1pt","h1mc_L1pt",60,0.,30.);
  TH1F *h1mc_L2pt   = new TH1F("h1mc_L2pt","h1mc_L2pt",40,0.,20.);
  TH1F *h1mc_Bpt    = new TH1F("h1mc_Bpt","h1mc_Bpt",80,0.,80.);
  TH1F *h1mc_Kpt    = new TH1F("h1mc_Kpt","h1mc_Kpt",40,0.,20.);
  TH1F *h1mc_LKdz   = new TH1F("h1mc_LKdz","h1mc_LKdz",20,0.,1.);
  TH1F *h1mc_L1L2dr = new TH1F("h1mc_L1L2dr","h1mc_L1L2dr",20,0.,2.);
  TH1F *h1mc_LKdr   = new TH1F("h1mc_LKdr","h1mc_LKdr",40,0.,3.);
  TH1F *h1mc_L1iso  = new TH1F("h1mc_L1iso","h1mc_L1iso",30,0.,30.);
  TH1F *h1mc_Kiso   = new TH1F("h1mc_Kiso","h1mc_Kiso",30,0.,30.);
  TH1F *h1mc_Kip3d  = new TH1F("h1mc_Kip3d","h1mc_Kip3d",14,-0.07,0.07);
  TH1F *h1mc_Passymetry = new TH1F("h1mc_Passymetry","h1mc_Passymetry",20,-1.,1.);
  TH1F *h1mc_KLmassD0   = new TH1F("h1mc_KLmassD0",  "h1mc_KLmassD0",  30, 0.,6.);
  h1mc_mll    -> Sumw2(); 
  h1mc_bmass  -> Sumw2(); 
  h1mc_xgb    -> Sumw2(); 
  h1mc_L1id   -> Sumw2();
  h1mc_L2id   -> Sumw2();
  h1mc_Bprob  -> Sumw2();
  h1mc_BsLxy  -> Sumw2();
  h1mc_Bcos   -> Sumw2();
  h1mc_L1pt   -> Sumw2();
  h1mc_L2pt   -> Sumw2();
  h1mc_Bpt    -> Sumw2();
  h1mc_Kpt    -> Sumw2();
  h1mc_LKdz   -> Sumw2();
  h1mc_L1L2dr -> Sumw2();
  h1mc_LKdr   -> Sumw2();
  h1mc_L1iso  -> Sumw2();
  h1mc_Kiso   -> Sumw2();
  h1mc_Kip3d  -> Sumw2();
  h1mc_Passymetry -> Sumw2();
  h1mc_KLmassD0   -> Sumw2();

  // Loop over entries
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  cout << "nentries = " << nentries << endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;

    if (jentry%1000==0) cout << jentry << endl;

    // Same selection as applied in data
    if (Bmass>5.7 || Bmass<4.7) continue;
    if (Mll>q2sup || Mll<q2inf) continue;
    if (xgb<bdtCut)   continue;
    if (KLmassD0<2.0) continue;      

    // HLT and MC-match should be applied already at ntuple level

    // Fill all histos
    h1mc_mll         -> Fill(Mll);
    h1mc_bmass       -> Fill(Bmass);
    h1mc_xgb         -> Fill(xgb);
    h1mc_L1id        -> Fill(L1id);
    h1mc_L2id        -> Fill(L2id);   
    h1mc_Bprob       -> Fill(Bprob);  
    h1mc_BsLxy       -> Fill(BsLxy);    
    h1mc_Bcos        -> Fill(Bcos);     
    h1mc_L1pt        -> Fill(L1pt);     
    h1mc_L2pt        -> Fill(L2pt);  
    h1mc_Bpt         -> Fill(Bpt);  
    h1mc_Kpt         -> Fill(Kpt);      
    h1mc_LKdz        -> Fill(LKdz);  
    h1mc_L1L2dr      -> Fill(L1L2dr);   
    h1mc_LKdr        -> Fill(LKdr);  
    h1mc_L1iso       -> Fill(L1iso);  
    h1mc_Kiso        -> Fill(Kiso);  
    h1mc_Passymetry  -> Fill(Passymetry);  
    h1mc_Kip3d       -> Fill(Kip3d);  
    h1mc_KLmassD0    -> Fill(KLmassD0);

  } // Loop over entries

  // Save outputs
  TFile myfile("myFileMC.root","RECREATE");
  myfile.cd();
  h1mc_mll         -> Write();
  h1mc_bmass       -> Write();
  h1mc_xgb         -> Write();
  h1mc_L1id        -> Write();  
  h1mc_L2id        -> Write();  
  h1mc_Bprob       -> Write();  
  h1mc_BsLxy       -> Write();  
  h1mc_Bcos        -> Write();  
  h1mc_L1pt        -> Write();  
  h1mc_L2pt        -> Write();  
  h1mc_Bpt         -> Write();  
  h1mc_Kpt         -> Write();  
  h1mc_LKdz        -> Write();  
  h1mc_L1L2dr      -> Write();  
  h1mc_LKdr        -> Write();  
  h1mc_L1iso       -> Write();  
  h1mc_Kiso        -> Write();  
  h1mc_Passymetry  -> Write();  
  h1mc_Kip3d       -> Write();  
  h1mc_KLmassD0    -> Write();
  myfile.Close();
}
