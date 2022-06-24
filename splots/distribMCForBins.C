#define distribMCForBins_cxx
#include "distribMCForBins.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>   
#include <iostream>  

// MC distributions to compare the BDT output in different q2 bins

void distribMCForBins::Loop(int q2bin)
{
  // Define q^2 ranges
  float q2inf = 0;
  float q2sup = 1000;
  if (q2bin==0) { q2inf=1.05; q2sup=2.45; }
  if (q2bin==1) { q2inf=2.90; q2sup=3.20; }
  if (q2bin==2) { q2inf=3.55; q2sup=3.80; }
  if (q2bin==3) { q2inf=4.00; q2sup=4.80; }

  if (fChain == 0) return;

  // Prepare histos
  TH1F *h1mc_xgb      = new TH1F("h1mc_xgb",      "h1mc_xgb",      30,  0.0, 15.);
  TH1F *h1mc_mll      = new TH1F("h1mc_mll",      "h1mc_mll",      60,  1.0,  5.);
  TH1F *h1mc_Bpt      = new TH1F("h1mc_Bpt",      "h1mc_Bpt",      80,  0.0, 80.); 
  TH1F *h1mc_KLmassD0 = new TH1F("h1mc_KLmassD0", "h1mc_KLmassD0", 30,  0.0,  6.); 
  TH1F *h1mc_Ele1pt   = new TH1F("h1mc_Ele1pt",   "h1mc_Ele1pt",   80,  0.0, 40.); 
  TH1F *h1mc_Ele2pt   = new TH1F("h1mc_Ele2pt",   "h1mc_Ele2pt",   40,  0.0, 20.); 
  TH1F *h1mc_ElesDr   = new TH1F("h1mc_ElesDr",   "h1mc_ElesDr",   20,  0.0, 2.); 

  // Sumw2
  h1mc_xgb -> Sumw2();
  h1mc_mll -> Sumw2();
  h1mc_Bpt -> Sumw2();
  h1mc_KLmassD0 -> Sumw2();
  h1mc_Ele1pt -> Sumw2();
  h1mc_Ele2pt -> Sumw2();
  h1mc_ElesDr -> Sumw2();

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

    //if (KLmassD0<2.0) continue;
    //if (xgb<8)        continue;


    // Fill all histos
    h1mc_mll -> Fill(Mll);
    h1mc_xgb -> Fill(xgb);
    h1mc_Bpt -> Fill(Bpt);
    h1mc_KLmassD0 -> Fill(KLmassD0);
    h1mc_Ele1pt -> Fill(L1pt);
    h1mc_Ele2pt -> Fill(L2pt);
    h1mc_ElesDr -> Fill(L1L2dr);

  } // Loop over entries


  // Save outputs
  TFile myfile("myFileMC.root","RECREATE");
  myfile.cd();
  h1mc_mll -> Write();
  h1mc_xgb -> Write();
  h1mc_Bpt -> Write();
  h1mc_KLmassD0 -> Write();
  h1mc_Ele1pt -> Write();       
  h1mc_Ele2pt -> Write();       
  h1mc_ElesDr -> Write();       
  myfile.Close();
}
