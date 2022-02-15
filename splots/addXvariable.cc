#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>
#include "TLorentzVector.h"

using namespace std;

void addXvariable(const char* filename) {

  cout << "Formatting " << filename << endl;  

  // Original ntuple   
  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(TString("./")+TString(filename)+TString(".root"));
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("mytreefit");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }
  
  fileOrig->cd();
  if (!treeOrig) {
    cout << "tree not existing !" << endl; 
    return;    
  }

  treeOrig->SetMakeClass(0);
  cout << "TreeOrig->Size = "<< treeOrig->GetEntries() << endl;

  // number of entries saved in the first tree
  int nentriesOrig = treeOrig->GetEntries();   

  // Tree for the final format
  TFile *fileNew = TFile::Open(TString("./")+TString(filename)+TString("_withX.root"),"recreate");
  fileNew->ls();
  fileNew->cd();
  TTree *treeNew = new TTree("mytreefit","Tree with x variable for splots");
  treeNew->SetAutoSave(-99999999999);
  treeNew->SetAutoFlush(-99999999999);

  std::vector<TTree*> trees; 
  trees.push_back(treeNew);

  // original tree leaves
  Double_t        xgb = 0.;
  Double_t        Bmass = 0.;
  Double_t        Mll = 0.;
  Double_t        Npv = 0.;
  Double_t        Bprob = 0.;
  Double_t        BsLxy = 0.;
  Double_t        L1pt = 0.;
  Double_t        L2pt = 0.;
  Double_t        Bpt = 0.;
  Double_t        Kpt = 0.;
  Double_t        Bcos = 0.;
  Double_t        LKdz = 0.;
  Double_t        L1L2dr = 0.;
  Double_t        LKdr = 0.;
  Double_t        L2id = 0.;
  Double_t        Kiso = 0.;
  Double_t        BBDphi = 0.;
  Double_t        BTrkdxy2 = 0.;
  Double_t        L1id = 0.;
  Double_t        L1iso = 0.;
  Double_t        KLmassD0 = 0.;
  Double_t        Passymetry = 0.;
  Double_t        Kip3d = 0.;
  Double_t        Kip3dErr = 0.;

  // List of branches - original tree
  TBranch        *b_xgb;   //!
  TBranch        *b_Bmass;   //!
  TBranch        *b_Mll;   //!
  TBranch        *b_Npv;   //!
  TBranch        *b_Bprob;   //!
  TBranch        *b_BsLxy;   //!
  TBranch        *b_L1pt;   //!
  TBranch        *b_L2pt;   //!
  TBranch        *b_Bpt;   //!
  TBranch        *b_Kpt;   //!
  TBranch        *b_Bcos;   //!
  TBranch        *b_LKdz;   //!
  TBranch        *b_L1L2dr;   //!
  TBranch        *b_LKdr;   //!
  TBranch        *b_L2id;   //!
  TBranch        *b_Kiso;   //!
  TBranch        *b_BBDphi;   //!
  TBranch        *b_BTrkdxy2;   //!
  TBranch        *b_L1id;   //!
  TBranch        *b_L1iso;   //!
  TBranch        *b_KLmassD0;   //!
  TBranch        *b_Passymetry;   //!
  TBranch        *b_Kip3d;   //!
  TBranch        *b_Kip3dErr;   //!

  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("xgb", &xgb, &b_xgb);
  treeOrig->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
  treeOrig->SetBranchAddress("Mll", &Mll, &b_Mll);
  treeOrig->SetBranchAddress("Npv", &Npv, &b_Npv);
  treeOrig->SetBranchAddress("Bprob", &Bprob, &b_Bprob);
  treeOrig->SetBranchAddress("BsLxy", &BsLxy, &b_BsLxy);
  treeOrig->SetBranchAddress("L1pt", &L1pt, &b_L1pt);
  treeOrig->SetBranchAddress("L2pt", &L2pt, &b_L2pt);
  treeOrig->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
  treeOrig->SetBranchAddress("Kpt", &Kpt, &b_Kpt);
  treeOrig->SetBranchAddress("Bcos", &Bcos, &b_Bcos);
  treeOrig->SetBranchAddress("LKdz", &LKdz, &b_LKdz);
  treeOrig->SetBranchAddress("L1L2dr", &L1L2dr, &b_L1L2dr);
  treeOrig->SetBranchAddress("LKdr", &LKdr, &b_LKdr);
  treeOrig->SetBranchAddress("L2id", &L2id, &b_L2id);
  treeOrig->SetBranchAddress("Kiso", &Kiso, &b_Kiso);
  treeOrig->SetBranchAddress("BBDphi", &BBDphi, &b_BBDphi);
  treeOrig->SetBranchAddress("BTrkdxy2", &BTrkdxy2, &b_BTrkdxy2);
  treeOrig->SetBranchAddress("L1id", &L1id, &b_L1id);
  treeOrig->SetBranchAddress("L1iso", &L1iso, &b_L1iso);
  treeOrig->SetBranchAddress("KLmassD0", &KLmassD0, &b_KLmassD0);
  treeOrig->SetBranchAddress("Passymetry", &Passymetry, &b_Passymetry);
  treeOrig->SetBranchAddress("Kip3d", &Kip3d, &b_Kip3d);
  treeOrig->SetBranchAddress("Kip3dErr", &Kip3dErr, &b_Kip3dErr);

  // New variable
  Double_t x;

  // Branches for the new tree
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];

    theTreeNew->Branch("xgb",    &xgb,    "xgb/D");
    theTreeNew->Branch("Bmass",  &Bmass,  "Bmass/D");
    theTreeNew->Branch("Mll",    &Mll,    "Mll/D");
    theTreeNew->Branch("Npv",    &Npv,    "Npv/D");
    theTreeNew->Branch("Bprob",  &Bprob,  "Bprob/D");
    theTreeNew->Branch("BsLxy",  &BsLxy,  "BsLxy/D");
    theTreeNew->Branch("L1pt",   &L1pt,   "L1pt/D");
    theTreeNew->Branch("L2pt",   &L2pt,   "L2pt/D");
    theTreeNew->Branch("Bpt",    &Bpt,    "Bpt/D");
    theTreeNew->Branch("Kpt",    &Kpt,    "Kpt/D");
    theTreeNew->Branch("Bcos",   &Bcos,   "Bcos/D");
    theTreeNew->Branch("LKdz",   &LKdz,   "LKdz/D");
    theTreeNew->Branch("L1L2dr", &L1L2dr, "L1L2dr/D");
    theTreeNew->Branch("LKdr",   &LKdr,   "LKdr/D");
    theTreeNew->Branch("L2id",   &L2id,   "L2id/D");
    theTreeNew->Branch("Kiso",   &Kiso,   "Kiso/D");
    theTreeNew->Branch("BBDphi", &BBDphi, "BBDphi/D");
    theTreeNew->Branch("BTrkdxy2", &BTrkdxy2, "BTrkdxy2/D");
    theTreeNew->Branch("L1id",  &L1id,    "L1id/D");
    theTreeNew->Branch("L1iso", &L1iso,   "L1iso/D");
    theTreeNew->Branch("KLmassD0",   &KLmassD0, "KLmassD0/D");
    theTreeNew->Branch("Passymetry", &Passymetry, "Passymetry/D");
    theTreeNew->Branch("Kip3d",    &Kip3d, "Kip3d/D");
    theTreeNew->Branch("Kip3dErr", &Kip3dErr, "Kip3dErr/D");
    theTreeNew->Branch("x", &x, "x/D");
  }
  
  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    x = Bmass;     

    treeNew->Fill();
  }
  
  // new format
  treeNew->Write();
  fileNew->Close();
  fileNew->ls();
  
  fileOrig->cd();
  fileOrig->Close();  
}

