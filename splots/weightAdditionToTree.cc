#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <iostream>
#include <vector>
#include <TRandom.h>
#include "TLorentzVector.h"

using namespace std;

double GetEle1ptWeight(double ele1pt);
double GetEle2ptWeight(double ele2pt);
double GetElesDrWeight(double elesdr);

// For weights application
std::vector<Double_t> anweightEle1pt_;
std::vector<Double_t> anweightEle2pt_;
std::vector<Double_t> anweightElesDr_;
std::vector<Double_t> lowedgeEle1pt_;
std::vector<Double_t> lowedgeEle2pt_;
std::vector<Double_t> lowedgeElesDr_;

void weightAdditionToTree(const char* dataFilename, const char* ratioFilename, int theVariable) {

  cout << "Formatting " << dataFilename << endl;  
  if (theVariable==0) cout << "adding weight corresponding to ele1 pT" << endl;
  if (theVariable==1) cout << "adding weight corresponding to ele2 pT" << endl;
  if (theVariable==2) cout << "adding weight corresponding to dR(ele1, ele2)" << endl;
  

  // Loading JPsi/LowQ2 ratios to compute weights
  TFile *ratioFileMCBins = 0;
  ratioFileMCBins = TFile::Open(TString("./")+TString(ratioFilename));
  TH1F *ratioEle1pt_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_Ele1Pt-JPsiToLowQ2");
  TH1F *ratioEle2pt_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_Ele2Pt-JPsiToLowQ2");
  TH1F *ratioElesDr_fromFile = (TH1F*)ratioFileMCBins->Get("ratio_ElesDr-JPsiToLowQ2");


  // Original ntuple   
  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(TString("./")+TString(dataFilename));
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("mytreefit");
  } else {
    cout << "File " << dataFilename << " not existing !" << endl;
    return;
  }
  
  fileOrig->cd();
  if (!treeOrig) {
    cout << "mytreefit not existing !" << endl; 
    return;    
  }

  treeOrig->SetMakeClass(0);
  cout << "TreeOrig->Size = "<< treeOrig->GetEntries() << endl;
  
  // number of entries saved in the first tree
  int nentriesOrig = treeOrig->GetEntries();   
  
  // Tree for the final format
  TFile *fileNew = TFile::Open(TString("./Formatted_Check_")+TString(dataFilename),"recreate");
  fileNew->ls();
  fileNew->cd();
  TTree *treeNew = new TTree("mytreefit","mytreefit with weight");
  treeNew->SetAutoSave(-99999999999);
  treeNew->SetAutoFlush(-99999999999);

  std::vector<TTree*> trees; 
  trees.push_back(treeNew);

  // original tree leaves
  Double_t         xgb = 0.;
  Double_t         Bmass = 0.;
  Double_t         Mll = 0.;
  Double_t         Npv = 0.;
  Double_t         Bprob = 0.;
  Double_t         BsLxy = 0.;
  Double_t         L1pt = 0.;
  Double_t         L2pt = 0.;
  Double_t         Kpt = 0.;
  Double_t         Bcos = 0.;
  Double_t         LKdz = 0.;
  Double_t         L1L2dr = 0.;
  Double_t         LKdr = 0.;
  Double_t         L2id = 0.;
  Double_t         Kiso = 0.;
  Double_t         BBDphi = 0.;
  Double_t         BTrkdxy2 = 0.;
  Double_t         L1id = 0.;
  Double_t         L1iso = 0.;
  Double_t         KLmassD0 = 0.;
  Double_t         Passymetry = 0.;
  Double_t         Kip3d = 0.;
  Double_t         Kip3dErr = 0.;
  Double_t         Bpt = 0.;
  Double_t         Mu7_IP4 = 0.;
  Double_t         Mu8_IP3 = 0.;
  Double_t         Mu8_IP5 = 0.;
  Double_t         Mu8_IP6 = 0.;
  Double_t         Mu9_IP5 = 0.;
  Double_t         Mu9_IP6 = 0.;
  Double_t         Mu12_IP6 = 0.;

  // List of branches - original tree
  TBranch        *b_xgb;   //!
  TBranch        *b_Bmass;   //!
  TBranch        *b_Mll;   //!
  TBranch        *b_Npv;   //!
  TBranch        *b_Bprob;   //!
  TBranch        *b_BsLxy;   //!
  TBranch        *b_L1pt;   //!
  TBranch        *b_L2pt;   //!
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
  TBranch        *b_Bpt;   //!
  TBranch        *b_Mu7_IP4;   //!
  TBranch        *b_Mu8_IP3;   //!
  TBranch        *b_Mu8_IP5;   //!
  TBranch        *b_Mu8_IP6;   //!
  TBranch        *b_Mu9_IP5;   //!
  TBranch        *b_Mu9_IP6;   //!
  TBranch        *b_Mu12_IP6;   //!

  // Set branch addresses and branch pointers 
  treeOrig->SetBranchAddress("xgb", &xgb, &b_xgb);
  treeOrig->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
  treeOrig->SetBranchAddress("Mll", &Mll, &b_Mll);
  treeOrig->SetBranchAddress("Npv", &Npv, &b_Npv);
  treeOrig->SetBranchAddress("Bprob", &Bprob, &b_Bprob);
  treeOrig->SetBranchAddress("BsLxy", &BsLxy, &b_BsLxy);
  treeOrig->SetBranchAddress("L1pt", &L1pt, &b_L1pt);
  treeOrig->SetBranchAddress("L2pt", &L2pt, &b_L2pt);
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
  treeOrig->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
  treeOrig->SetBranchAddress("Mu7_IP4", &Mu7_IP4, &b_Mu7_IP4);
  treeOrig->SetBranchAddress("Mu8_IP3", &Mu8_IP3, &b_Mu8_IP3);
  treeOrig->SetBranchAddress("Mu8_IP5", &Mu8_IP5, &b_Mu8_IP5);
  treeOrig->SetBranchAddress("Mu8_IP6", &Mu8_IP6, &b_Mu8_IP6);
  treeOrig->SetBranchAddress("Mu9_IP5", &Mu9_IP5, &b_Mu9_IP5);
  treeOrig->SetBranchAddress("Mu9_IP6", &Mu9_IP6, &b_Mu9_IP6);
  treeOrig->SetBranchAddress("Mu12_IP6", &Mu12_IP6, &b_Mu12_IP6);

  // New variables
  Double_t weight;

  // New branches
  for(int i=0; i<(int)trees.size();i++) {
    TTree *theTreeNew = trees[i];
    theTreeNew->Branch("weight", &weight, "weight/D");
    theTreeNew->Branch("xgb", &xgb, "xgb/D");
    theTreeNew->Branch("Bmass", &Bmass, "Bmass/D");
    theTreeNew->Branch("Mll", &Mll, "Mll/D");
    theTreeNew->Branch("Npv", &Npv, "Npv/D");
    theTreeNew->Branch("Bprob", &Bprob, "Bprob/D");
    theTreeNew->Branch("BsLxy", &BsLxy, "BsLxy/D");
    theTreeNew->Branch("L1pt", &L1pt, "L1pt/D");
    theTreeNew->Branch("L2pt", &L2pt, "L2pt/D");
    theTreeNew->Branch("Kpt", &Kpt, "Kpt/D");
    theTreeNew->Branch("Bcos", &Bcos, "Bcos/D");
    theTreeNew->Branch("LKdz", &LKdz, "LKdz/D");
    theTreeNew->Branch("L1L2dr", &L1L2dr, "L1L2dr/D");
    theTreeNew->Branch("LKdr", &LKdr, "LKdr/D");
    theTreeNew->Branch("L2id", &L2id, "L2id/D");
    theTreeNew->Branch("Kiso", &Kiso, "Kiso/D");
    theTreeNew->Branch("BBDphi", &BBDphi, "BBDphi/D");
    theTreeNew->Branch("BTrkdxy2", &BTrkdxy2, "BTrkdxy2/D");
    theTreeNew->Branch("L1id", &L1id, "L1id/D");
    theTreeNew->Branch("L1iso", &L1iso, "L1iso/D");
    theTreeNew->Branch("KLmassD0", &KLmassD0, "KLmassD0/D");
    theTreeNew->Branch("Passymetry", &Passymetry, "Passymetry/D");
    theTreeNew->Branch("Kip3d", &Kip3d, "Kip3d/D");
    theTreeNew->Branch("Kip3dErr", &Kip3dErr, "Kip3dErr/D");
    theTreeNew->Branch("Bpt", &Bpt, "Bpt/D");
    theTreeNew->Branch("Mu7_IP4", &Mu7_IP4, "Mu7_IP4/D");
    theTreeNew->Branch("Mu8_IP3", &Mu8_IP3, "Mu8_IP3/D");
    theTreeNew->Branch("Mu8_IP5", &Mu8_IP5, "Mu8_IP5/D");
    theTreeNew->Branch("Mu8_IP6", &Mu8_IP6, "Mu8_IP6/D");
    theTreeNew->Branch("Mu9_IP5", &Mu9_IP5, "Mu9_IP5/D");
    theTreeNew->Branch("Mu9_IP6", &Mu9_IP6, "Mu9_IP6/D");
    theTreeNew->Branch("Mu12_IP6", &Mu12_IP6, "Mu12_IP6/D");
  }


  // ----------------------------------------------------------------------
  // Load weights: ele1 pt
  for (int i = 0; i<ratioEle1pt_fromFile->GetNbinsX(); i++) {

    double myweightEle1pt=ratioEle1pt_fromFile->GetBinContent(i+1);
    anweightEle1pt_.push_back(myweightEle1pt);              

    double mylowedgeEle1pt=ratioEle1pt_fromFile->GetBinLowEdge(i+1); 
    lowedgeEle1pt_.push_back(mylowedgeEle1pt);
  }

  // Load weights: ele2 pt
  for (int i = 0; i<ratioEle2pt_fromFile->GetNbinsX(); i++) {

    double myweightEle2pt=ratioEle2pt_fromFile->GetBinContent(i+1);
    anweightEle2pt_.push_back(myweightEle2pt);              

    double mylowedgeEle2pt=ratioEle2pt_fromFile->GetBinLowEdge(i+1); 
    lowedgeEle2pt_.push_back(mylowedgeEle2pt);
  }

  // Load weights: dR(ele1, ele2(
  for (int i = 0; i<ratioElesDr_fromFile->GetNbinsX(); i++) {

    double myweightElesDr=ratioElesDr_fromFile->GetBinContent(i+1);
    anweightElesDr_.push_back(myweightElesDr);              

    double mylowedgeElesDr=ratioElesDr_fromFile->GetBinLowEdge(i+1); 
    lowedgeElesDr_.push_back(mylowedgeElesDr);
  }


  // Loop over events
  cout << "Now preparing the new tree" << endl;
  for(int i=0; i<nentriesOrig; i++) {
    
    if (i%10000 == 0) std::cout << ">>> Event # " << i << " / " << nentriesOrig << " entries" << std::endl; 
    treeOrig->GetEntry(i);

    // compute and add the weight
    double OneOverWeightEle1pt = GetEle1ptWeight(L1pt);
    double OneOverWeightEle2pt = GetEle2ptWeight(L2pt);
    double OneOverWeightElesDr = GetElesDrWeight(L1L2dr);
    double thisEventWeightEle1pt, thisEventWeightEle2pt, thisEventWeightElesDr;
    if (OneOverWeightEle1pt!=0) 
      thisEventWeightEle1pt = 1./OneOverWeightEle1pt;
    else
      thisEventWeightEle1pt = 0;

    if (OneOverWeightEle2pt!=0) 
      thisEventWeightEle2pt = 1./OneOverWeightEle2pt;
    else
      thisEventWeightEle2pt = 0;

    if (OneOverWeightElesDr!=0) 
      thisEventWeightElesDr = 1./OneOverWeightElesDr;
    else
      thisEventWeightElesDr = 0;

    if (theVariable==0) weight = thisEventWeightEle1pt;
    if (theVariable==1) weight = thisEventWeightEle2pt;
    if (theVariable==2) weight = thisEventWeightElesDr;

    treeNew->Fill();
  }

  // new format
  treeNew->Write();
  fileNew->Close();
  fileNew->ls();
  
  fileOrig->cd();
  fileOrig->Close();  
}

double GetEle1ptWeight(double ele1pt) {

  int thesizem1 = lowedgeEle1pt_.size()-1;
  
  double weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedgeEle1pt_[i]<=ele1pt && lowedgeEle1pt_[i+1]>ele1pt) weight = anweightEle1pt_[i];  
  }
  
  return weight;
}

double GetEle2ptWeight(double ele2pt) {

  int thesizem1 = lowedgeEle2pt_.size()-1;
  
  double weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedgeEle2pt_[i]<=ele2pt && lowedgeEle2pt_[i+1]>ele2pt) weight = anweightEle2pt_[i];  
  }
  
  return weight;
}

double GetElesDrWeight(double elesdr) {

  int thesizem1 = lowedgeElesDr_.size()-1;
  
  double weight=1;
  for (int i = 0; i<thesizem1; i++) {   
    if (lowedgeElesDr_[i]<=elesdr && lowedgeElesDr_[i+1]>elesdr) weight = anweightElesDr_[i];  
  }
  
  return weight;
}



