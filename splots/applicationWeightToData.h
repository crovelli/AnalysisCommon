//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 23 12:28:39 2022 by ROOT version 6.24/06
// from TTree mytreefit/mytreefit with weight
// found on file: forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_data_mvaCut0_withElesDrWeight.root
//////////////////////////////////////////////////////////

#ifndef applicationWeightToData_h
#define applicationWeightToData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class applicationWeightToData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        weight;
   Double_t        xgb;
   Double_t        Bmass;
   Double_t        Mll;
   Double_t        Npv;
   Double_t        Bprob;
   Double_t        BsLxy;
   Double_t        L1pt;
   Double_t        L2pt;
   Double_t        Kpt;
   Double_t        Bcos;
   Double_t        LKdz;
   Double_t        L1L2dr;
   Double_t        LKdr;
   Double_t        L2id;
   Double_t        Kiso;
   Double_t        BBDphi;
   Double_t        BTrkdxy2;
   Double_t        L1id;
   Double_t        L1iso;
   Double_t        KLmassD0;
   Double_t        Passymetry;
   Double_t        Kip3d;
   Double_t        Kip3dErr;
   Double_t        Bpt;
   Double_t        Mu7_IP4;
   Double_t        Mu8_IP3;
   Double_t        Mu8_IP5;
   Double_t        Mu8_IP6;
   Double_t        Mu9_IP5;
   Double_t        Mu9_IP6;
   Double_t        Mu12_IP6;

   // List of branches
   TBranch        *b_weight;   //!
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

   applicationWeightToData(TTree *tree=0);
   virtual ~applicationWeightToData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int theVariable);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef applicationWeightToData_cxx
applicationWeightToData::applicationWeightToData(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withEle1PtWeight.root");
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withEle2PtWeight.root");
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withElesDrWeight.root");
    if (!f || !f->IsOpen()) {
      //f = new TFile("forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withEle1PtWeight.root");
      //f = new TFile("forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withEle2PtWeight.root");
      f = new TFile("forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withElesDrWeight.root");
    }
    f->GetObject("mytreefit",tree);
    
  }
  Init(tree);
}

applicationWeightToData::~applicationWeightToData()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t applicationWeightToData::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t applicationWeightToData::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
     fCurrent = fChain->GetTreeNumber();
     Notify();
   }
   return centry;
}

void applicationWeightToData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("xgb", &xgb, &b_xgb);
   fChain->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
   fChain->SetBranchAddress("Mll", &Mll, &b_Mll);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("Bprob", &Bprob, &b_Bprob);
   fChain->SetBranchAddress("BsLxy", &BsLxy, &b_BsLxy);
   fChain->SetBranchAddress("L1pt", &L1pt, &b_L1pt);
   fChain->SetBranchAddress("L2pt", &L2pt, &b_L2pt);
   fChain->SetBranchAddress("Kpt", &Kpt, &b_Kpt);
   fChain->SetBranchAddress("Bcos", &Bcos, &b_Bcos);
   fChain->SetBranchAddress("LKdz", &LKdz, &b_LKdz);
   fChain->SetBranchAddress("L1L2dr", &L1L2dr, &b_L1L2dr);
   fChain->SetBranchAddress("LKdr", &LKdr, &b_LKdr);
   fChain->SetBranchAddress("L2id", &L2id, &b_L2id);
   fChain->SetBranchAddress("Kiso", &Kiso, &b_Kiso);
   fChain->SetBranchAddress("BBDphi", &BBDphi, &b_BBDphi);
   fChain->SetBranchAddress("BTrkdxy2", &BTrkdxy2, &b_BTrkdxy2);
   fChain->SetBranchAddress("L1id", &L1id, &b_L1id);
   fChain->SetBranchAddress("L1iso", &L1iso, &b_L1iso);
   fChain->SetBranchAddress("KLmassD0", &KLmassD0, &b_KLmassD0);
   fChain->SetBranchAddress("Passymetry", &Passymetry, &b_Passymetry);
   fChain->SetBranchAddress("Kip3d", &Kip3d, &b_Kip3d);
   fChain->SetBranchAddress("Kip3dErr", &Kip3dErr, &b_Kip3dErr);
   fChain->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
   fChain->SetBranchAddress("Mu7_IP4", &Mu7_IP4, &b_Mu7_IP4);
   fChain->SetBranchAddress("Mu8_IP3", &Mu8_IP3, &b_Mu8_IP3);
   fChain->SetBranchAddress("Mu8_IP5", &Mu8_IP5, &b_Mu8_IP5);
   fChain->SetBranchAddress("Mu8_IP6", &Mu8_IP6, &b_Mu8_IP6);
   fChain->SetBranchAddress("Mu9_IP5", &Mu9_IP5, &b_Mu9_IP5);
   fChain->SetBranchAddress("Mu9_IP6", &Mu9_IP6, &b_Mu9_IP6);
   fChain->SetBranchAddress("Mu12_IP6", &Mu12_IP6, &b_Mu12_IP6);
   Notify();
}

Bool_t applicationWeightToData::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void applicationWeightToData::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t applicationWeightToData::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef applicationWeightToData_cxx
