//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  3 16:15:10 2022 by ROOT version 6.24/06
// from TTree mytreefit/mytreefit
// found on file: ottoPFPF/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_0_MCres.root
//////////////////////////////////////////////////////////

#ifndef distribMC_h
#define distribMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class distribMC {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  Double_t        xgb;
  Float_t         Bmass;
  Float_t         Mll;
  Float_t         Npv;
  Float_t         Bprob;
  Float_t         BsLxy;
  Float_t         L1pt;
  Float_t         L2pt;
  Float_t         Bpt;
  Float_t         Kpt;
  Float_t         Bcos;
  Float_t         LKdz;
  Float_t         L1L2dr;
  Float_t         LKdr;
  Float_t         L2id;
  Float_t         Kiso;
  Float_t         BBDphi;
  Float_t         BTrkdxy2;
  Float_t         L1id;
  Float_t         L1iso;
  Float_t         KLmassD0;
  Float_t         Passymetry;
  Float_t         Kip3d;
  Float_t         Kip3dErr;
  Float_t         Mu8_IP3;
  Float_t         Mu8_IP5;
  Float_t         Mu8_IP6;
  Float_t         Mu9_IP5;
  Float_t         Mu9_IP6;
  Float_t         Mu12_IP6;
  
  // List of branches
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
  TBranch        *b_Mu8_IP3;   //!
  TBranch        *b_Mu8_IP5;   //!
  TBranch        *b_Mu8_IP6;   //!
  TBranch        *b_Mu9_IP5;   //!
  TBranch        *b_Mu9_IP6;   //!
  TBranch        *b_Mu12_IP6;   //!
  
  distribMC(TTree *tree=0);
  virtual ~distribMC();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(int q2bin, float bdtCut, int isPFPF);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef distribMC_cxx
distribMC::distribMC(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {          
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCres.root");      // JPsi   PFPF
    //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC.root");         // Low-q2 PFPF
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCPsi2S.root");    // Psi2s  PFPF
    // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MCres.root");       // JPsi PFLP
    if (!f || !f->IsOpen()) {
      //f = new TFile("ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCres.root");           // JPsi PFPF   
      //f = new TFile("ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC.root");              // Low-q2 PFPF     
      f = new TFile("ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCPsi2S.root");         // Psi2s  PFPF   
      // f = new TFile("ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MCres.root");            // JPsi PFLP  
    }
    f->GetObject("mytreefit",tree);
    
  }
  Init(tree);
}

distribMC::~distribMC()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t distribMC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t distribMC::LoadTree(Long64_t entry)
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

void distribMC::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("xgb", &xgb, &b_xgb);
  fChain->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
  fChain->SetBranchAddress("Mll", &Mll, &b_Mll);
  fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
  fChain->SetBranchAddress("Bprob", &Bprob, &b_Bprob);
  fChain->SetBranchAddress("BsLxy", &BsLxy, &b_BsLxy);
  fChain->SetBranchAddress("L1pt", &L1pt, &b_L1pt);
  fChain->SetBranchAddress("L2pt", &L2pt, &b_L2pt);
  fChain->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
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
  fChain->SetBranchAddress("Mu8_IP3", &Mu8_IP3, &b_Mu8_IP3);
  fChain->SetBranchAddress("Mu8_IP5", &Mu8_IP5, &b_Mu8_IP5);
  fChain->SetBranchAddress("Mu8_IP6", &Mu8_IP6, &b_Mu8_IP6);
  fChain->SetBranchAddress("Mu9_IP5", &Mu9_IP5, &b_Mu9_IP5);
  fChain->SetBranchAddress("Mu9_IP6", &Mu9_IP6, &b_Mu9_IP6);
  fChain->SetBranchAddress("Mu12_IP6", &Mu12_IP6, &b_Mu12_IP6);
  Notify();
}

Bool_t distribMC::Notify()
{
  return kTRUE;
}

void distribMC::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t distribMC::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef distribMC_cxx
