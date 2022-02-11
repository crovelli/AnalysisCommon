// Based on ROOT/tutorials/roostats/rs301_splot.C

// To be run on data tnp formatted ntuples 
// - Fit B mass as a reference distribution
// - Extract analysis BDT distributions for signal and background with the sPlots technique

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsArg.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooKeysPdf.h"
#include "RooFormulaVar.h"
#include "RooStats/SPlot.h"
#include "DoubleSidedCB2.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

#include <iostream>

// use this order for safety on library loading
using namespace RooFit;
using namespace std;

// see below for implementation
void DoSPlot(RooWorkspace*);                         
void MakePlots(RooWorkspace*);                       
void MakeHistos(RooWorkspace*, float);
void getDataSet(const char *, RooWorkspace*, float, float, float);   

void lowq2Bin(const char* workspacefile, const char* roottreefile, float bdtCut)
{
  // Define q^2 ranges
  float q2inf = 1.05;
  float q2sup = 2.45;
  // if (q2bin==0) { q2inf=1.05; q2sup=2.45; }
  // if (q2bin==1) { q2inf=2.90; q2sup=3.20; }
  // if (q2bin==2) { q2inf=3.55; q2sup=3.80; }
  // if (q2bin==3) { q2inf=4.00; q2sup=4.80; }

  // Load the existing workspace 
  cout << endl;
  cout << endl;
  TFile *wsfile = new TFile(workspacefile);
  RooWorkspace* ws = (RooWorkspace*)wsfile->Get("wspace");
  ws->Print();

  // Add dataset from converted root tree. 
  // Should be identical to the one in the workspace + extra variables needed for splots
  getDataSet(roottreefile, ws, bdtCut, q2inf, q2sup);

  // make a new dataset with sWeights added for every event.
  DoSPlot(ws);

  // Make some plots showing the discriminating variable and
  // the control variable after unfolding.
  MakePlots(ws);
  
  // Save variables in histos
  MakeHistos(ws, bdtCut);

  // cleanup
  delete ws;
}

// Load everything from WS, compute s-plots and make the final fit
void DoSPlot(RooWorkspace* wspace){

  cout << endl;
  cout << "=================== Load dataset (reading from WS) ==================" << endl;   
  RooAbsData *mydata = wspace->data("mydata");
  RooDataSet *mydatads = (RooDataSet*)mydata;
  mydata->Print(); 
  mydatads->Print(); 
  
  cout << endl;
  cout << "=================== AddModel (reading from WS) ==================" << endl;
  cout << endl;


  // make a RooRealVar for the observable
  RooRealVar *Bmass = wspace->var("Bmass");
  RooRealVar *x = wspace->var("x");


  // --------------------------------------
  // signal model
  std::cout << "Signal model" << std::endl;
  //
  RooRealVar *mean   = wspace->var("mean");  
  RooRealVar *width  = wspace->var("width");  
  RooRealVar *alpha1 = wspace->var("alpha1");  
  RooRealVar *alpha2 = wspace->var("alpha2");  
  RooRealVar *n1 = wspace->var("n1");  
  RooRealVar *n2 = wspace->var("n2");  
  DoubleSidedCB2 mysig("mysig", "mysig", *x, *mean, *width, *alpha1, *n1, *alpha2, *n2);
  RooRealVar *nsig = wspace->var("nsig");


  // --------------------------------------
  // background model
  std::cout << "OtherB model" << std::endl;
  //
  RooRealVar *nbkg = wspace->var("nbkg");
  RooRealVar *exp_alpha = wspace->var("exp_alpha");  
  RooAbsPdf  *bkg = wspace->pdf("bkg");


  // --------------------------------------
  // background model: K*Jpsi
  std::cout << "K*Jpsi background model" << std::endl;
  //
  RooRealVar *nkjpsi = wspace->var("nkjpsi");
  RooAbsPdf  *kjpsi  = wspace->pdf("kjpsi");


  // --------------------------------------
  // Combined model
  std::cout << "Full model" << std::endl; 
  RooAddPdf mymodel("mymodel","signal+background models",
		    RooArgList(mysig, *bkg, *kjpsi),
		    RooArgList(*nsig, *nbkg, *nkjpsi)); 
  mymodel.Print();

  wspace->import(mymodel);


  // -----------------------------------------
  cout << endl;
  cout << "=================== Now sPlots ==================" << endl;   
  cout << endl;

  // fit the model to the data.
  mymodel.fitTo(*mydata, Extended() );
  RooMsgService::instance().setSilentMode(false);


  // The sPlot technique requires that we fix the parameters
  // of the model that are not yields after doing the fit.
  mean->setConstant();
  width->setConstant();
  alpha1->setConstant();
  alpha2->setConstant();
  n1->setConstant();
  n2->setConstant();
  exp_alpha->setConstant();


  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  RooStats::SPlot* sData = new RooStats::SPlot("sData","A SPlot", *mydatads, &mymodel, RooArgList(*nsig, *nbkg, *nkjpsi) );

  // Check that our weights have the desired properties
  std::cout << "Check SWeights:" << std::endl;
  
  std::cout << std::endl <<  "Signal yield is " << nsig->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("nsig") << std::endl;
  
  std::cout << std::endl <<  "OtherB yield is " << nbkg->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("nbkg") << std::endl;
  
  std::cout << std::endl <<  "K*JPsi background fraction is " << nkjpsi->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("nkjpsi") << std::endl;
  
  std::cout << std::endl;
  
  // import this new dataset with sWeights
  std::cout << "import new dataset with sWeights" << std::endl;
  wspace->import(*mydata, Rename("mydataWithSWeights"));
}

// Control plots
void MakePlots(RooWorkspace* ws){
  
  cout << endl;
  cout << "====================== MakePlots ==================" << endl;
  cout << endl;

  // get what we need out of the workspace
  RooAbsPdf* mymodel = ws->pdf("mymodel");
  mymodel->Print();

  RooAbsPdf* mysig   = ws->pdf("mysig");
  RooAbsPdf* bkg     = ws->pdf("bkg");
  RooAbsPdf* kjpsi   = ws->pdf("kjpsi");
  RooRealVar *x      = ws->var("x");  
  RooRealVar *xgb    = ws->var("xgb");
  x->SetTitle("m(K^{+}e^{+}e^{-}) [GeV]");
  xgb->SetTitle("BDT output");

  // note, we get the dataset with sWeights
  RooDataSet* mydata = (RooDataSet*) ws->data("mydataWithSWeights");


  // ------------------------------------------
  // Fit plot
  TCanvas* cdata = new TCanvas("cdata","cdata", 800, 600);
  cdata->SetGridx();
  cdata->SetGridy();
  RooPlot* frame = x->frame() ;
  mydata->plotOn(frame, RooFit::Binning(20)) ;
  mymodel->plotOn(frame, LineColor(kRed)) ;
  mymodel->plotOn(frame, Components(*mysig), Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(4),  LineWidth(2), LineStyle(kDashed));
  mymodel->plotOn(frame, Components(*bkg),   Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(49), LineWidth(2), LineStyle(kDashed));
  mymodel->plotOn(frame, Components(*kjpsi), Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(30),  LineWidth(2), LineStyle(kDashed));
  frame->SetTitle("Fit of model to discriminating variable");
  cdata->cd(1);
  frame->Draw();
  cdata->SaveAs("Fit.png");


  // ------------------------------------------
  // Now use the sWeights to show our variable distribution for B and background.
   
  // create weighted data set
  RooDataSet * dataw_signal = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"nsig_sw") ;
  RooDataSet * dataw_bkg    = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"nbkg_sw") ;
  RooDataSet * dataw_kjpsi  = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"nkjpsi_sw") ;

  TCanvas* csignalMass = new TCanvas("csignalMass","sPlot demo mass, signal", 800, 600);  
  RooPlot* fsignalMass = x->frame() ;   
  dataw_signal->plotOn(fsignalMass, DataError(RooAbsData::SumW2) ) ;
  fsignalMass->SetTitle("Mass");
  csignalMass->cd(1);
  fsignalMass->Draw() ;
  csignalMass->SaveAs("sPlotMassSignal.png");

  TCanvas* cbkgMass = new TCanvas("cbkgMass","sPlot demo mass, bkg", 800, 600);  
  RooPlot* fbkgMass = x->frame() ;   
  dataw_bkg->plotOn(fbkgMass, DataError(RooAbsData::SumW2) ) ;
  fbkgMass->SetTitle("Mass");
  cbkgMass->cd(1);
  fbkgMass->Draw() ;
  cbkgMass->SaveAs("sPlotMassBkg.png");

  TCanvas* ckjpsiMass = new TCanvas("ckjpsiMass","sPlot demo mass, kjpsi", 800, 600);  
  RooPlot* fkjpsiMass = x->frame() ;  
  dataw_kjpsi->plotOn(fkjpsiMass, DataError(RooAbsData::SumW2) ) ;
  fkjpsiMass->SetTitle("Mass");
  ckjpsiMass->cd(1);
  fkjpsiMass->Draw() ;
  ckjpsiMass->SaveAs("sPlotMassKjpsi.png");

  TCanvas* csignalXgb = new TCanvas("csignalXgb","sPlot demo bdt, signal", 800, 600);  
  RooPlot* fsignalXgb = xgb->frame() ;   
  dataw_signal->plotOn(fsignalXgb, DataError(RooAbsData::SumW2) ) ;
  fsignalXgb->SetTitle("Analysis BDT");
  csignalXgb->cd(1);
  fsignalXgb->Draw() ;
  csignalXgb->SaveAs("sPlotXgbSignal.png");

  TCanvas* cbkgXgb = new TCanvas("cbkgXgb","sPlot demo bdt, otherb", 800, 600);  
  RooPlot* fbkgXgb = xgb->frame() ;   
  dataw_bkg->plotOn(fbkgXgb, DataError(RooAbsData::SumW2) ) ;
  fbkgXgb->SetTitle("Analysis BDT");
  cbkgXgb->cd(1);
  fbkgXgb->Draw() ;
  cbkgXgb->SaveAs("sPlotXgbOtherB.png");

  TCanvas* ckjpsiXgb = new TCanvas("ckjpsiXgb","sPlot demo bdt, kjpsi", 800, 600);  
  RooPlot* fkjpsiXgb = xgb->frame() ;  
  dataw_kjpsi->plotOn(fkjpsiXgb, DataError(RooAbsData::SumW2) ) ;
  fkjpsiXgb->SetTitle("Analysis BDT");
  ckjpsiXgb->cd(1);
  fkjpsiXgb->Draw() ;
  ckjpsiXgb->SaveAs("sPlotXgbKjpsi.png");
}


// Control plots
void MakeHistos(RooWorkspace* ws, float bdtCut){
  
  gStyle->SetOptStat(0);

  std::cout << std::endl;
  std::cout << "Save histos" << std::endl;

  RooRealVar* xgb      = ws->var("xgb");                
  RooRealVar* L2id     = ws->var("L2id");            
  RooRealVar* L1id     = ws->var("L1id");            
  RooRealVar* Bprob    = ws->var("Bprob");           
  RooRealVar* BsLxy    = ws->var("BsLxy");           
  RooRealVar* Bcos     = ws->var("Bcos");            
  RooRealVar* L1pt     = ws->var("L1pt");         
  RooRealVar* L2pt     = ws->var("L2pt");         
  RooRealVar* Kpt      = ws->var("Kpt");          
  RooRealVar* LKdz     = ws->var("LKdz");            
  RooRealVar* L1L2dr   = ws->var("L1L2dr");          
  RooRealVar* LKdr     = ws->var("LKdr");            
  RooRealVar* L1iso    = ws->var("L1iso");       
  RooRealVar* Kiso     = ws->var("Kiso");        
  RooRealVar* Passymetry = ws->var("Passymetry");                
  RooRealVar* Kip3d    = ws->var("Kip3d");         

  // note, we get the dataset with sWeights
  RooDataSet* mydata = (RooDataSet*) ws->data("mydataWithSWeights");

  // create weighted data set
  RooDataSet * mydataw_sgn = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"nsig_sw") ;
  mydataw_sgn->Print();
  
  // convert to TH1
  float theinf = bdtCut-2;
  float thedelta = 12.5-theinf;
  int thebin = thedelta/0.5;
  TH1 *h1_xgb    = mydataw_sgn->createHistogram("h1_xgb",*xgb,Binning(thebin,theinf,12.5)); 
  TH1 *h1_L1id   = mydataw_sgn->createHistogram("h1_L1id",*L1id, Binning(22,-4.,7.)); 
  TH1 *h1_L2id   = mydataw_sgn->createHistogram("h1_L2id",*L2id, Binning(22,-4.,7.)); 
  TH1 *h1_Bprob  = mydataw_sgn->createHistogram("h1_Bprob",*Bprob, Binning(10,0.,1.)); 
  TH1 *h1_BsLxy  = mydataw_sgn->createHistogram("h1_BsLxy",*BsLxy, Binning(10,0.,100.)); 
  TH1 *h1_Bcos   = mydataw_sgn->createHistogram("h1_Bcos",*Bcos, Binning(10,0.99,1.)); 
  TH1 *h1_L1pt   = mydataw_sgn->createHistogram("h1_L1pt",*L1pt, Binning(60,0.,30)); 
  TH1 *h1_L2pt   = mydataw_sgn->createHistogram("h1_L2pt",*L2pt, Binning(40,0.,20)); 
  TH1 *h1_Kpt    = mydataw_sgn->createHistogram("h1_Kpt",*Kpt, Binning(40,0.,20)); 
  TH1 *h1_LKdz   = mydataw_sgn->createHistogram("h1_LKdz",*LKdz, Binning(20,0.,1.)); 
  TH1 *h1_L1L2dr = mydataw_sgn->createHistogram("h1_L1L2dr",*L1L2dr, Binning(20,0.,2.)); 
  TH1 *h1_LKdr   = mydataw_sgn->createHistogram("h1_LKdr",*LKdr, Binning(40,0.,1.)); 
  TH1 *h1_L1iso  = mydataw_sgn->createHistogram("h1_L1iso",*L1iso, Binning(30,0.,30.)); 
  TH1 *h1_Kiso   = mydataw_sgn->createHistogram("h1_Kiso",*Kiso, Binning(30,0.,30.)); 
  TH1 *h1_Passymetry = mydataw_sgn->createHistogram("h1_Passymetry",*Passymetry, Binning(20,-1.,1.)); 
  TH1 *h1_Kip3d  = mydataw_sgn->createHistogram("h1_Kip3d",*Kip3d, Binning(40,-0.2,0.2)); 

  TFile myFileSPlots("myFileSPlots.root","RECREATE");
  myFileSPlots.cd();
  h1_xgb   -> Write();
  h1_L1id  -> Write();
  h1_L2id  -> Write();
  h1_Bprob -> Write();
  h1_BsLxy -> Write();
  h1_Bcos  -> Write();
  h1_L1pt  -> Write();
  h1_L2pt  -> Write();
  h1_Kpt   -> Write();
  h1_LKdz  -> Write();
  h1_L1L2dr -> Write();
  h1_LKdr  -> Write();
  h1_L1iso -> Write();
  h1_Kiso  -> Write();
  h1_Passymetry -> Write();
  h1_Kip3d -> Write();
}

// Convert ROOT tree in RooDataset
void getDataSet(const char *rootfile, RooWorkspace *ws, float bdtCut, float q2min, float q2max) {    

  cout << endl;  
  cout << "============= getDataSet =================" << endl;
  cout << "roofitting file " << rootfile << endl;
  
  // fit variables
  RooRealVar Bmass("Bmass", "M_B", 4.7, 5.7,"GeV");
  RooRealVar x("x", "x", 4.7, 5.7,"GeV");
  // 
  // discriminating variables and BDT output   
  RooRealVar Mll("Mll", "Mll", q2min, q2max, "");           
  RooRealVar KLmassD0("KLmassD0", "KLmassD0", 2., 100000., "");
  //
  // BDT output
  RooRealVar xgb("xgb", "xgb", -15., 15., "");           

  // BDT inputs
  RooRealVar L2id("L2id", "L2id", -20.5, 20.5, "");      
  RooRealVar L1id("L1id", "L1id", -20.5, 20.5, "");      
  RooRealVar Bprob("Bprob","Bprob",-1.,1.,"");           
  RooRealVar BsLxy("BsLxy","BsLxy",0.,500.,"");          
  RooRealVar Bcos("Bcos","Bcos",-1.,1.,"");         
  RooRealVar L1pt("L1pt","L1pt", 0., 5000.,"");       
  RooRealVar L2pt("L2pt","L2pt", 0., 5000.,"");       
  RooRealVar Kpt("Kpt","Kpt",   0., 5000.,"");        
  RooRealVar LKdz("LKdz", "LKdz",     0., 50.,"");                 
  RooRealVar L1L2dr("L1L2dr","L1L2dr", 0., 50.,"");                 
  RooRealVar LKdr("LKdr", "LKdr",     0., 50, "");                 
  RooRealVar L1iso("L1iso","L1iso", 0., 10000.,"");      
  RooRealVar Kiso("Kiso", "Kiso",  0., 10000.,"");      
  RooRealVar Passymetry("Passymetry", "Passymetry", -1., 1., "");
  RooRealVar Kip3d("Kip3d", "Kip3d", -1000, 1000, ""); 

  RooArgSet setall(Bmass,x,Mll,L2id,L1id,xgb,KLmassD0);
  setall.add(Bprob);
  setall.add(BsLxy);
  setall.add(Bcos);
  setall.add(L1pt);
  setall.add(L2pt);
  setall.add(Kpt);
  setall.add(LKdz);
  setall.add(L1L2dr);
  setall.add(LKdr);
  setall.add(L1iso);
  setall.add(Kiso);
  setall.add(Passymetry);
  setall.add(Kip3d);

  TFile *file = TFile::Open(rootfile);
  TTree *tree = (TTree*)file->Get("mytreefit");
  
  RooDataSet *mydata = new RooDataSet("mydata","mydata",tree,setall,0); 

  // Inclusive
  string theBdtCut = to_string(bdtCut);
  string theQ2minCut = to_string(q2min);
  string theQ2maxCut = to_string(q2max);
  TString theSel ="Bmass>4.7 && Bmass<5.7 && Mll>"+theQ2minCut+" && Mll<"+theQ2maxCut+" && xgb>"+theBdtCut+ " && KLmassD0 > 2.0";

  cout << endl;
  cout << "theSel = " << theSel << endl;
  mydata = (RooDataSet*)mydata->reduce(theSel);
  cout << endl;
  mydata->Print();

  ws->import(*mydata);
  ws->Print();
}

