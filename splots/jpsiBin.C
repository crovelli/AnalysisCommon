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
#include "./SPlotJPsi.h"

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

void jpsiBin(const char* workspacefile, const char* roottreefile, float bdtCut)
{
  // Define q^2 ranges
  float q2inf = 2.90;
  float q2sup = 3.20;
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
  RooRealVar *n1 = wspace->var("n1");  
  RooRealVar *gauss_mean  = wspace->var("gauss_mean");  
  RooRealVar *gauss_width = wspace->var("gauss_width");  
  RooRealVar *frac = wspace->var("frac");  
  //
  RooAbsPdf *cb_signal = wspace->pdf("cb_signal");
  RooAbsPdf *g_signal  = wspace->pdf("g_signal");
  RooAbsPdf *signal    = wspace->pdf("signal");
  // 
  RooRealVar *nsignal = wspace->var("nsignal");


  // --------------------------------------
  // background model: otherB
  std::cout << "OtherB model" << std::endl;
  //
  RooRealVar *notherB = wspace->var("notherB");
  RooRealVar *exp_alpha_otherb = wspace->var("exp_alpha_otherb");  
  RooAbsPdf  *exp_otherb = wspace->pdf("exp_otherb");


  // --------------------------------------
  // background model: combinatorial
  std::cout << "Combinatorial background model" << std::endl;
  //
  RooRealVar *ncomb = wspace->var("ncomb");
  RooRealVar *exp_alpha_comb = wspace->var("exp_alpha_comb");
  RooAbsPdf  *exp_comb = wspace->pdf("exp_comb");


  // --------------------------------------
  // background model: K*Jpsi
  std::cout << "K*Jpsi background model" << std::endl;
  //
  RooRealVar *frac_partial  = wspace->var("frac_partial");
  RooAbsPdf  *kstarjpsi     = wspace->pdf("kstarjpsi");

    
  // --------------------------------------
  // Combined model
  std::cout << "Full model" << std::endl; 
  RooAbsPdf *model = wspace->pdf("model");
  model->Print();

  // -----------------------------------------
  cout << endl;
  cout << "=================== Now sPlots ==================" << endl;   
  cout << endl;


  // fit the model to the data.
  model->fitTo(*mydata, Extended() );
  RooMsgService::instance().setSilentMode(false);

  // The sPlot technique requires that we fix the parameters
  // of the model that are not yields after doing the fit.
  mean->setConstant(); 
  width->setConstant();
  alpha1->setConstant();
  n1->setConstant();
  gauss_mean->setConstant();
  gauss_width->setConstant();
  frac->setConstant();
  exp_alpha_otherb->setConstant();  
  exp_alpha_comb->setConstant();

  // Now we use the SPlot class to add SWeights to our data set
  // based on our model and our yield variables
  SPlotJPsi* sData = new SPlotJPsi("sData","A SPlot", *mydatads, model, RooArgList(*nsignal, *notherB, *ncomb, *frac_partial) );
  
  // Check that our weights have the desired properties
  std::cout << "Check SWeights:" << std::endl;
  
  std::cout << std::endl <<  "Signal yield is " << nsignal->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("nsignal") << std::endl;
  
  std::cout << std::endl <<  "OtherB yield is " << notherB->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("notherB") << std::endl;
  
  std::cout << std::endl <<  "Combinatorial background yield is " << ncomb->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("ncomb") << std::endl;
  
  std::cout << std::endl <<  "K*JPsi background fraction is " << frac_partial->getVal() 
	    << ".  From sWeights it is " << sData->GetYieldFromSWeight("frac_partial") << std::endl;
  
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
  RooAbsPdf* model      = ws->pdf("model");
  RooAbsPdf* signal     = ws->pdf("signal");
  RooAbsPdf* exp_otherb = ws->pdf("exp_otherb");
  RooAbsPdf* exp_comb   = ws->pdf("exp_comb");
  RooAbsPdf* kstarjpsi  = ws->pdf("kstarjpsi");
  RooRealVar *x         = ws->var("x");
  RooRealVar *xgb       = ws->var("xgb");
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
  mydata->plotOn(frame, RooFit::Binning(50)) ;
  model->plotOn(frame, LineColor(kRed)) ;
  model->plotOn(frame, Components(*signal),     Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(4),  LineWidth(2), LineStyle(kDashed));
  model->plotOn(frame, Components(*exp_otherb), Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(30), LineWidth(2), LineStyle(kDashed));
  model->plotOn(frame, Components(*exp_comb),   Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(49), LineWidth(2), LineStyle(kDashed));
  model->plotOn(frame, Components(*kstarjpsi),  Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(1),  LineWidth(2), LineStyle(kDashed));
  frame->SetTitle("Fit of model to discriminating variable");
  cdata->cd(1);
  frame->Draw();
  cdata->SaveAs("Fit.png");

  
  // ------------------------------------------
  // Now use the sWeights to show our variable distribution for B and background.
   
  // create weighted data set
  RooDataSet * dataw_signal  = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"nsignal_sw") ;
  RooDataSet * dataw_otherb  = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"notherB_sw") ;
  RooDataSet * dataw_comb    = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"ncomb_sw") ;
  RooDataSet * dataw_partial = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"frac_partial_sw") ;

  TCanvas* csignalMass = new TCanvas("csignalMass","sPlot demo mass, signal", 800, 600);  
  RooPlot* fsignalMass = x->frame() ;   
  dataw_signal->plotOn(fsignalMass, DataError(RooAbsData::SumW2) ) ;
  fsignalMass->SetTitle("Analysis BDT");
  csignalMass->cd(1);
  fsignalMass->Draw() ;
  csignalMass->SaveAs("sPlotMassSignal.png");

  TCanvas* cotherbMass = new TCanvas("cotherbMass","sPlot demo mass, otherb", 800, 600);  
  RooPlot* fotherbMass = x->frame() ;   
  dataw_otherb->plotOn(fotherbMass, DataError(RooAbsData::SumW2) ) ;
  fotherbMass->SetTitle("Analysis BDT");
  cotherbMass->cd(1);
  fotherbMass->Draw() ;
  cotherbMass->SaveAs("sPlotMassOtherB.png");

  TCanvas* ccombMass = new TCanvas("ccombMass","sPlot demo mass, comb", 800, 600);  
  RooPlot* fcombMass = x->frame() ;     
  dataw_comb->plotOn(fcombMass, DataError(RooAbsData::SumW2) ) ;
  fcombMass->SetTitle("Analysis BDT");
  ccombMass->cd(1);
  fcombMass->Draw() ;
  ccombMass->SaveAs("sPlotMassComb.png");

  TCanvas* cpartialMass = new TCanvas("cpartialMass","sPlot demo mass, partial", 800, 600);  
  RooPlot* fpartialMass = x->frame() ;  
  dataw_partial->plotOn(fpartialMass, DataError(RooAbsData::SumW2) ) ;
  fpartialMass->SetTitle("Analysis BDT");
  cpartialMass->cd(1);
  fpartialMass->Draw() ;
  cpartialMass->SaveAs("sPlotMassPartial.png");

  TCanvas* csignalXgb = new TCanvas("csignalXgb","sPlot demo bdt, signal", 800, 600);  
  RooPlot* fsignalXgb = xgb->frame() ;   
  dataw_signal->plotOn(fsignalXgb, DataError(RooAbsData::SumW2) ) ;
  fsignalXgb->SetTitle("Analysis BDT");
  csignalXgb->cd(1);
  fsignalXgb->Draw() ;
  csignalXgb->SaveAs("sPlotXgbSignal.png");

  TCanvas* cotherbXgb = new TCanvas("cotherbXgb","sPlot demo bdt, otherb", 800, 600);  
  RooPlot* fotherbXgb = xgb->frame() ;   
  dataw_otherb->plotOn(fotherbXgb, DataError(RooAbsData::SumW2) ) ;
  fotherbXgb->SetTitle("Analysis BDT");
  cotherbXgb->cd(1);
  fotherbXgb->Draw() ;
  cotherbXgb->SaveAs("sPlotXgbOtherB.png");

  TCanvas* ccombXgb = new TCanvas("ccombXgb","sPlot demo bdt, comb", 800, 600);  
  RooPlot* fcombXgb = xgb->frame() ;     
  dataw_comb->plotOn(fcombXgb, DataError(RooAbsData::SumW2) ) ;
  fcombXgb->SetTitle("Analysis BDT");
  ccombXgb->cd(1);
  fcombXgb->Draw() ;
  ccombXgb->SaveAs("sPlotXgbComb.png");

  TCanvas* cpartialXgb = new TCanvas("cpartialXgb","sPlot demo bdt, partial", 800, 600);  
  RooPlot* fpartialXgb = xgb->frame() ;  
  dataw_partial->plotOn(fpartialXgb, DataError(RooAbsData::SumW2) ) ;
  fpartialXgb->SetTitle("Analysis BDT");
  cpartialXgb->cd(1);
  fpartialXgb->Draw() ;
  cpartialXgb->SaveAs("sPlotXgbPartial.png");
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
  RooRealVar* KLmassD0 = ws->var("KLmassD0");

  // note, we get the dataset with sWeights
  RooDataSet* mydata = (RooDataSet*) ws->data("mydataWithSWeights");

  // create weighted data set
  RooDataSet * mydataw_sgn = new RooDataSet(mydata->GetName(),mydata->GetTitle(),mydata,*mydata->get(),0,"nsignal_sw") ;
  mydataw_sgn->Print();

  // convert to TH1
  float theinf = bdtCut-2;
  float thedelta = 15.-theinf;
  int thebin = thedelta/0.5;
  TH1 *h1_xgb    = mydataw_sgn->createHistogram("h1_xgb",*xgb,Binning(thebin,theinf,15)); 
  TH1 *h1_L1id   = mydataw_sgn->createHistogram("h1_L1id",*L1id, Binning(33,-4.,7.)); 
  TH1 *h1_L2id   = mydataw_sgn->createHistogram("h1_L2id",*L2id, Binning(33,-4.,7.)); 
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
  TH1 *h1_Kip3d      = mydataw_sgn->createHistogram("h1_Kip3d",*Kip3d, Binning(40,-0.2,0.2)); 
  TH1 *h1_KLmassD0   = mydataw_sgn->createHistogram("h1_KLmassD0",*KLmassD0, Binning(30,0.,6.));

  TFile myFileSPlots("myFileSPlots.root","RECREATE");
  myFileSPlots.cd();
  h1_xgb   -> Write();
  h1_L1id  -> Write();
  h1_L2id  -> Write();
  h1_Bprob -> Write();
  h1_BsLxy -> Write();
  h1_Bcos -> Write();
  h1_L1pt -> Write();
  h1_L2pt -> Write();
  h1_Kpt  -> Write();
  h1_LKdz -> Write();
  h1_L1L2dr -> Write();
  h1_LKdr   -> Write();
  h1_L1iso  -> Write();
  h1_Kiso   -> Write();
  h1_Passymetry -> Write();
  h1_Kip3d    -> Write();
  h1_KLmassD0 -> Write();
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

