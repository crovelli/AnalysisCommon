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
#include "./SPlotRK.h"

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
void MakeHistos(RooWorkspace*, int, int);
void getDataSet(const char *, RooWorkspace*, float, float, float);   

void B_splot_fitnew(const char* workspacefile, const char* roottreefile, int q2bin, float bdtCut)
{
  // Define q^2 ranges
  float q2inf = 0;
  float q2sup = 1000;
  if (q2bin==0) { q2inf=1.05; q2sup=2.45; }
  if (q2bin==1) { q2inf=2.90; q2sup=3.20; }
  if (q2bin==2) { q2inf=3.55; q2sup=3.80; }
  if (q2bin==3) { q2inf=4.00; q2sup=4.80; }

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
  // MakeHistos(wspace, wantPFPF);

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
  SPlotRK* sData = new SPlotRK("sData","A SPlot", *mydatads, model, RooArgList(*nsignal, *notherB, *ncomb, *frac_partial) );
  
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
/*
void MakeHistos(RooWorkspace* ws, int wantPFPF){
  
  gStyle->SetOptStat(0);

  std::cout << std::endl;
  std::cout << "save histos" << std::endl;

  RooRealVar* theAnalysisBdtO = ws->var("theAnalysisBdtO");                
  RooRealVar* theAnalysisBdtG = ws->var("theAnalysisBdtG");                

  // note, we get the dataset with sWeights
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

  // create weighted data set
  RooDataSet * dataw_b   = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bYield_sw") ;
  RooDataSet * dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;
  
  // convert to TH1
  TH1 *h1_theAnalysisBdtO_B = dataw_b->createHistogram("h1_theAnalysisBdtO_B",*theAnalysisBdtO,Binning(30));         
  TH1 *h1_theAnalysisBdtO_bkg  = dataw_bkg->createHistogram("h1_theAnalysisBdtO_bkg",*theAnalysisBdtO,Binning(30));           
  TH1 *h1_theAnalysisBdtG_B = dataw_b->createHistogram("h1_theAnalysisBdtG_B",*theAnalysisBdtG,Binning(30));         
  TH1 *h1_theAnalysisBdtG_bkg  = dataw_bkg->createHistogram("h1_theAnalysisBdtG_bkg",*theAnalysisBdtG,Binning(30));           

  TFile myFileSPlots("myFileSPlots.root","RECREATE");
  myFileSPlots.cd();

  h1_theAnalysisBdtO_B->Write();        
  h1_theAnalysisBdtO_bkg->Write();        
  h1_theAnalysisBdtG_B->Write();        
  h1_theAnalysisBdtG_bkg->Write();        

  if (wantPFPF==1 || checkO==1) {
    TCanvas* ch1 = new TCanvas("ch1","ch1", 1);
    h1_theAnalysisBdtO_B->SetLineWidth(2);
    h1_theAnalysisBdtO_bkg->SetLineWidth(2);
    h1_theAnalysisBdtO_B->SetLineColor(6);
    h1_theAnalysisBdtO_bkg->SetLineColor(4);
    h1_theAnalysisBdtO_B->SetTitle("");
    h1_theAnalysisBdtO_bkg->SetTitle("");
    h1_theAnalysisBdtO_B->DrawNormalized("hist");
    h1_theAnalysisBdtO_bkg->DrawNormalized("samehist");
    ch1->SaveAs("BdtOH.png");
  }

  if (wantPFPF==1 || checkO==0){
    TCanvas* ch2 = new TCanvas("ch2","ch2", 1);
    h1_theAnalysisBdtG_B->SetLineWidth(2);
    h1_theAnalysisBdtG_bkg->SetLineWidth(2);
    h1_theAnalysisBdtG_B->SetLineColor(6);
    h1_theAnalysisBdtG_bkg->SetLineColor(4);
    h1_theAnalysisBdtG_B->SetTitle("");
    h1_theAnalysisBdtG_bkg->SetTitle("");
    h1_theAnalysisBdtG_B->DrawNormalized("hist");
    h1_theAnalysisBdtG_bkg->DrawNormalized("samehist");
    ch2->SaveAs("BdtGH.png");
  }
}
*/

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
  RooRealVar L2id("L2id", "L2id", -20.5, 20.5, "");      // range ok: -4, 7          
  RooRealVar L1id("L1id", "L1id", -20.5, 20.5, "");      // -4, 7      
  RooRealVar Bprob("Bprob","Bprob",-1.,1.,"");           // 0-1
  RooRealVar BsLxy("BsLxy","BsLxy",0.,500.,"");          ///// 0-100
  RooRealVar Bcos("Bcos","Bcos",-1.,1.,"");         // 0.9-1
  RooRealVar L1ptrel("L1pt/Bmass","L1pt/Bmass", 0., 500.,"");          // 0-10
  RooRealVar L2ptrel("L2pt/Bmass","L2pt/Bmass", 0., 500.,"");          // 0-10
  RooRealVar Kptrel("Kpt/Bmass","Kpt/Bmass",   0., 500.,"");           // 0-10
  RooRealVar LKdz("LKdz",  "LKdz",     0., 50.,"");                    // 0-1     
  RooRealVar L1L2dr("L1L2dr","L1L2dr", 0., 50.,"");                    // 0-3
  RooRealVar LKdr("LKdr",  "LKdr",     0., 50, "");                    // 0-4
  RooRealVar L1isorel("L1iso/L1pt","L1iso/L1pt", 0., 700.,"");         // 0-100
  RooRealVar Kisorel("Kiso/Kpt",   "Kiso/Kpt",   0., 700.,"");         // 0-50
  RooRealVar BBDphi("BBDphi", "BBDphi", -4., 4., "");                  
  RooRealVar BTrkdxy2("BTrkdxy2", "BTrkdxy2", 0., 30., "");             // 0-2  
  RooRealVar Passymetry("Passymetry", "Passymetry", -1., 1., "");
  RooRealVar Kiprel("Kip3d/Kip3dErr", "Kip3d/Kip3dErr", -50, 50, "");  // -7,7

  RooArgSet setall(Bmass,x,Mll,L2id,L1id,xgb,KLmassD0);
  setall.add(Bprob);
  setall.add(BsLxy);
  setall.add(Bcos);
  setall.add(L1ptrel);
  setall.add(L2ptrel);
  setall.add(Kptrel);
  setall.add(LKdz);
  setall.add(L1L2dr);
  setall.add(LKdr);
  setall.add(L1isorel);
  setall.add(Kisorel);
  setall.add(BBDphi);
  setall.add(BTrkdxy2);
  setall.add(Passymetry);
  setall.add(Kiprel);

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

