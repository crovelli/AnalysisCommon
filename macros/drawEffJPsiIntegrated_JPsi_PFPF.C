#include "TString.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TVirtualFitter.h"
#include <iostream>

const TString lumiString = "";

// ----------------------------------------------------------  
// binning
const int nEtaBins = 1;
const int nPtBins = 1;
const double ptBinLimits[nPtBins+1]  = { 0., 100.};
const double ptBinCenters[nPtBins]   = { 50. };
const double ptBinHalfWidth[nPtBins] = { 50. };


// ----------------------------------------------------------
// Signal yield in data and MC: hardcoded inputs
double dataPassingFull[nEtaBins][nPtBins] = {
  { 4.7608e+03 }      
};

double dataPassingAntiD0WrtNominalBdt[nEtaBins][nPtBins] = {
  { 4.7608e+03 }      
};

double dataPassingBdtWrtAntiD0[nEtaBins][nPtBins] = {
  { 4.7608e+03 }      
};

double dataPassingAntiD0WrtBdt0[nEtaBins][nPtBins] = {
  { 9.1590e+03 }      
};

double dataPassingErrFull[nEtaBins][nPtBins] = {
  { 7.95e+01 }     
};

double dataPassingErrAntiD0WrtNominalBdt[nEtaBins][nPtBins] = {
  { 7.95e+01 }     
};

double dataPassingErrBdtWrtAntiD0[nEtaBins][nPtBins] = {
  { 7.95e+01 }     
};
double dataPassingErrAntiD0WrtBdt0[nEtaBins][nPtBins] = {
  { 351. }     
};

double dataTotalFull[nEtaBins][nPtBins] = {
  { 9.6654e+03 }      
};

double dataTotalAntiD0WrtNominalBdt[nEtaBins][nPtBins] = {
  { 5.2079e+03 }      
};

double dataTotalBdtWrtAntiD0[nEtaBins][nPtBins] = {
  { 9.1590e+03 }      
};

double dataTotalAntiD0WrtBdt0[nEtaBins][nPtBins] = {
  { 9.6654e+03 }      
};

double dataTotalErrFull[nEtaBins][nPtBins] = {
  { 254. }      
};

double dataTotalErrAntiD0WrtNominalBdt[nEtaBins][nPtBins] = {
  { 84.2 }      
};

double dataTotalErrBdtWrtAntiD0[nEtaBins][nPtBins] = {
  { 351. }      
};

double dataTotalErrAntiD0WrtBdt0[nEtaBins][nPtBins] = {
  { 254. }      
};

double mcPassingFull[nEtaBins][nPtBins] = {
  { 42021 }      
};

double mcPassingAntiD0WrtNominalBdt[nEtaBins][nPtBins] = {
  { 42021 }      
};

double mcPassingBdtWrtAntiD0[nEtaBins][nPtBins] = {
  { 42021 }      
};

double mcPassingAntiD0WrtBdt0[nEtaBins][nPtBins] = {
  { 79251 }      
};

double mcTotalFull[nEtaBins][nPtBins] = {
  { 90329 }      
};

double mcTotalAntiD0WrtNominalBdt[nEtaBins][nPtBins] = {
  { 46986 }      
};

double mcTotalBdtWrtAntiD0[nEtaBins][nPtBins] = {
  { 79251 }      
};

double mcTotalAntiD0WrtBdt0[nEtaBins][nPtBins] = {
  { 90329 }      
};




void drawResults( int whichEff){

  double dataPassing[nEtaBins][nPtBins]; 
  double dataTotal[nEtaBins][nPtBins]; 
  double dataPassingErr[nEtaBins][nPtBins]; 
  double dataTotalErr[nEtaBins][nPtBins]; 
  double mcPassing[nEtaBins][nPtBins]; 
  double mcTotal[nEtaBins][nPtBins]; 
  //
  if(whichEff==1) {
    for (int ii=0; ii<nEtaBins; ii++) {
      for (int jj=0; jj<nPtBins; jj++) {
	dataPassing[ii][jj]    = dataPassingFull[ii][jj];
	dataTotal[ii][jj]      = dataTotalFull[ii][jj];
	dataPassingErr[ii][jj] = dataPassingErrFull[ii][jj];
	dataTotalErr[ii][jj]   = dataTotalErrFull[ii][jj];
	mcPassing[ii][jj]      = mcPassingFull[ii][jj];
	mcTotal[ii][jj]        = mcTotalFull[ii][jj];
      }}
  }
  //
  if(whichEff==3) {
    for (int ii=0; ii<nEtaBins; ii++) {
      for (int jj=0; jj<nPtBins; jj++) {
	dataPassing[ii][jj]    = dataPassingAntiD0WrtNominalBdt[ii][jj];
	dataTotal[ii][jj]      = dataTotalAntiD0WrtNominalBdt[ii][jj];
	dataPassingErr[ii][jj] = dataPassingErrAntiD0WrtNominalBdt[ii][jj];
	dataTotalErr[ii][jj]   = dataTotalErrAntiD0WrtNominalBdt[ii][jj];
	mcPassing[ii][jj]      = mcPassingAntiD0WrtNominalBdt[ii][jj];
	mcTotal[ii][jj]        = mcTotalAntiD0WrtNominalBdt[ii][jj];
      }}
  }
  //
  if(whichEff==4) {
    for (int ii=0; ii<nEtaBins; ii++) {
      for (int jj=0; jj<nPtBins; jj++) {
	dataPassing[ii][jj]    = dataPassingBdtWrtAntiD0[ii][jj];
	dataTotal[ii][jj]      = dataTotalBdtWrtAntiD0[ii][jj];
	dataPassingErr[ii][jj] = dataPassingErrBdtWrtAntiD0[ii][jj];
	dataTotalErr[ii][jj]   = dataTotalErrBdtWrtAntiD0[ii][jj];
	mcPassing[ii][jj]      = mcPassingBdtWrtAntiD0[ii][jj];
	mcTotal[ii][jj]        = mcTotalBdtWrtAntiD0[ii][jj];
      }}
  }
  //
  if(whichEff==5) {
    for (int ii=0; ii<nEtaBins; ii++) {
      for (int jj=0; jj<nPtBins; jj++) {
	dataPassing[ii][jj]    = dataPassingAntiD0WrtBdt0[ii][jj];
	dataTotal[ii][jj]      = dataTotalAntiD0WrtBdt0[ii][jj];
	dataPassingErr[ii][jj] = dataPassingErrAntiD0WrtBdt0[ii][jj];
	dataTotalErr[ii][jj]   = dataTotalErrAntiD0WrtBdt0[ii][jj];
	mcPassing[ii][jj]      = mcPassingAntiD0WrtBdt0[ii][jj];
	mcTotal[ii][jj]        = mcTotalAntiD0WrtBdt0[ii][jj];
      }}
  }


  // Now data efficiency with statistical errors
  double theDataEff[nEtaBins][nPtBins];
  double theDataEffErr[nEtaBins][nPtBins];
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      theDataEff[ii][jj] = dataPassing[ii][jj]/dataTotal[ii][jj];
      theDataEffErr[ii][jj] = sqrt(theDataEff[ii][jj]*(1-theDataEff[ii][jj])/dataTotal[ii][jj]);
    }}

  // Now mc efficiency with statistical errors 
  double mcPassingErr[nEtaBins][nPtBins];
  double mcTotalErr[nEtaBins][nPtBins];
  double theMcEff[nEtaBins][nPtBins];
  double theMcEffErr[nEtaBins][nPtBins];

  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      mcPassingErr[ii][jj] = sqrt(mcPassing[ii][jj]);
      mcTotalErr[ii][jj]   = sqrt(mcTotal[ii][jj]);
      theMcEff[ii][jj]     = mcPassing[ii][jj]/mcTotal[ii][jj];
      theMcEffErr[ii][jj]  = sqrt(theMcEff[ii][jj]*(1-theMcEff[ii][jj])/mcTotal[ii][jj]);
    }}

  // -----------------------------------------
  std::cout << "Now scale factors:" << std::endl;  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;  

  double sf[nEtaBins][nPtBins];
  double sfErrTot[nEtaBins][nPtBins];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iPt=0; iPt<nPtBins; iPt++){ 
      sf[iEta][iPt] = theDataEff[iEta][iPt]/theMcEff[iEta][iPt];
      float sigmaDoD   = theDataEffErr[iEta][iPt]/theDataEff[iEta][iPt];
      float sigmaMCoMC = theMcEffErr[iEta][iPt]/theMcEff[iEta][iPt];
      sfErrTot[iEta][iPt] = sf[iEta][iPt]*sqrt( (sigmaDoD*sigmaDoD) + (sigmaMCoMC*sigmaMCoMC) );
    }
  }

  // -----------------------------------------
  std::cout << "data: passing = " << dataPassing[0][0] << ", total =" << dataTotal[0][0] << std::endl;
  std::cout << "mc: passing = "   << mcPassing[0][0]   << ", total =" << mcTotal[0][0]   << std::endl;
  std::cout << "dataEff = " << theDataEff[0][0] << " +-" << theDataEffErr[0][0] << std::endl;
  std::cout << "mcEff = " << theMcEff[0][0] << " +-" << theMcEffErr[0][0] << std::endl;
  std::cout << "SF = " << sf[0][0] << " +/- " << sfErrTot[0][0] << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << std::endl;
  std::cout << std::endl;
  
  for(int ieta = 0; ieta<nEtaBins; ieta++){

    TString cname = "sfEff";
    TCanvas *c1 = new TCanvas(cname, cname, 10,10,700,700);
    c1->SetFillColor(kWhite);
    c1->Draw();
    TPad *pad1 = new TPad("main","",0, 0.3, 1.0, 1.0);
    pad1->SetTopMargin(0.20);
    pad1->SetBottomMargin(0.02);
    pad1->SetGrid();
    TPad *pad2 = new TPad("ratio", "", 0, 0, 1.0, 0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.30);
    pad2->SetGrid();

    pad1->Draw();
    pad2->Draw();

    // Create and fill arrays for graphs for this eta bin
    double *dataSlice    = new double[nPtBins];
    double *dataSliceErr = new double[nPtBins];
    double *mcSlice      = new double[nPtBins];
    double *mcSliceErr   = new double[nPtBins];
    double *sfSlice      = new double[nPtBins];
    double *sfSliceErr   = new double[nPtBins];
    for(int ipt = 0; ipt<nPtBins; ipt++){
      dataSlice   [ipt] = theDataEff[ieta][ipt];
      dataSliceErr[ipt] = theDataEffErr[ieta][ipt];
      mcSlice     [ipt] = theMcEff[ieta][ipt];
      mcSliceErr  [ipt] = theMcEffErr[ieta][ipt];
      sfSlice     [ipt] = sf[ieta][ipt];
      sfSliceErr  [ipt] = sfErrTot[ieta][ipt];
    }

    // Create and configure the graphs   
    TGraphErrors *grData = new TGraphErrors(nPtBins, ptBinCenters, dataSlice, ptBinHalfWidth, dataSliceErr);
    TGraphErrors *grMc   = new TGraphErrors(nPtBins, ptBinCenters, mcSlice, ptBinHalfWidth, mcSliceErr);
    TGraphErrors *grSf   = new TGraphErrors(nPtBins, ptBinCenters, sfSlice, ptBinHalfWidth, sfSliceErr);
    
    grData->SetLineColor(kBlack);
    grData->SetMarkerColor(kBlack);
    grData->SetMarkerStyle(20);
    grData->SetMarkerSize(1.);

    int ci = TColor::GetColor("#99ccff");
    grMc->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grMc->SetLineColor(kGreen+4);
    grMc->SetMarkerStyle(22);
    grMc->SetMarkerColor(kGreen+4);
    grMc->SetMarkerSize(1.);

    ci = TColor::GetColor("#99ccff");
    grSf->SetFillColor(kGreen-8);
    ci = TColor::GetColor("#3399ff");
    grSf->SetLineColor(kGreen+4);
    grSf->SetMarkerStyle(20);
    grSf->SetMarkerColor(kGreen+4);
    grSf->SetMarkerSize(1.);
    // grSf->Fit("pol0","","",20,350);

    // Create and configure the dummy histograms on which to draw the graphs
    TH2F *h1 = new TH2F("dummy1","", 100, 0, 100.5, 100, 0., 1.1);
    h1->GetYaxis()->SetTitle("Efficiency");
    h1->SetStats(0);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetDecimals();
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetYaxis()->SetTitleSize(0.05);
    TH2F *h2 = new TH2F("dummy2","", 100, 0, 100.5, 100, 0.2, 1.5);
    h2->GetXaxis()->SetTitle("B p_{T} [GeV]");
    h2->GetYaxis()->SetTitle("Scale Factor");
    h2->GetXaxis()->SetTitleOffset(1.0);
    h2->GetXaxis()->SetTitleSize(0.1);
    h2->GetYaxis()->SetTitleOffset(0.4);
    h2->GetYaxis()->SetTitleSize(0.1);
    h2->GetXaxis()->SetLabelSize(0.08);
    h2->GetYaxis()->SetLabelSize(0.08);
    h2->GetYaxis()->SetNdivisions(505);
    h2->GetYaxis()->SetDecimals();
    h2->SetStats(0);

    TLegend *leg = new TLegend(0.65,0.1,0.9,0.25);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(grData, "Data", "pl");
    leg->AddEntry(grMc, "MC (cut&count)", "pFlE");

    TLatex *latLumi = new TLatex(0, 1.15, lumiString);

    //TLatex *latEta = new TLatex(60.0, 0.5, etaLimitsStringArray[ieta]);


    // --------------------------------------
    // Draw the efficiencies
    pad1->cd();
    h1->Draw();
    grMc  ->Draw("2same");
    grMc  ->Draw("pZ,same");
    grData->Draw("PEZ,same");
    leg->Draw("same");
    //latEta->Draw("same");
    latLumi->Draw("same");
    // Draw the scale factors
    pad2->cd();
    h2->Draw();
    grSf  ->Draw("2same");
    grSf  ->Draw("pEZ,same");
    // Save into a file
    TString fname = cname += ".png";
    c1->Print(fname);
  }

}



