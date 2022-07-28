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
const int nPartialBins = 10;
const double ptBinLimits[nPartialBins+1]  = { 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
const double ptBinCenters[nPartialBins]   = { 5, 15, 25, 35, 45, 55, 65, 75, 85, 95 };
const double ptBinHalfWidth[nPartialBins] = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };

// ----------------------------------------------------------
// Signal yield in data and MC: hardcoded inputs
double dataPassing[nEtaBins][nPartialBins] = {
  { 4.7608e+03, 4.7608e+03, 4.7608e+03, 4.7608e+03, 4.7608e+03, 4.7607e+03, 4.7608e+03, 4.7608e+03, 4.7555e+03, 4.6956e+03 }
};

double dataPassingErr[nEtaBins][nPartialBins] = {
  { 7.96e+01, 7.96e+01, 7.96e+01, 7.96e+01, 7.96e+01, 7.96e+01, 7.96e+01, 7.96e+01, 7.93e+01, 6.90e+01 }
};

double dataTotal[nEtaBins][nPartialBins] = {   
  { 9.6654e+03, 9.6655e+03, 9.6572e+03, 9.6367e+03, 9.6155e+03, 9.5904e+03, 9.5632e+03, 9.5344e+03, 9.5016e+03, 9.4677e+03 }
};

double dataTotalErr[nEtaBins][nPartialBins] = {  
  { 2.54e+02, 2.78e+02, 2.58e+02, 2.50e+02, 2.45e+02, 2.41e+02, 2.37e+02, 2.33e+02, 2.31e+02, 2.29e+02 } 
};

double mcPassing[nEtaBins][nPartialBins] = {
  { 42021, 42021, 42021, 42021, 42021, 42021, 42021, 42021, 42021, 42021 }
};

double mcTotal[nEtaBins][nPartialBins] = {
  { 90329, 90329, 90329, 90329, 90329, 90329, 90329, 90329, 90329, 90329 }
};


void drawResults(){

  // Now data efficiency with statistical errors
  double theDataEff[nEtaBins][nPartialBins];
  double theDataEffErr[nEtaBins][nPartialBins];
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPartialBins; jj++) {
      theDataEff[ii][jj] = dataPassing[ii][jj]/dataTotal[ii][jj];
      theDataEffErr[ii][jj] = sqrt(theDataEff[ii][jj]*(1-theDataEff[ii][jj])/dataTotal[ii][jj]);
    }}

  // Now mc efficiency with statistical errors 
  double mcPassingErr[nEtaBins][nPartialBins];
  double mcTotalErr[nEtaBins][nPartialBins];
  double theMcEff[nEtaBins][nPartialBins];
  double theMcEffErr[nEtaBins][nPartialBins];

  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPartialBins; jj++) {
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

  double sf[nEtaBins][nPartialBins];
  double sfErrTot[nEtaBins][nPartialBins];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iPt=0; iPt<nPartialBins; iPt++){ 
      sf[iEta][iPt] = theDataEff[iEta][iPt]/theMcEff[iEta][iPt];
      float sigmaDoD   = theDataEffErr[iEta][iPt]/theDataEff[iEta][iPt];
      float sigmaMCoMC = theMcEffErr[iEta][iPt]/theMcEff[iEta][iPt];
      sfErrTot[iEta][iPt] = sf[iEta][iPt]*sqrt( (sigmaDoD*sigmaDoD) + (sigmaMCoMC*sigmaMCoMC) );
      std::cout << sf[iEta][iPt] << " +/- " << sfErrTot[iEta][iPt] << std::endl;
    }
  }

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
    double *dataSlice    = new double[nPartialBins];
    double *dataSliceErr = new double[nPartialBins];
    double *mcSlice      = new double[nPartialBins];
    double *mcSliceErr   = new double[nPartialBins];
    double *sfSlice      = new double[nPartialBins];
    double *sfSliceErr   = new double[nPartialBins];
    for(int ipt = 0; ipt<nPartialBins; ipt++){
      dataSlice   [ipt] = theDataEff[ieta][ipt];
      dataSliceErr[ipt] = theDataEffErr[ieta][ipt];
      mcSlice     [ipt] = theMcEff[ieta][ipt];
      mcSliceErr  [ipt] = theMcEffErr[ieta][ipt];
      sfSlice     [ipt] = sf[ieta][ipt];
      sfSliceErr  [ipt] = sfErrTot[ieta][ipt];
    }

    // Create and configure the graphs   
    TGraphErrors *grData = new TGraphErrors(nPartialBins, ptBinCenters, dataSlice, ptBinHalfWidth, dataSliceErr);
    TGraphErrors *grMc   = new TGraphErrors(nPartialBins, ptBinCenters, mcSlice, ptBinHalfWidth, mcSliceErr);
    TGraphErrors *grSf   = new TGraphErrors(nPartialBins, ptBinCenters, sfSlice, ptBinHalfWidth, sfSliceErr);
    
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
    TH2F *h1 = new TH2F("dummy1","", 100, 0, 100., 100, 0., 1.1);
    h1->GetYaxis()->SetTitle("Efficiency");
    h1->SetStats(0);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetDecimals();
    h1->GetYaxis()->SetTitleOffset(0.8);
    h1->GetYaxis()->SetTitleSize(0.05);
    TH2F *h2 = new TH2F("dummy2","", 100, 0, 100., 100, 1.0, 1.1);
    h2->GetXaxis()->SetTitle("K*JPsi / KJPsi [%]");
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


    // --------------------------------------
    // Draw the efficiencies
    pad1->cd();
    h1->Draw();
    grMc  ->Draw("2same");
    grMc  ->Draw("pZ,same");
    grData->Draw("PEZ,same");
    leg->Draw("same");
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



