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

// JPsi bin
double dataPassingBdtWrtAntiD0JPsi[nEtaBins][nPtBins] = {
  // { 4.7608e+03 }      // PFPF
  { 2.0326e+03 }         // PFLP
};

double dataPassingErrBdtWrtAntiD0JPsi[nEtaBins][nPtBins] = {
  // { 79.5 }        // PFPF
  { 51.9 }           // PFLP
};

double dataTotalBdtWrtAntiD0JPsi[nEtaBins][nPtBins] = {
  //{ 8780.5 }        // PFPF
  { 6720.5 }          // PFLP
};

double dataTotalErrBdtWrtAntiD0JPsi[nEtaBins][nPtBins] = {
  // { 233 }          // PFPF
  { 255 }             // PFLP 
};

double mcPassingBdtWrtAntiD0JPsi[nEtaBins][nPtBins] = {
  // { 42021 }        // PFPF
  { 18390 }           // PFLP
};

double mcTotalBdtWrtAntiD0JPsi[nEtaBins][nPtBins] = {
  // { 79251 }        // PFPF 
  { 57948 }           // PFLP
};

// --------------------------------
double dataPassingBdtWrtAntiD0LowQ2[nEtaBins][nPtBins] = {
  // { 4.6070e+03 }        // PFPF
  { 2.1171e+03 }           // PFLP
};

double dataPassingErrBdtWrtAntiD0LowQ2[nEtaBins][nPtBins] = {
  // { 7.82e+01 }          // PFPF
  { 53.2 }                 // PFLP
};

double dataTotalBdtWrtAntiD0LowQ2[nEtaBins][nPtBins] = {
  // { 8.8406e+03  }       // PFPF
  { 6.7295e+03  }         // PFLP
};

double dataTotalErrBdtWrtAntiD0LowQ2[nEtaBins][nPtBins] = {
  // { 2.43e+02 }       // PFPF
  { 227 }
};

double mcPassingBdtWrtAntiD0LowQ2[nEtaBins][nPtBins] = {
  //{ 14938 }       // PFPF  
  { 7678 }          // PFLP
};

double mcTotalBdtWrtAntiD0LowQ2[nEtaBins][nPtBins] = {
  // { 28832 }      // PFPF
  { 23022 }         // PFLP
};

double mcPassingBdtWrtAntiD0JPsiScaled[nEtaBins][nPtBins] = {
  // { 40767.1 }       // PFPF
  { 18840.1 }        // PFLP
};

double mcTotalBdtWrtAntiD0JPsiScaled[nEtaBins][nPtBins] = {
  // { 79355.6 }         // PFPF
  { 57905.7 }        // PFLP
};



void relativeError(){

  double dataPassingJPsi[nEtaBins][nPtBins]; 
  double dataTotalJPsi[nEtaBins][nPtBins]; 
  double dataPassingErrJPsi[nEtaBins][nPtBins]; 
  double dataTotalErrJPsi[nEtaBins][nPtBins]; 
  double mcPassingJPsi[nEtaBins][nPtBins]; 
  double mcTotalJPsi[nEtaBins][nPtBins]; 
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      dataPassingJPsi[ii][jj]    = dataPassingBdtWrtAntiD0JPsi[ii][jj];
      dataTotalJPsi[ii][jj]      = dataTotalBdtWrtAntiD0JPsi[ii][jj];
      dataPassingErrJPsi[ii][jj] = dataPassingErrBdtWrtAntiD0JPsi[ii][jj];
      dataTotalErrJPsi[ii][jj]   = dataTotalErrBdtWrtAntiD0JPsi[ii][jj];
      mcPassingJPsi[ii][jj]      = mcPassingBdtWrtAntiD0JPsi[ii][jj];
      mcTotalJPsi[ii][jj]        = mcTotalBdtWrtAntiD0JPsi[ii][jj];
    }}
  //
  double dataPassingLowQ2[nEtaBins][nPtBins]; 
  double dataTotalLowQ2[nEtaBins][nPtBins]; 
  double dataPassingErrLowQ2[nEtaBins][nPtBins]; 
  double dataTotalErrLowQ2[nEtaBins][nPtBins]; 
  double mcPassingLowQ2[nEtaBins][nPtBins]; 
  double mcTotalLowQ2[nEtaBins][nPtBins]; 
  double mcPassingJPsiScaled[nEtaBins][nPtBins]; 
  double mcTotalJPsiScaled[nEtaBins][nPtBins]; 
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      dataPassingLowQ2[ii][jj]    = dataPassingBdtWrtAntiD0LowQ2[ii][jj];
      dataTotalLowQ2[ii][jj]      = dataTotalBdtWrtAntiD0LowQ2[ii][jj];
      dataPassingErrLowQ2[ii][jj] = dataPassingErrBdtWrtAntiD0LowQ2[ii][jj];
      dataTotalErrLowQ2[ii][jj]   = dataTotalErrBdtWrtAntiD0LowQ2[ii][jj];
      mcPassingLowQ2[ii][jj]      = mcPassingBdtWrtAntiD0LowQ2[ii][jj];
      mcTotalLowQ2[ii][jj]        = mcTotalBdtWrtAntiD0LowQ2[ii][jj];
      mcPassingJPsiScaled[ii][jj] = mcPassingBdtWrtAntiD0JPsiScaled[ii][jj];
      mcTotalJPsiScaled[ii][jj]   = mcTotalBdtWrtAntiD0JPsiScaled[ii][jj];
    }}



  // Now data efficiency with statistical errors
  double theDataEffJPsi[nEtaBins][nPtBins];
  double theDataEffErrJPsi[nEtaBins][nPtBins];
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      theDataEffJPsi[ii][jj] = dataPassingJPsi[ii][jj]/dataTotalJPsi[ii][jj];
      theDataEffErrJPsi[ii][jj] = sqrt(theDataEffJPsi[ii][jj]*(1-theDataEffJPsi[ii][jj])/dataTotalJPsi[ii][jj]);
    }}
  //
  double theDataEffLowQ2[nEtaBins][nPtBins];
  double theDataEffErrLowQ2[nEtaBins][nPtBins];
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      theDataEffLowQ2[ii][jj] = dataPassingLowQ2[ii][jj]/dataTotalLowQ2[ii][jj];
      theDataEffErrLowQ2[ii][jj] = sqrt(theDataEffLowQ2[ii][jj]*(1-theDataEffLowQ2[ii][jj])/dataTotalLowQ2[ii][jj]);
    }}

  // Now mc efficiency with statistical errors 
  double mcPassingErrJPsi[nEtaBins][nPtBins];
  double mcTotalErrJPsi[nEtaBins][nPtBins];
  double theMcEffJPsi[nEtaBins][nPtBins];
  double theMcEffErrJPsi[nEtaBins][nPtBins];

  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      mcPassingErrJPsi[ii][jj] = sqrt(mcPassingJPsi[ii][jj]);
      mcTotalErrJPsi[ii][jj]   = sqrt(mcTotalJPsi[ii][jj]);
      theMcEffJPsi[ii][jj]     = mcPassingJPsi[ii][jj]/mcTotalJPsi[ii][jj];
      theMcEffErrJPsi[ii][jj]  = sqrt(theMcEffJPsi[ii][jj]*(1-theMcEffJPsi[ii][jj])/mcTotalJPsi[ii][jj]);
    }}

  double mcPassingErrLowQ2[nEtaBins][nPtBins];
  double mcTotalErrLowQ2[nEtaBins][nPtBins];
  double theMcEffLowQ2[nEtaBins][nPtBins];
  double theMcEffErrLowQ2[nEtaBins][nPtBins];
  double theMcEffJPsiScaled[nEtaBins][nPtBins];
  double theMcEffErrJPsiScaled[nEtaBins][nPtBins];

  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      mcPassingErrLowQ2[ii][jj]     = sqrt(mcPassingLowQ2[ii][jj]);
      mcTotalErrLowQ2[ii][jj]       = sqrt(mcTotalLowQ2[ii][jj]);
      theMcEffLowQ2[ii][jj]         = mcPassingLowQ2[ii][jj]/mcTotalLowQ2[ii][jj];
      theMcEffErrLowQ2[ii][jj]      = sqrt(theMcEffLowQ2[ii][jj]*(1-theMcEffLowQ2[ii][jj])/mcTotalLowQ2[ii][jj]);
      theMcEffJPsiScaled[ii][jj]    = mcPassingJPsiScaled[ii][jj]/mcTotalJPsiScaled[ii][jj];
      theMcEffErrJPsiScaled[ii][jj] = sqrt(theMcEffJPsiScaled[ii][jj]*(1-theMcEffJPsiScaled[ii][jj])/mcTotalJPsiScaled[ii][jj]);
    }}


  // -----------------------------------------
  std::cout << "Now scale factors:" << std::endl;  
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;  

  double sfJPsi[nEtaBins][nPtBins];
  double sfErrTotJPsi[nEtaBins][nPtBins];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iPt=0; iPt<nPtBins; iPt++){ 
      sfJPsi[iEta][iPt] = theDataEffJPsi[iEta][iPt]/theMcEffJPsi[iEta][iPt];
      float sigmaDoDJPsi   = theDataEffErrJPsi[iEta][iPt]/theDataEffJPsi[iEta][iPt];
      float sigmaMCoMCJPsi = theMcEffErrJPsi[iEta][iPt]/theMcEffJPsi[iEta][iPt];
      sfErrTotJPsi[iEta][iPt] = sfJPsi[iEta][iPt]*sqrt( (sigmaDoDJPsi*sigmaDoDJPsi) + (sigmaMCoMCJPsi*sigmaMCoMCJPsi) );
    }
  }

  double sfLowQ2[nEtaBins][nPtBins];
  double sfErrTotLowQ2[nEtaBins][nPtBins];
  for (int iEta=0; iEta<nEtaBins; iEta++){ 
    for (int iPt=0; iPt<nPtBins; iPt++){ 
      sfLowQ2[iEta][iPt] = theDataEffLowQ2[iEta][iPt]/theMcEffLowQ2[iEta][iPt];
      float sigmaDoDLowQ2   = theDataEffErrLowQ2[iEta][iPt]/theDataEffLowQ2[iEta][iPt];
      float sigmaMCoMCLowQ2 = theMcEffErrLowQ2[iEta][iPt]/theMcEffLowQ2[iEta][iPt];
      sfErrTotLowQ2[iEta][iPt] = sfLowQ2[iEta][iPt]*sqrt( (sigmaDoDLowQ2*sigmaDoDLowQ2) + (sigmaMCoMCLowQ2*sigmaMCoMCLowQ2) );
    }
  }

  // -----------------------------------------
  std::cout << "================================" << std::endl;  
  std::cout << "JPsi bin" << std::endl;
  std::cout << "data: passing = " << dataPassingJPsi[0][0] << ", total =" << dataTotalJPsi[0][0] << std::endl;
  std::cout << "mc: passing = "   << mcPassingJPsi[0][0]   << ", total =" << mcTotalJPsi[0][0]   << std::endl;
  std::cout << "dataEff = " << theDataEffJPsi[0][0] << " +-" << theDataEffErrJPsi[0][0] << std::endl;
  std::cout << "mcEff = " << theMcEffJPsi[0][0] << " +-" << theMcEffErrJPsi[0][0] << std::endl;
  std::cout << "SF = " << sfJPsi[0][0] << " +/- " << sfErrTotJPsi[0][0] << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << std::endl;
  std::cout << "LowQ2 bin" << std::endl;
  std::cout << "data: passing = " << dataPassingLowQ2[0][0] << ", total =" << dataTotalLowQ2[0][0] << std::endl;
  std::cout << "mc: passing = "   << mcPassingLowQ2[0][0]   << ", total =" << mcTotalLowQ2[0][0]   << std::endl;
  std::cout << "mc jpsi scaled: passing = "   << mcPassingJPsiScaled[0][0]   << ", total =" << mcTotalJPsiScaled[0][0]   << std::endl;
  std::cout << "dataEff = " << theDataEffLowQ2[0][0] << " +-" << theDataEffErrLowQ2[0][0] << std::endl;
  std::cout << "mcEff = " << theMcEffLowQ2[0][0] << " +-" << theMcEffErrLowQ2[0][0] << std::endl;
  std::cout << "mcEff jpsi scaled = " << theMcEffJPsiScaled[0][0] << " +-" << theMcEffErrJPsiScaled[0][0] << std::endl;
  std::cout << "SF = " << sfLowQ2[0][0] << " +/- " << sfErrTotLowQ2[0][0] << std::endl;
  std::cout << "(MC eff LowQ2 - MC eff JPsi scaled)/(MC eff  LowQ2) = " << (theMcEffLowQ2[0][0]-theMcEffJPsiScaled[0][0])/theMcEffLowQ2[0][0] << std::endl;
  std::cout << "================================" << std::endl;  
  std::cout << std::endl;
  

  // SF as error
  double mcEffErrJPsi[nEtaBins][nPtBins]; 
  double mcEffErrLowQ2[nEtaBins][nPtBins]; 
  double mcEffRelErrSFLowQ2[nEtaBins][nPtBins]; 
  double mcEffAbsErrSFLowQ2[nEtaBins][nPtBins]; 
  double mcEffAbsErrWeightLowQ2[nEtaBins][nPtBins]; 
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      //mcEffErrJPsi[ii][jj]     = fabs(1-sfJPsi[ii][jj]);    // 
      //mcEffRelErrLowQ2[ii][jj] = fabs(1-sfLowQ2[ii][jj]);   // 
      /*
      mcEffErrJPsi[ii][jj]       = 0.02;       // da scale factor PFPF
      mcEffRelErrSFLowQ2[ii][jj] = 0.01;       // da scale factor PFPF
      */
      mcEffErrJPsi[ii][jj]       = 0.05;       // da scale factor PFLP
      mcEffRelErrSFLowQ2[ii][jj] = 0.06;       // da scale factor PFLP
      //
      mcEffAbsErrSFLowQ2[ii][jj]     = mcEffRelErrSFLowQ2[ii][jj]*theMcEffLowQ2[ii][jj];
      mcEffAbsErrWeightLowQ2[ii][jj] = theMcEffLowQ2[ii][jj]-theMcEffJPsiScaled[ii][jj];
      mcEffErrLowQ2[ii][jj] = (sqrt(mcEffAbsErrSFLowQ2[ii][jj]*mcEffAbsErrSFLowQ2[ii][jj]  + mcEffAbsErrWeightLowQ2[ii][jj]*mcEffAbsErrWeightLowQ2[ii][jj]))/theMcEffLowQ2[ii][jj];
    }}

  std::cout << std::endl;
  std::cout << "JPsi bin :  mcEff relative error (SF) = "     << mcEffErrJPsi[0][0] << std::endl;
  std::cout << "LowQ2 bin : mcEff relative error (SF) = "     << mcEffRelErrSFLowQ2[0][0] << std::endl;
  std::cout << "LowQ2 bin : mcEff abs error (weight) =      " << mcEffAbsErrWeightLowQ2[0][0] << std::endl;
  std::cout << "LowQ2 bin : mcEff relative error (tot) = "    << mcEffErrLowQ2[0][0] << std::endl;

  // +- 1 sigma
  double mcEffPlus1sJPsi[nEtaBins][nPtBins]; 
  double mcEffMinus1sJPsi[nEtaBins][nPtBins]; 
  double mcEffPlus1sLowQ2[nEtaBins][nPtBins]; 
  double mcEffMinus1sLowQ2[nEtaBins][nPtBins]; 
  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      mcEffPlus1sJPsi[ii][jj]   = theMcEffJPsi[ii][jj]*(1+mcEffErrJPsi[ii][jj]);
      mcEffMinus1sJPsi[ii][jj]  = theMcEffJPsi[ii][jj]*(1-mcEffErrJPsi[ii][jj]);
      mcEffPlus1sLowQ2[ii][jj]  = theMcEffLowQ2[ii][jj]*(1+mcEffErrLowQ2[ii][jj]);
      mcEffMinus1sLowQ2[ii][jj] = theMcEffLowQ2[ii][jj]*(1-mcEffErrLowQ2[ii][jj]);
    }}


  std::cout << std::endl;
  std::cout << "JPsi bin" << std::endl;
  std::cout << "mcEff: " << mcEffMinus1sJPsi[0][0] << ", " << theMcEffJPsi[0][0] << ", " << mcEffPlus1sJPsi[0][0] << endl;
  std::cout << "LowQ2 bin" << std::endl;
  std::cout << "mcEff: " << mcEffMinus1sLowQ2[0][0] << ", " << theMcEffLowQ2[0][0] << ", " << mcEffPlus1sLowQ2[0][0] << endl;
  std::cout << "================================" << std::endl;  


  // Efficiency ratio
  double effRatio[nEtaBins][nPtBins]; 
  double effRatioMinus1s[nEtaBins][nPtBins];
  double effRatioPlus1s[nEtaBins][nPtBins];
  double deltaEffRatioMinus1s[nEtaBins][nPtBins];
  double deltaEffRatioPlus1s[nEtaBins][nPtBins];
  double deltaEffRatioMax[nEtaBins][nPtBins];
  double deltaEffRatioMaxOverRatio[nEtaBins][nPtBins];

  for (int ii=0; ii<nEtaBins; ii++) {
    for (int jj=0; jj<nPtBins; jj++) {
      effRatio[ii][jj] = theMcEffLowQ2[ii][jj]/theMcEffJPsi[ii][jj];
      effRatioMinus1s[ii][jj] = mcEffMinus1sLowQ2[ii][jj]/mcEffMinus1sJPsi[ii][jj];
      effRatioPlus1s[ii][jj]  = mcEffPlus1sLowQ2[ii][jj]/mcEffPlus1sJPsi[ii][jj];
      deltaEffRatioMinus1s[ii][jj] = fabs(effRatioMinus1s[ii][jj] - effRatio[ii][jj]);
      deltaEffRatioPlus1s[ii][jj]  = fabs(effRatioPlus1s[ii][jj] - effRatio[ii][jj]);
      deltaEffRatioMax[ii][jj] = deltaEffRatioMinus1s[ii][jj];
      if (deltaEffRatioPlus1s[ii][jj]>deltaEffRatioMinus1s[ii][jj]) deltaEffRatioMax[ii][jj] = deltaEffRatioPlus1s[ii][jj];
      deltaEffRatioMaxOverRatio[ii][jj] = deltaEffRatioMax[ii][jj]/effRatio[ii][jj];
    }}

  std::cout << std::endl;
  std::cout << "Eff ratio = " << effRatio[0][0] << endl;
  std::cout << "Eff ratio +1s = " << effRatioPlus1s[0][0] << endl;
  std::cout << "Eff ratio -1s = " << effRatioMinus1s[0][0] << endl;
  std::cout << "Delta, +1s = " << deltaEffRatioPlus1s[0][0] << endl;
  std::cout << "Delta, -1s = " << deltaEffRatioMinus1s[0][0] << endl;
  std::cout << "Delta, max = " << deltaEffRatioMax[0][0] << endl;
  std::cout << "Final : delta max / eff ratio = " << deltaEffRatioMax[0][0]/effRatio[0][0] << std::endl;
  std::cout << "================================" << std::endl;  

}




