#ifndef SPlot_JPsi
#define SPlot_JPsi

#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooDataSet.h"

class RooAbsReal;
class RooAbsPdf;
class RooFitResult;
class RooRealVar;
class RooSimultaneous;

class SPlotJPsi : public TNamed {

 public:

  SPlotJPsi();
  SPlotJPsi(const SPlotJPsi &other);
  SPlotJPsi(const char* name, const char* title);
  SPlotJPsi(const char* name, const char* title, const RooDataSet &data);
  SPlotJPsi(const char* name, const char* title,RooDataSet& data, RooAbsPdf* pdf,
	  const RooArgList &yieldsList,const RooArgSet &projDeps=RooArgSet(),
	  bool includeWeights=kTRUE, bool copyDataSet = kFALSE, const char* newName = "");

  ~SPlotJPsi();

  RooDataSet* SetSData(RooDataSet* data);
  
  RooDataSet* GetSDataSet() const;
  
  RooArgList GetSWeightVars() const;
  
  Int_t GetNumSWeightVars() const;
  
  void AddSWeight(RooAbsPdf* pdf, const RooArgList &yieldsTmp,
		  const RooArgSet &projDeps=RooArgSet(), bool includeWeights=kTRUE);
  
  Double_t GetSumOfEventSWeight(Int_t numEvent) const;
  
  Double_t GetYieldFromSWeight(const char* sVariable) const;
  
  Double_t GetSWeight(Int_t numEvent, const char* sVariable) const;
  
  
 protected:
  
  enum {
    kOwnData = BIT(20)
  };
  
  RooArgList fSWeightVars;
  
  RooDataSet* fSData;
  
  ClassDef(SPlotJPsi,1)   // Class used for making sPlots
    
  };

#endif
