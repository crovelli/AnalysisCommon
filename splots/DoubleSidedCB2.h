#ifndef Double_Sided_CB2
#define Double_Sided_CB2

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsPdf;

class DoubleSidedCB2 : public RooAbsPdf {

 public:

  DoubleSidedCB2();
  
  DoubleSidedCB2(const char *name, const char *title, RooAbsReal& _m, RooAbsReal& _mean, RooAbsReal& _width, RooAbsReal& _alpha1, RooAbsReal& _n1, RooAbsReal& _alpha2, RooAbsReal& _n2);
  
  DoubleSidedCB2(const DoubleSidedCB2& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new DoubleSidedCB2(*this, newname); }

  ~DoubleSidedCB2();
  
 protected:
  
  RooRealProxy m;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;

  Double_t evaluate() const;

 private:

  ClassDef(DoubleSidedCB2,1) // Double Crystal Ball lineshape PDF
    };

#endif

