#include "RooFit.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"

#include "Riostream.h"
#include <math.h>
#include "TMath.h"

#include "DoubleSidedCB2.h"

ClassImp(DoubleSidedCB2); ;

using namespace std;

DoubleSidedCB2::~DoubleSidedCB2()
{ } 

DoubleSidedCB2::DoubleSidedCB2(const char *name, const char *title,
			       RooAbsReal& _m, RooAbsReal& _mean, 
			       RooAbsReal& _width, 
			       RooAbsReal& _alpha1, RooAbsReal& _n1, 
			       RooAbsReal& _alpha2, RooAbsReal& _n2) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  mean("mean", "mean", this, _mean),
  width("width", "width", this, _width),
  alpha1("alpha1", "alpha1", this, _alpha1),
  n1("n1", "n1", this, _n1),
  alpha2("alpha2", "alpha2", this, _alpha2),
  n2("n2", "n2", this, _n2)
{ }


DoubleSidedCB2::DoubleSidedCB2(const DoubleSidedCB2& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), mean("mean", this, other.mean), width("width", this, other.width), alpha1("alpha1", this, other.alpha1), n1("n1", this, other.n1), alpha2("alpha2", this, other.alpha2), n2("n2", this, other.n2)
{ }


Double_t DoubleSidedCB2::evaluate() const {

  Double_t arg = m - mean;
  
  if (arg < 0.0)
    {
      Double_t t = (m-mean)/width; //t < 0
      
      Double_t absAlpha = fabs((Double_t)alpha1); 

      if (t >= -absAlpha) { //-absAlpha <= t < 0
        return exp(-0.5*t*t);
      }
      else {
        Double_t a = TMath::Power(n1/absAlpha,n1)*exp(-0.5*absAlpha*absAlpha);
        Double_t b = n1/absAlpha - absAlpha;

        return a/TMath::Power(b - t, n1); //b - t
      }
    }
  else 
    {
      Double_t t = (m-mean)/width; //t > 0

      Double_t absAlpha = fabs((Double_t)alpha2);

      if (t <= absAlpha) { //0 <= t <= absAlpha
        return exp(-0.5*t*t);
      }
      else {
        Double_t a = TMath::Power(n2/absAlpha,n2)*exp(-0.5*absAlpha*absAlpha);
        Double_t b = n2/absAlpha - absAlpha;
	
        return a/TMath::Power(b + t, n2); //b + t
      }
    }
  
}
