import ROOT
from ROOT import RooFit
import math
from roofit_helper import *
#ROOT.gROOT.ProcessLine(open('roofit_models.h').read())
#from ROOT import DoubleSidedCB
#from ROOT import ROOT_DoubleSidedCB
from roofit_models import root_function_DoubleSidedCB
ROOT.gInterpreter.Declare(root_function_DoubleSidedCB)
#ROOT.gStyle.SetOptFit(0000);
ROOT.gROOT.SetBatch(True);
ROOT.gROOT.SetStyle("Plain");
msgservice = ROOT.RooMsgService.instance()
msgservice.setGlobalKillBelow(RooFit.FATAL)
import numpy as np
import csv
import os.path
import atexit

nbin_data = 50

def residuals(xframe, var,name):
   hresid = xframe.residHist()
   xframe2 = var.frame()
   xframe2.addPlotable(hresid,"P")
   c2=canvas_create(xframe2,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',False)
   c2.SaveAs(name+'_residual.png')
   hpull = xframe.pullHist()
   xframe3 = var.frame()
   xframe3.addPlotable(hpull,"P")
   c3=canvas_create(xframe3,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',False)
   c3.SaveAs(name+'_pull.png')


def define_workspace_bmass_data(wspace_name,mB_branch,tree,Bmass_min=4.7):
   wspace = ROOT.RooWorkspace(wspace_name)
   fitvars = ROOT.RooArgSet()
   bMass = ROOT.RooRealVar(mB_branch, "m(K^{+}e^{+}e^{-})", Bmass_min, 5.7, "GeV")
   fitvars.add(bMass)
   dataset = ROOT.RooDataSet('data','data',tree, ROOT.RooArgSet(fitvars))
   theBMassfunc = ROOT.RooFormulaVar("x", "x", "@0", ROOT.RooArgList(bMass) )
   theBMass     = dataset.addColumn(theBMassfunc) ;
   theBMass.setRange(Bmass_min,5.7);
   fitvars.add(theBMass)
   getattr(wspace, "import")(dataset, RooFit.Rename('data'))
   # When a RooWorkspace is set as an attribute of a class, it can trigger a memory error
   # This is solved in the latest ROOT version > 6.14.08
   atexit.register(wspace.Delete)
   return wspace,dataset,bMass,theBMass


def define_workspace_bmass_data_withweight(wspace_name,mB_branch,tree,Bmass_min=4.7):
   wspace = ROOT.RooWorkspace(wspace_name)
   fitvars = ROOT.RooArgSet()
   bMass = ROOT.RooRealVar(mB_branch, "m(K^{+}e^{+}e^{-})", Bmass_min, 5.7, "GeV")
   weight = ROOT.RooRealVar("weight", "weight", -1, 100, "")  
   fitvars.add(bMass)
   fitvars.add(weight)
   dataset = ROOT.RooDataSet('data','data',tree, ROOT.RooArgSet(fitvars), '', 'weight')
   dataset.Print("V")
   print "is weighted = ", dataset.isWeighted()
   theBMassfunc = ROOT.RooFormulaVar("x", "x", "@0", ROOT.RooArgList(bMass) )
   theBMass     = dataset.addColumn(theBMassfunc) ;
   theBMass.setRange(Bmass_min,5.7);
   fitvars.add(theBMass)
   getattr(wspace, "import")(dataset, RooFit.Rename('data'))
   # When a RooWorkspace is set as an attribute of a class, it can trigger a memory error
   # This is solved in the latest ROOT version > 6.14.08
   atexit.register(wspace.Delete)
   return wspace,dataset,bMass,theBMass

def get_visible_yield(obs, results, pdf, amplitude):
   intgral_pdf = pdf.createIntegral(obs,obs,"window")
   visible = amplitude.getVal() * intgral_pdf.getVal()
   return visible

def get_visible_yield_error(obs, results, pdf, amplitude):
   intgral_pdf = pdf.createIntegral(obs,obs,"window")
   intgral_pdf_err = intgral_pdf.getPropagatedError(results, obs)  
   visible = amplitude.getVal() * intgral_pdf.getVal()
   visible_err = visible * math.sqrt(pow(amplitude.getError()/amplitude.getVal(), 2) + pow(intgral_pdf_err/intgral_pdf.getVal(), 2)) if amplitude.getVal() != 0.0 else 0.0
   return visible, visible_err



##################### signal fit == Double sided CB ######################
def signal_fit(tree, outputfile, branches, SavePlot=True):
   print "Signal"
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_signal",branches[0],tree)
   # signal
   wspace.factory('mean[5.2873e+00, 5.25e+00, 5.3e+00]')
   wspace.factory('width[5.0642e-02, 1.0e-6, 5.0e-1]')
   wspace.factory('alpha1[8.1430e-01, 0.0, 10.0]')
   wspace.factory('n1[9.8615e+01, 0.0, 100.0]')
   wspace.factory('CBShape::sigcb(x,mean,width,alpha1,n1)')
   wspace.factory('mean2[5.2274e+00, 5.0e+00, 5.30e+00]')
   wspace.factory('width2[9.3738e-02, 1.0e-6, 5.0e-1]')
   wspace.factory('RooGaussian::sigg(x,mean2,width2)')
   wspace.factory('frac[4.4875e-01, 0, 1.0]')
   wspace.factory('SUM::sig(sigcb,frac*sigg)')

   sgnframe=theBMass.frame()
   wspace.factory('nsig[1000,0,100000000]')
   wspace.factory('RooExtendPdf::esig(sig,nsig)')
   sig=wspace.pdf('sig')
   nsig=wspace.var('nsig')
   esig=wspace.pdf('esig')
   results = esig.fitTo(dataset,RooFit.Extended(True),RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1))

   results.Print()
   dataset.plotOn(sgnframe,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   esig.plotOn(sgnframe, RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected), RooFit.LineColor(ROOT.kBlue), RooFit.LineWidth(2) )
   n_param = results.floatParsFinal().getSize()  
   print "chi2",sgnframe.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()
  
   c1=canvas_create(sgnframe,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   #c1, top, bottom =canvas_create_pull(sgnframe,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',theBMass)
   #top.cd()
   CMS_lumi()
   if SavePlot:
     c1.SaveAs('sgn_eek_'+outputfile+'.png')
     residuals(sgnframe,theBMass,'sgn_eek_'+outputfile)
   params=esig.getParameters(ROOT.RooArgSet(bMass))
   #return {"mean":1}
   return {"mean":params.getRealValue('mean'),"width": params.getRealValue('width'),"alpha1":params.getRealValue('alpha1'),"n1":params.getRealValue('n1'),"frac":params.getRealValue('frac'),"gauss_mean":params.getRealValue('mean2'),"gauss_width":params.getRealValue('width2')}

############################ BKG fit #################################

###############################  other B fit==Expo ##############################
def otherB_fit(tree, outputfile, branches, SavePlot=True):
   print "otherB"
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_otherB_bkg",branches[0],tree)
   wspace.factory('exp_alpha_otherb[-6.7, -100.0, -1.e-4]')
   exp_alpha_otherb = wspace.var('exp_alpha_otherb')
   wspace.factory('Exponential::exp_otherb(x,exp_alpha_otherb)')
   wspace.factory('notherB[1000,0,10e+6]')
   wspace.factory('RooExtendPdf::eexp_otherb(exp_otherb,notherB)')
   notherB=wspace.var('notherB')
   eexp_otherb=wspace.pdf('eexp_otherb')
   results = eexp_otherb.fitTo(dataset,RooFit.Extended(True),RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1))
   results.Print()
   otherb_frame=theBMass.frame()
   dataset.plotOn(otherb_frame,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   eexp_otherb.plotOn(otherb_frame,RooFit.Name("eexpOB"),RooFit.LineColor(30),RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected),RooFit.LineWidth(3))
   c1=canvas_create(otherb_frame,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   #c1, top, bottom =canvas_create_pull(otherb_frame,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',theBMass)
   #top.cd()
   CMS_lumi()
   if SavePlot:
     c1.SaveAs('bkg_otherB_'+outputfile+'.png')
     residuals(otherb_frame,theBMass,'bkg_otherB_'+outputfile)
   params=eexp_otherb.getParameters(ROOT.RooArgSet(bMass))
   n_param = results.floatParsFinal().getSize()
   print "chi2",otherb_frame.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()

   return {"exp_alpha_otherb":params.getRealValue('exp_alpha_otherb')}

############################ KDE fit #######################
def kde_fit(tree, outputfile, branches, pdfname, SavePlot=True, par=1.5, Bmass_min=4.7):
   print "KDE - {}".format(pdfname)
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_kde",branches[0],tree, Bmass_min=Bmass_min)
   kde_frame=theBMass.frame()
   wspace.factory('KeysPdf::{0}(x,data,MirrorLeft,2.0)'.format(pdfname))
   wspace.factory('KeysPdf::{0}(x,data,MirrorLeft,{1})'.format(pdfname, par))
   kde = wspace.pdf(pdfname)

   dataset.plotOn(kde_frame,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   kde.plotOn(kde_frame, RooFit.LineColor(ROOT.kBlue), RooFit.LineWidth(2) )

   #wf = ROOT.TFile('ws_bkg_kde_'+outputfile+'_{}.root'.format(pdfname), "RECREATE")
   #wspace.Write()
   #wf.Close()

   c1=canvas_create(kde_frame,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   #c1, top, bottom =canvas_create_pull(kde_frame,Bmass_min,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',theBMass)
   #top.cd()
   CMS_lumi()
   if SavePlot:
     c1.SaveAs('bkg_kde_'+outputfile+'_{}.png'.format(pdfname))
     residuals(kde_frame,theBMass,'bkg_kde_'+outputfile+'_'+pdfname)
   return kde


############################ Combinatorial fit =Expo ###########################
def bkg_fit(tree, outputfile, branches, SavePlot=True):
   print "combinatorial"
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_comb_bkg",branches[0],tree)
   wspace.factory('exp_alpha_comb[-1.5, -100.0, -1.e-4]')
   wspace.factory('Exponential::exp_comb(x,exp_alpha_comb)')
   wspace.factory('ncomb[1000,0,10e+6]')
   wspace.factory('RooExtendPdf::eexp_comb(exp_comb,ncomb)')
   ncomb=wspace.var('ncomb')
   eexp_comb=wspace.pdf('eexp_comb')

   results = eexp_comb.fitTo(dataset,RooFit.Extended(True),RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1) )
   results.Print()
   combframe=theBMass.frame()
   dataset.plotOn(combframe,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   eexp_comb.plotOn(combframe,RooFit.Name("eexp_comb"), RooFit.LineColor(49),RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected),RooFit.LineWidth(3))
   params=eexp_comb.getParameters(ROOT.RooArgSet(bMass))
   c1=canvas_create(combframe,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   # c1, top, bottom =canvas_create_pull(combframe,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',theBMass)
   #top.cd()
   CMS_lumi()
   if SavePlot: 
     c1.SaveAs('bkg_comb_'+outputfile+'.png')
   n_param = results.floatParsFinal().getSize()
   print "chi2",combframe.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()
   if SavePlot:
     residuals(combframe,theBMass,'bkg_comb_'+outputfile)

   return {"exp_alpha_comb":params.getRealValue('exp_alpha_comb')}



############################# total fit ##############################
def total_fit(tree, outputfile, branches, signal_parameters=None,  otherB_pdf=None, KstarJpsi_pdf=None, KstarPlusJpsi_pdf=None, comb_parameters=None,Significance_range=None, partial_ratio=None, partial_ratio_kstarplus=None, mvacut="",log='log.csv'):

   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_total",branches[0],tree)
   #wspace,dataset,bMass,theBMass = define_workspace_bmass_data_withweight("wspace_total",branches[0],tree)
   print "Total"
   #amplitudes
   wspace.factory('nsignal[10000.0, 0.0, 1000000.0]' )
   wspace.factory('ncomb[500.0, 0.0, 1000000.0]')
   wspace.factory('notherB[1000.0, 0.0, 1000000.0]')
   wspace.factory('frac_partial[0.35, 0.0, 10.0]')
   wspace.factory('prod::nKstarJpsi(frac_partial,nsignal)')
   wspace.factory('frac_kstarplus[0.1, 0.0, 10.0]')
   wspace.factory('prod::nKstarPlusJpsi(frac_kstarplus,nKstarJpsi)')

   # signal
   wspace.factory('mean[5.278e+00, 5.22e+00, 5.5e+00]')
   wspace.factory('width[5.8851e-02, 1.0e-6, 5.0e-1]')
   wspace.factory('alpha1[1.85, 0.0, 10.0]')
   wspace.factory('n1[19.9999, 0.0, 2000.0]')
   wspace.factory('CBShape::cb_signal(x,mean,width,alpha1,n1)')
   wspace.factory('gauss_mean[5.19e+00, 5.0e+00, 5.30e+00]')
   wspace.factory('gauss_width[1.3367e-01, 1.0e-6, 5.0e-1]')
   wspace.factory('RooGaussian::g_signal(x,gauss_mean,gauss_width)')
   wspace.factory('frac[0.5, 0, 1.0]')
   wspace.factory('SUM::signal(cb_signal,frac*g_signal)')

   # other B - bkg
   getattr(wspace, "import")(otherB_pdf, RooFit.Rename('otherb'))

   # K*Jpsi - bkg
   getattr(wspace, "import")(KstarJpsi_pdf, RooFit.Rename('kstarjpsi'))

   # K*+Jpsi - bkg
   getattr(wspace, "import")(KstarPlusJpsi_pdf, RooFit.Rename('kstarplusjpsi'))

   # combinatorial - bkg
   wspace.factory('exp_alpha_comb[-1.0, -10.0, -1.e-4]')
   wspace.factory('Exponential::exp_comb(x,exp_alpha_comb)')

   #sum
   wspace.factory('SUM::model(nsignal*signal,ncomb*exp_comb,nKstarJpsi*kstarjpsi,nKstarPlusJpsi*kstarplusjpsi,notherB*otherb)')
 
   model = wspace.pdf('model');   
   signal = wspace.pdf('signal');        
   exp_comb = wspace.pdf('exp_comb')
   otherb = wspace.pdf('otherb')
   kstarjpsi = wspace.pdf('kstarjpsi')
   kstarplusjpsi = wspace.pdf('kstarplusjpsi')
   nsignal = wspace.var('nsignal');      
   ncomb = wspace.var('ncomb')
   nKstarJpsi = wspace.obj('nKstarJpsi')
   nKstarPlusJpsi = wspace.obj('nKstarPlusJpsi')
   notherB = wspace.var('notherB')
   frac_partial = wspace.var('frac_partial')
   frac_kstarplus = wspace.var('frac_kstarplus')

   for par in signal_parameters.keys():
     (wspace.var(par)).setVal(signal_parameters[par])
     if par not in ['mean', 'gauss_mean']:
       (wspace.var(par)).setConstant(True)
   
   for par in comb_parameters.keys():
     (wspace.var(par)).setVal(comb_parameters[par])
     (wspace.var(par)).setConstant(True)
     
   if partial_ratio is not None:
     wspace.var('frac_partial').setVal(partial_ratio)          
     # wspace.var('frac_partial').setConstant(True)                # chiara! in questo modo si fissa il rapporto relativo

   # chiara: this is to set a limit on frac_partial
   #frac_partial_inf = partial_ratio-0.9*partial_ratio;
   #frac_partial_sup = 2*partial_ratio
   #print ("chiara: inf, sup, mean = ", frac_partial_inf, frac_partial_sup, partial_ratio)
   #if partial_ratio is not None:
   #  wspace.var('frac_partial').setVal(partial_ratio)          
   #  wspace.var('frac_partial').setRange(frac_partial_inf,frac_partial_sup)       

   if partial_ratio_kstarplus is not None:
     wspace.var('frac_kstarplus').setVal(partial_ratio_kstarplus)  
     wspace.var('frac_kstarplus').setConstant(True)

   print('frac_partial = ',frac_partial.getVal())
   print('frac_kstarplus = ',frac_kstarplus.getVal())

   results = model.fitTo(dataset, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1))
   print results.Print()
   xframe=theBMass.frame(RooFit.Title(""))
    
   norm=1.
   dataset.plotOn(xframe,RooFit.Binning(nbin_data), RooFit.Name("datas")) 

   model.plotOn(xframe,RooFit.Name("exp_comb"),RooFit.Components("exp_comb"),RooFit.Range("Full"),RooFit.DrawOption("L"),RooFit.VLines(),RooFit.FillColor(49),RooFit.LineColor(49),RooFit.LineStyle(2),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineWidth(3))
   model.plotOn(xframe,RooFit.Name("otherb"),RooFit.Components("otherb"),RooFit.Range("Full"),RooFit.FillColor(30),RooFit.LineColor(30),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineStyle(2), RooFit.LineWidth(3),RooFit.DrawOption("L"),RooFit.MoveToBack())
   model.plotOn(xframe,RooFit.Name("kstarjpsi"),RooFit.Components("kstarjpsi"),RooFit.Range("Full"),RooFit.FillColor(30),RooFit.LineColor(12),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineStyle(2), RooFit.LineWidth(3),RooFit.DrawOption("L"),RooFit.MoveToBack())
   model.plotOn(xframe,RooFit.Name("kstarplusjpsi"),RooFit.Components("kstarplusjpsi"),RooFit.Range("Full"),RooFit.FillColor(30),RooFit.LineColor(46),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineStyle(2), RooFit.LineWidth(3),RooFit.DrawOption("L"),RooFit.MoveToBack())
   model.plotOn(xframe,RooFit.Name("signal"),RooFit.Components("signal"),RooFit.Range("Full"),RooFit.DrawOption("L"),RooFit.LineStyle(2),RooFit.LineColor(ROOT.kBlue),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected), RooFit.LineWidth(3))
   model.plotOn(xframe, RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineColor(ROOT.kRed) )

   wspace.defineSet('obs', 'x')
   obs  = wspace.set('obs') 

   params=model.getParameters(ROOT.RooArgSet(bMass))
   if Significance_range==None:
     Significance_range={"min":params["mean"]-2*params["width"],"max":params["mean"]+2*params["width"]}
   theBMass.setRange("window",Significance_range["min"],Significance_range["max"])

   #theBMass.setRange("window",5.0,5.4)
   obs2= ROOT.RooRealVar("obs2","obs2",Significance_range["min"],Significance_range["max"])
   nset = ROOT.RooArgSet(obs2)
   print nsignal.getVal(), nsignal.getVal(nset)
   print Significance_range
   ncomb_visible, ncomb_visible_err = get_visible_yield_error(obs, results, exp_comb, ncomb)
   nsig_visible, nsig_visible_err = get_visible_yield_error(obs, results, signal, nsignal)
   nKstarJpsi_visible = get_visible_yield(obs, results, kstarjpsi , nKstarJpsi)
   nKstarJpsi_visible_err = 0.0
   nKstarPlusJpsi_visible = get_visible_yield(obs, results, kstarplusjpsi , nKstarPlusJpsi)
   nKstarPlusJpsi_visible_err = 0.0
   notherB_visible, notherB_visible_err = get_visible_yield_error(obs, results, otherb, notherB)
   nbkg_visible = nKstarJpsi_visible+nKstarPlusJpsi_visible+ncomb_visible+notherB_visible

   n_param = results.floatParsFinal().getSize()
   print "chi2",xframe.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()

   c1=canvas_create(xframe,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   #c1, top, bottom =canvas_create_pull(xframe,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]',theBMass)
   #top.cd()

   legend = ROOT.TLegend(0.65,0.65,0.92,0.85)
   legend.AddEntry(xframe.findObject("exp_comb"),"Combinatorial","l");
   legend.AddEntry(xframe.findObject("kstarjpsi"),"B -> J/#psiK*","l");
   legend.AddEntry(xframe.findObject("kstarplusjpsi"),"B -> J/#psiK*+","l");
   legend.AddEntry(xframe.findObject("otherb"),"Other B","l");
   legend.AddEntry(xframe.findObject("signal"),"B -> J/#psiK","l");
   legend.SetLineColor(ROOT.kWhite)
   legend.SetTextFont(42);
   legend.SetTextSize(0.04);
   legend.AddEntry(xframe.findObject("datas"),"Data","lpe");
   legend.Draw();
   pt=pt_create(mvacut,nsig_visible,nsig_visible_err,nbkg_visible)
   pt.Draw()
   CMS_lumi()
   c1.cd()
   c1.Update()
   c1.SaveAs('total_fit_'+outputfile+'.png')
   c1.SaveAs('total_fit_'+outputfile+'.root')
   print "signal",nsig_visible,"+/-", nsig_visible_err,"comb",ncomb_visible,"+/-", ncomb_visible_err,"K* J/psi", nKstarJpsi_visible, "+/-", nKstarJpsi_visible_err, "K*+ J/psi", nKstarPlusJpsi_visible, "+/-", nKstarPlusJpsi_visible_err, "otherB", notherB_visible, "+/-",notherB_visible_err
   residuals(xframe,theBMass,outputfile+"_poutana")

   saveWS = True

   if saveWS:
      # get likelihood
     ''' 
     nll = model.createNLL(dataset)
     ROOT.RooMinuit(nll).migrad()
     nll_frame = nsignal.frame(RooFit.Bins(30),RooFit.Range(7500.0,10000.0),RooFit.Title("LL and profileLL in nsignal")) ;
     pll = nll.createProfile(ROOT.RooArgSet(nsignal))
     pll.plotOn(nll_frame,RooFit.LineColor(2),RooFit.ShiftToZero())
     nll_frame.SetMinimum(0);
     nll_frame.SetMaximum(10)
     c2 = canvas_create(nll_frame,0.0,0.0,0.0,'nsig')
     CMS_lumi()
     c2.cd()
     c2.Update()
     c2.SaveAs('nll_'+outputfile+'.png')
     '''

     wspace_output = ROOT.RooWorkspace('wspace')
     getattr(wspace_output, "import")(model)
     getattr(wspace_output, "import")(dataset)
     params = model.getParameters(dataset)
     wspace_output.saveSnapshot("nominal_values",params)
     wspace_output.Print("V")
     wspace_output.writeToFile('wspace_'+outputfile+'.root')


   csv_header = ['cut', 'nsig_total', 'nsig_total_unc', 'nKstarJpsi_total', 'nKstarJpsi_total_unc', 'nKstarPlusJpsi_total', 'nKstarPlusJpsi_total_unc', 'nsig', 'nbkg', 'ncomb', 'nKstarJpsi', 'nKstarPlusJpsi', 'notherB', 'snr', 'chi2']
   df = {}
   df['cut'] = mvacut
   df['nsig_total'] = nsignal.getVal()
   df['nsig_total_unc'] = nsignal.getError()
   df['nKstarJpsi_total'] = nKstarJpsi.getVal()
   df['nKstarJpsi_total_unc'] = nKstarJpsi.getPropagatedError(results)
   df['nKstarPlusJpsi_total'] = nKstarPlusJpsi.getVal()
   df['nKstarPlusJpsi_total_unc'] = nKstarPlusJpsi.getPropagatedError(results)
   df['nsig'] = nsig_visible
   df['nbkg'] = nbkg_visible
   df['ncomb'] = ncomb_visible
   df['nKstarJpsi'] = nKstarJpsi_visible
   df['nKstarPlusJpsi'] = nKstarPlusJpsi_visible
   df['notherB'] = notherB_visible
   df['snr'] = nsig_visible / np.sqrt(nsig_visible + nbkg_visible)
   df['chi2'] = xframe.chiSquare(n_param)
   csv_outputfile = log
   file_exists = os.path.isfile(csv_outputfile)
   with open (csv_outputfile, 'a+') as filedata:                            
     writer = csv.DictWriter(filedata, delimiter=',', fieldnames=csv_header)
     if not file_exists:
       writer.writeheader()
     writer.writerow(df) 

   return (nsig_visible, ncomb_visible, nKstarJpsi_visible, notherB_visible, nsignal.getVal() )


#################################### main ################################
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Unbinned likelihood fit')
    parser.add_argument("-i", "--inputfile", dest="inputfile", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_v5.1_12B_Mu9_*_data.root", help="Input data file")
    parser.add_argument("-o", "--outputfile", dest="outputfile", default="test", help="Output name")
    parser.add_argument("--mvacut", dest="mva",default=6.2,type=float)
    parser.add_argument("--isgn", "--isgnfile", dest="isgnfile", default="mc_files/reg/BParkingNANO_2021Mar05_BuToKJpsi_Toee_v2_BToKEEAnalyzer_2021Mar12_HLTMu9IP6_mc_tighterPreselectionWODmass_mva_pf.root", help="Input signal file")
    parser.add_argument("--ibkg", dest="ibkg", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_v5.1_12B_Mu9_*_samesign.root", help="Input combinatorial BKG file")
    parser.add_argument("--iotherB", dest="iotherB_BKG", default="OtherB_BKGtree_KJpsiEE.root", help="Input combinatorial BKG file")
    parser.add_argument("--iKstarJpsi_BKG", dest="iKstarJpsi_BKG", default="mc_files/reg/BParkingNANO_2021Mar05_BdToKstarJpsi_Toee_v2_BToKEEAnalyzer_2021Mar12_HLTMu9IP6_mc_tighterPreselectionWODmass_kaon_mva_pf.root", help="Input combinatorial BKG file")
    parser.add_argument("--iKstarPlusJpsi_BKG", dest="iKstarPlusJpsi_BKG", default="mc_files/reg/BParkingNANO_2021Mar05_BdToKstarJpsi_Toee_v2_BToKEEAnalyzer_2021Mar12_HLTMu9IP6_mc_tighterPreselectionWODmass_kaon_mva_pf.root", help="Input combinatorial BKG file")
    parser.add_argument("--fit_primitive", dest="fit_primtv", default=False, action='store_true', help="primitive fits for fixing params or use defaults")
    parser.add_argument("--skip_realfit", dest="skip_realfit", default=False, action='store_true', help="does not perform the final fit. useful for defining parameters, tests on pdfs")
    parser.add_argument("--sel_primitive", dest="sel_primtv", default=None,  help="runs only the selected primitive fits. Options: sgn bkg_comb bkg_otherb bkg_kstar_kee bkg_kstar_piee. they can be combined with ',' in strings. No spaces")
    parser.add_argument("--minx", dest="minx", default=-1.0,type=float,  help="minx for integral")
    parser.add_argument("--maxx", dest="maxx", default=-1.0,type=float,  help="maxx for integral")
    parser.add_argument("--log", dest="log", default="log.csv", help="log of the fitting results")
    parser.add_argument("--partial_ratio", dest="partial_ratio", default=None, type=float, help="fixing the partially reco. yield")
    parser.add_argument("--partial_ratio_kstarplus", dest="partial_ratio_kstarplus", default=None, type=float, help="fixing the partially reco. yield")
    args = parser.parse_args()

    branches=["Bmass","Mll","xgb","KLmassD0","Mu9_IP6"]


    ### chiara - nominal
    cuts = branches[2]+">"+str(args.mva)+" && 2.9<"+branches[1]+" && "+branches[1]+"<3.2" + " && " + branches[3] + ">2.0" #+ " && " + branches[4]
    cuts_otherb = branches[2]+">"+str(args.mva)+" && 2.9<"+branches[1]+" && "+branches[1]+"<3.2" + " && " + branches[3] + ">2.0" #+ " && " + branches[0] + " < 5.2"
    # 
    if args.mva < 6.0:
      cuts_samesign = branches[2]+">"+str(args.mva)+" && 2.9<"+branches[1]+" && "+branches[1]+"<3.2"
    else:
      cuts_samesign = branches[2]+"> 6.0 "+" && 2.9<"+branches[1]+" && "+branches[1]+"<3.2" 


    ## chiara: denominator  
    #cuts = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && "+branches[2]+"> 0" 
    #cuts_otherb = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && "+branches[2]+"> 0" 
    #cuts_samesign = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && "+branches[2]+"> 0"
    
    ## chiara: failing
    #cuts = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && (" + branches[2]+"<"+str(args.mva) + " || " + branches[3] + "<2.0)"
    #cuts_otherb = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && (" + branches[2]+"<"+str(args.mva) + " || " + branches[3] + "<2.0)" 
    #if args.mva < 6.0:
    #  cuts_samesign = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && " + branches[2]+"<"+str(args.mva) 
    #else:
    #  cuts_samesign = branches[2]+"< 6.0 "+" && 2.9<"+branches[1]+" && "+branches[1]+"<3.2" 

    ## chiara: denominator BDT on top of antiDO
    #cuts = "2.9<"+branches[1]+" && "+branches[1]+"<3.2" + " && " + branches[3] + ">2.0"
    #cuts_otherb = "2.9<"+branches[1]+" && "+branches[1]+"<3.2" + " && " + branches[3] + ">2.0" 
    #cuts_samesign = "2.9<"+branches[1]+" && "+branches[1]+"<3.2"

    ## chiara: denominator antoDo on top of BDT 
    #cuts = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && " + branches[2]+">"+str(args.mva)
    #cuts_otherb = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && " + branches[2]+">"+str(args.mva)
    #cuts_samesign = "2.9<"+branches[1]+" && "+branches[1]+"<3.2 && " + branches[2]+"> 6.0"


    print "cuts = ", cuts  
    print "cuts_otherb = ", cuts_otherb  
    print "cuts_samesign = ", cuts_samesign  




    args.outputfile+="_wp"+str(args.mva)

    comb_parameters={'exp_alpha_comb':-1.50615}

    print "start"
    if args.fit_primtv:
      if args.sel_primtv!= None:
         args.sel_primtv = args.sel_primtv.split(",")
      else:
         args.sel_primtv = ["sgn","bkg_comb","bkg_otherb","bkg_kstarjpsi","bkg_kstarplusjpsi"]
      print "primitive params"
      if "sgn" in args.sel_primtv:
        tree_sgn = ROOT.TChain('mytreefit')
        tree_sgn.Add(args.isgnfile)
        tree_sgn_cut=tree_sgn.CopyTree(cuts)
        print tree_sgn_cut.GetEntries()
        signal_parameters = signal_fit(tree_sgn_cut, args.outputfile+"_sgnMC", branches)
        print "parameters SGN", signal_parameters['mean'], signal_parameters['width'],signal_parameters['alpha1'],signal_parameters['n1'],signal_parameters['frac'],signal_parameters['gauss_mean'],signal_parameters['gauss_width']
      
      if "bkg_otherb" in args.sel_primtv:
        tree_otherB = ROOT.TChain('mytreefit')
        tree_otherB.Add(args.iotherB_BKG)
        tree_otherB_cut=tree_otherB.CopyTree(cuts_otherb)
        otherB_pdf = kde_fit(tree_otherB_cut, args.outputfile+"_otherBMC", branches, 'otherb', par=2.0, Bmass_min=4.5)   
        print "finished kde fit for other B"
    
      if "bkg_comb" in args.sel_primtv:
        tree_bkg = ROOT.TChain('mytreefit')
        tree_bkg.Add(args.ibkg)
        tree_bkg_cut=tree_bkg.CopyTree(cuts_samesign)
        comb_parameters = bkg_fit(tree_bkg_cut, args.outputfile+"_SameSign", branches)
        print "parameters Combinatorial BKG", comb_parameters['exp_alpha_comb']

      if "bkg_kstarjpsi" in args.sel_primtv:
        tree_KstarJpsi = ROOT.TChain('mytreefit')
        tree_KstarJpsi.Add(args.iKstarJpsi_BKG)
        tree_KstarJpsi_cut=tree_KstarJpsi.CopyTree(cuts)
        KstarJpsi_pdf = kde_fit(tree_KstarJpsi_cut, args.outputfile+"_KstarJpsiMC", branches, 'kstarjpsi')
        print "finished kde fit for K* J/psi"

      if "bkg_kstarplusjpsi" in args.sel_primtv:
        tree_KstarPlusJpsi = ROOT.TChain('mytreefit')
        tree_KstarPlusJpsi.Add(args.iKstarPlusJpsi_BKG)
        tree_KstarPlusJpsi_cut=tree_KstarPlusJpsi.CopyTree(cuts)
        KstarPlusJpsi_pdf = kde_fit(tree_KstarPlusJpsi_cut, args.outputfile+"_KstarPlusJpsiMC", branches, 'kstarplusjpsi', par=2.0)
        print "finished kde fit for K*+ J/psi"

    else:
      print "Need to provide MC template"
 
    if not args.skip_realfit:
      tree = ROOT.TChain('mytreefit')
      tree.Add(args.inputfile)
      tree_cut=tree.CopyTree(cuts)
      if args.minx==-1: args.minx=signal_parameters["mean"]-2*signal_parameters["width"]
      if args.maxx==-1: args.maxx=signal_parameters["mean"]+2*signal_parameters["width"]

      nsig, nbkg,  nKstarJpsi, notherB, nsig_total =total_fit(tree_cut, args.outputfile, branches, signal_parameters,  otherB_pdf, KstarJpsi_pdf, KstarPlusJpsi_pdf, comb_parameters, {"min":args.minx,"max":args.maxx}, args.partial_ratio, args.partial_ratio_kstarplus, str(args.mva), args.log)
      print "sig",nsig,"comb", nbkg,"K* J/psi",  nKstarJpsi,"otherB", notherB,"all sig", nsig_total
      # combinatorial BKG parameters set but not fixed.
#      print "sigma",float(nsig)/math.sqrt(nsig+nkjpsi+nbkg),"nsig",float(nsig),"nbkg",float(nbkg),"Kjpsi leak",float(nkjpsi)



