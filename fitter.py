#!/usr/bin/env python
from collections import OrderedDict
import root_numpy
import ROOT
from ROOT import RooFit
import numpy as np
ROOT.gROOT.ProcessLine(open('models.cc').read())
from ROOT import DoubleCBFast
from helper import *
#ROOT.gErrorIgnoreLevel=ROOT.kError
#ROOT.RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt
import atexit

class fitter(object):
  def __init__(self):
    self.sig_low = BLIND_LOW
    self.sig_up = BLIND_UP
    self.fom_low = B_FOM_LOW
    self.fom_up = B_FOM_UP
    self.fit_low = FIT_LOW
    self.fit_up = FIT_UP
    self.nbin_data = 50
    self.fit_init = False

  def init_fit_mc(self, drawSNR=False, 
        sigName="B^{+}#rightarrow K^{+} J/#psi(#rightarrow e^{+}e^{-})", 
        fit_var='BToKEE_fit_mass', paramOn=True, **kwargs):

    init_dict = {}
    init_dict['isMC'] = True
    init_dict['drawSNR'] = drawSNR
    init_dict['mvaCut'] = None
    init_dict['sigName'] = sigName
    init_dict['fit_var'] = fit_var
    init_dict['paramOn'] = paramOn
    init_dict['fit_init'] = True
    init_dict['var'] = {}
    init_dict['pdf'] = {}
    kwargs.update(init_dict)
    self.__dict__.update(**kwargs)

  def init_fit_data(self, drawSNR=False, params={}, blinded=False, expS=0.0, 
        sigName="B^{+}#rightarrow K^{+} J/#psi(#rightarrow e^{+}e^{-})",
        fit_var='BToKEE_fit_mass', partialfit={}, mvaCut=0.0, plotBkgUnc=False, **kwargs):

    init_dict = {}
    init_dict['isMC'] = False
    init_dict['drawSNR'] = drawSNR
    init_dict['params'] = params
    init_dict['blinded'] = blinded
    init_dict['expS'] = expS
    init_dict['sigName'] = sigName
    init_dict['fit_var'] = fit_var
    init_dict['partialfit'] = partialfit
    init_dict['mvaCut'] = mvaCut
    init_dict['plotBkgUnc'] = plotBkgUnc
    init_dict['fit_init'] = True
    init_dict['var'] = {}
    init_dict['pdf'] = {}
    kwargs.update(init_dict)
    self.__dict__.update(**kwargs)

  def init_fit_kde(self, pdfname="partial", fit_var='BToKEE_fit_mass', **kwargs):
    init_dict = {}
    init_dict['isMC'] = True
    init_dict['pdfname'] = pdfname
    init_dict['fit_var'] = fit_var
    init_dict['paramOn'] = False
    init_dict['drawSNR'] = False
    init_dict['fit_init'] = True
    init_dict['var'] = {}
    init_dict['pdf'] = {}
    kwargs.update(init_dict)
    self.__dict__.update(**kwargs)

  def load_tree(self, tree):
    self.wspace = ROOT.RooWorkspace('myWorkSpace')
    thevars = ROOT.RooArgSet()
    Mass = ROOT.RooRealVar(self.fit_var, "m(K^{+}e^{+}e^{-})", self.fit_low, self.fit_up, "GeV")

    thevars.add(Mass)

    self.fulldata = ROOT.RooDataSet('fulldata', 'fulldata', tree, ROOT.RooArgSet(thevars))
    theMassfunc = ROOT.RooFormulaVar("x", "x", "@0", ROOT.RooArgList(Mass))
    self.mass = self.fulldata.addColumn(theMassfunc) 
    self.mass.setRange(self.fit_low, self.fit_up)
    thevars.add(self.mass)

    cut = ''

    self.data = self.fulldata.reduce(thevars, cut)
    getattr(self.wspace,'import')(self.data, RooFit.Rename("data"))
    # When a RooWorkspace is set as an attribute of a class, it can trigger a memory error
    # This is solved in the latest ROOT version > 6.14.08
    atexit.register(self.wspace.Delete)

  def setup_model_sig(self): 
    self.wspace.factory('nsig[5000.0, 0.0, 1000000.0]' )

    # Double-sided Crystal-ball
    self.wspace.factory('mean[{}, {}, {}]'.format((self.fit_up+self.fit_low)/2.0, self.fit_low, self.fit_up))
    self.wspace.factory('width[4.1858e-02, 1.0e-6, 5.0e-1]')
    self.wspace.factory('alpha1[1.0, 0.0, 10.0]')
    self.wspace.factory('n1[1.0, 1.0, 20.0]')
    self.wspace.factory('alpha2[1.0, 0.0, 10.0]')
    self.wspace.factory('n2[1.0, 1.0, 20.0]')
    self.wspace.factory('GenericPdf::sig("DoubleCBFast(x,mean,width,alpha1,n1,alpha2,n2)", {x,mean,width,alpha1,n1,alpha2,n2})')

    self.var['mean'] = self.wspace.var('mean')
    self.var['width'] = self.wspace.var('width')
    self.var['alpha1'] = self.wspace.var('alpha1')
    self.var['n1'] = self.wspace.var('n1')
    self.var['alpha2'] = self.wspace.var('alpha2')
    self.var['n2'] = self.wspace.var('n2')
    self.var['nsig'] = self.wspace.var('nsig')
    self.pdf['sig'] = self.wspace.pdf('sig')

  def setup_model_bkg(self):
    self.wspace.factory('nbkg[10000.0, 0.0, 1000000.0]')

    # Exponential
    self.wspace.factory('exp_alpha[-1.0, -100.0, -1.e-4]')
    self.wspace.factory('Exponential::bkg(x,exp_alpha)')

    # Partially reconstructed bkg
    for name, info in self.partialfit.items():
      wpf = ROOT.TFile(info['filename'], "READ")
      wp = wpf.Get("myPartialWorkSpace")
      partialPDF = wp.pdf(name)
      self.wspace.factory('n{}[{}, {}, {}]'.format(name, info['expected_yield'], \
          info['expected_yield'] - 4.0*np.sqrt(info['expected_yield']), \
          info['expected_yield'] + 4.0*np.sqrt(info['expected_yield'])) \
          if 'expected_yield' in info else 'n{}[10.0, 0.0, 100000.0]'.format(name))
      getattr(self.wspace, "import")(partialPDF, RooFit.Rename(name))
      self.var[name] = self.wspace.var('n'+name)

    self.var['exp_alpha'] = self.wspace.var('exp_alpha')
    self.var['nbkg'] = self.wspace.var('nbkg')
    self.pdf['bkg'] = self.wspace.pdf('bkg')
    for name in self.partialfit.keys():
      self.pdf[name] = self.wspace.pdf(name)

  def bkg_unc_propagation(self, outputfile):
    # define the set obs = (x)
    self.wspace.defineSet('obs', 'x')
    # make the set obs known to Python
    obs  = self.wspace.set('obs')

    nsig_interested_pdf = self.pdf['sig'].createIntegral(obs,obs,"fom_window") ;
    nsig_interested_pdf_err = nsig_interested_pdf.getPropagatedError(self.results, obs)
    self.nsig_interested = self.var['nsig'].getVal() * nsig_interested_pdf.getVal()
    self.nsig_interested_err = self.nsig_interested * np.sqrt(pow(self.var['nsig'].getError()/self.var['nsig'].getVal(), 2) + pow(nsig_interested_pdf_err/nsig_interested_pdf.getVal(), 2)) if self.var['nsig'].getVal() != 0.0 else 0.0
    nbkg_comb_pdf = self.pdf['bkg'].createIntegral(obs,obs,"fom_window")
    nbkg_comb_pdf_err = nbkg_comb_pdf.getPropagatedError(self.results, obs)
    nbkg_comb = self.var['nbkg'].getVal() * nbkg_comb_pdf.getVal()
    nbkg_comb_err = nbkg_comb * np.sqrt(pow(self.var['nbkg'].getError()/self.var['nbkg'].getVal(), 2) + pow(nbkg_comb_pdf_err/nbkg_comb_pdf.getVal(), 2)) if self.var['nbkg'].getVal() != 0.0 else 0.0
    self.nbkg_total = nbkg_comb
    print("*"*80)
    print("MVA Cut: {}".format(self.mvaCut))
    if not self.fitConverged:
      print("*"*20 + "NOT COVERGE" + "*"*20)
    print("Number of signals: {}".format(self.var['nsig'].getVal()))
    print("Number of signals in 3.0 sigma: {}, uncertainty: {}".format(self.nsig_interested, self.nsig_interested_err))
    print("Number of background - combinatorial: {}, uncertainty: {}".format(nbkg_comb, nbkg_comb_err))
    for name in partialfit.keys():
      nbkg_pdf_pdf = self.pdf[name].createIntegral(obs,obs,"fom_window")
      nbkg_partial = self.var[name].getVal() * nbkg_pdf_pdf.getVal()
      self.nbkg_total += nbkg_partial
      print("Number of background - {}: {}".format(name, nbkg_partial))
    
    # Calculate 1-sigma error band of the total bkg through linear error propagation
    bkgframe = self.mass.frame()
    self.data.plotOn(bkgframe, RooFit.Binning(self.nbin_data))

    nbinx = 1000
    xvar = np.linspace(self.fom_low, self.fom_up, nbinx)
    self.fit_params = self.pdf['model'].getVariables()
    ordered_fit_params = ['exp_alpha', 'nbkg'] + ['n'+name for name in self.partialfit.keys()]
    if not self.blinded:
      ordered_fit_params += ['nsig',]
    full_bkg = ['bkg',] + [name for name in self.partialfit.keys()]
    fit_params_info = OrderedDict()
    for name in ordered_fit_params:
      fit_params_info[name] = {'mean': self.fit_params.find(name).getVal(), 'error': self.fit_params.find(name).getError()}
    self.pdf['model'].plotOn(bkgframe,RooFit.Components(",".join(full_bkg)))
    model_curve = bkgframe.getCurve()
    model_cen = np.array([model_curve.interpolate(x) for x in xvar])
    bkgframe.remove(str(0),False)
    #self.results.covarianceMatrix().Print()
    #self.results.correlationMatrix().Print()
    covMatrix = root_numpy.matrix(self.results.covarianceMatrix())
    exp_event = self.pdf['model'].expectedEvents(self.fit_params)
    fa = []
    for name, info in fit_params_info.items():
      adjust_norm = info['error'] if (name in (['nsig', 'nbkg',] + ['n'+p for p in self.partialfit.keys()])) else 0.0
      self.fit_params.setRealValue(name, info['mean']+info['error'])

      self.pdf['model'].plotOn(bkgframe,RooFit.Components(",".join(full_bkg)),RooFit.Normalization(exp_event+adjust_norm, ROOT.RooAbsReal.NumEvent))
      model_curve = bkgframe.getCurve()
      fa_plus = np.array([model_curve.interpolate(x) for x in xvar])
      bkgframe.remove(str(0),False)
      self.fit_params.setRealValue(name, info['mean']-2.0*info['error'])

      self.pdf['model'].plotOn(bkgframe,RooFit.Components(",".join(full_bkg)),RooFit.Normalization(exp_event-adjust_norm, ROOT.RooAbsReal.NumEvent))
      model_curve = bkgframe.getCurve()
      fa_minus = np.array([model_curve.interpolate(x) for x in xvar])
      bkgframe.remove(str(0),False)
      if name == 'nsig':
        fa.append(np.zeros(nbinx))
      else:
        fa.append((fa_plus - fa_minus) / (2.0*info['error']))
      # reset the params matrix
      self.fit_params.setRealValue(name, info['mean'])

    fa = np.array(fa).T
    tmp = np.array([np.asarray(np.matmul(FA, covMatrix)).flatten() for FA in fa])
    bkg_unc = np.sqrt(np.array([np.dot(t, FA) for t, FA in zip(tmp, fa)]))
    self.nbkg_total_err = np.sqrt(np.trapz(bkg_unc*bkg_unc, x=xvar)) / ((self.fit_up-self.fit_low)/self.nbin_data)

    if self.plotBkgUnc:
      fig, ax = plt.subplots()
      ax.plot(xvar, model_cen, 'b-', label=r'$N_{{\rm bkg}}={0:.1f}\pm{1:.1f}$'.format(self.nbkg_total, self.nbkg_total_err))
      ax.fill_between(xvar, model_cen-bkg_unc, model_cen+bkg_unc, facecolor='red', alpha=0.5, linewidth=0.0, label=r'$1\sigma$')
      ax.set_xlabel(r'$m(K^{+}e^{+}e^{-}) [{\rm GeV}]$')
      ax.set_ylabel(r'a.u.')
      ax.set_ylim(bottom=0)
      ax.legend(loc='upper right')
      fig.savefig(outputfile.replace('.pdf','')+'_totalbkg_1sigma.pdf', bbox_inches='tight')
   
    self.SNR = self.nsig_interested/np.sqrt(self.nsig_interested + self.nbkg_total)

    print("Total number of background: {}, uncertainty: {}".format(self.nbkg_total, self.nbkg_total_err))
    print("S/sqrt(S+B): {}".format(self.SNR))
    print("*"*80)

  def plot(self, outputfile):
    ROOT.gStyle.SetOptFit(0000);
    ROOT.gROOT.SetBatch(True);
    ROOT.gROOT.SetStyle("Plain");
    ROOT.gStyle.SetGridStyle(3);
    ROOT.gStyle.SetOptStat(000000);
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.TH1.AddDirectory(False)

    #xframe = wspace.var('x').frame(RooFit.Title("PF electron"))
    xframe = self.mass.frame()

    if self.isMC:
      self.data.plotOn(xframe, RooFit.Binning(self.nbin_data), RooFit.Name("datapoint"))
      self.pdf['model'].plotOn(xframe,RooFit.Name("global"),RooFit.Range("Full"),RooFit.LineColor(2),RooFit.MoveToBack())
      if self.paramOn:
        self.pdf['model'].paramOn(xframe,RooFit.Layout(0.60,0.92,0.73))
        #self.pdf['model'].paramOn(xframe,RooFit.Layout(0.15,0.45,0.73))
        xframe.getAttText().SetTextSize(0.03)
      legend = ROOT.TLegend(0.65,0.75,0.92,0.85);
      pt = ROOT.TPaveText(0.72,0.38,0.92,0.50,"brNDC")
      #legend = ROOT.TLegend(0.15,0.75,0.42,0.85)
      #pt = ROOT.TPaveText(0.15,0.38,0.45,0.50,"brNDC")
      legend.AddEntry(xframe.findObject("global"),"Total Fit","l")

    else:
      if self.blinded:
        self.data.plotOn(xframe, RooFit.Binning(self.nbin_data), RooFit.CutRange("SB1,SB2"), RooFit.Name("datapoint"))
      else:
        self.data.plotOn(xframe, RooFit.Binning(self.nbin_data), RooFit.Name("datapoint"))
      self.pdf['model'].plotOn(xframe,RooFit.Name("global"),RooFit.Range("Full"),RooFit.LineColor(2),RooFit.MoveToBack()) 
      self.pdf['model'].plotOn(xframe,RooFit.Name("bkg"),RooFit.Components("bkg"),RooFit.Range("Full"),RooFit.DrawOption("F"),RooFit.VLines(),RooFit.FillColor(42),RooFit.LineColor(42),RooFit.LineWidth(1),RooFit.MoveToBack())
      plotted_partial = []
      for name, info in self.partialfit.items():
        self.pdf['model'].plotOn(xframe,RooFit.Name(name),RooFit.Components("bkg,"+",".join(plotted_partial)+",{}".format(name)),RooFit.Range("Full"),RooFit.DrawOption("F"),RooFit.VLines(),RooFit.FillColor(info['color']),RooFit.LineColor(info['color']),RooFit.LineWidth(1),RooFit.MoveToBack())
        plotted_partial.append(name)
      self.pdf['model'].plotOn(xframe,RooFit.Name("sig"),RooFit.Components("sig"),RooFit.Range("Full"),RooFit.DrawOption("L"),RooFit.LineStyle(2),RooFit.LineColor(1)) 
      legend = ROOT.TLegend(0.56,0.65,0.92,0.85) #if prefix == 'BToKEE' else ROOT.TLegend(0.46,0.70,0.92,0.85)
      legend.AddEntry(xframe.findObject("bkg"),"Combinatorial","f")
      for name, info in self.partialfit.items():
        legend.AddEntry(xframe.findObject(name),info['label'],"f")
      legend.AddEntry(xframe.findObject("sig"),self.sigName,"l")

    xframe.GetYaxis().SetTitleOffset(0.9)
    xframe.GetYaxis().SetTitleFont(42)
    xframe.GetYaxis().SetTitleSize(0.05)
    xframe.GetYaxis().SetLabelSize(0.04)
    xframe.GetYaxis().SetLabelFont(42)
    xframe.GetXaxis().SetTitleOffset(0.9)
    xframe.GetXaxis().SetTitleFont(42)
    xframe.GetXaxis().SetTitleSize(0.05)
    xframe.GetXaxis().SetLabelSize(0.04)
    xframe.GetXaxis().SetLabelFont(42)

    xframe.GetYaxis().SetTitle("Events / {0:.0f} MeV".format((self.fit_up - self.fit_low)/self.nbin_data*1000.))
    xtitle = "m(K^{+}e^{+}e^{-}) [GeV]" #if prefix == 'BToKEE' else "m(K^{+}K^{-}e^{+}e^{-}) [GeV]"
    xframe.GetXaxis().SetTitle(xtitle)
    xframe.SetStats(0)
    xframe.SetMinimum(0)


    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.AddEntry(xframe.findObject("datapoint"),"Data","lpe")

    if self.drawSNR:
      pt = ROOT.TPaveText(0.7,0.35,0.92,0.63,"brNDC")
      #pt = ROOT.TPaveText(0.72,0.30,0.92,0.63,"brNDC")
      pt.SetFillColor(0)
      pt.SetBorderSize(1)
      pt.SetTextFont(42);
      pt.SetTextSize(0.04);
      pt.SetTextAlign(12)
      if self.mvaCut:
        pt.AddText("MVA cut: {0:.2f}".format(self.mvaCut))
      pt.AddText("S_{{total}}: {0:.0f}#pm{1:.0f}".format(self.nsig_total, self.nsig_total_err))
      pt.AddText("S: {0:.0f}#pm{1:.0f}".format(self.nsig_interested, self.nsig_interested_err))
      if not self.isMC:
        pt.AddText("B: {0:.0f}#pm{1:.0f}".format(self.nbkg_total, self.nbkg_total_err))
        pt.AddText("S/#sqrt{{S+B}}: {0:.1f}".format(self.SNR))
        #pt.AddText("Punzi: {0:.1f}".format(Punzi(nbkgWindow, 2.0, 5.0)))
      if not self.fitConverged:
        pt.AddText("Fit is not converged")


    # Plot results of fit on a different frame
    c2 = ROOT.TCanvas('canvas', 'canvas', 800, 600)
    c2.SetGrid()
    c2.cd()
    ROOT.gPad.SetLeftMargin(0.10)
    ROOT.gPad.SetRightMargin(0.05)
    xframe.Draw()
    legend.Draw()
    if self.drawSNR: pt.Draw()
    CMS_lumi(self.isMC)
    c2.cd()
    c2.Update()
    c2.SaveAs(outputfile.replace('.pdf','')+'.pdf')
    print("="*80)


  def fit(self, tree, outputfile):
    if not self.fit_init: return
    msgservice = ROOT.RooMsgService.instance()
    msgservice.setGlobalKillBelow(RooFit.FATAL)
    self.load_tree(tree)
    self.setup_model_sig()
    if self.isMC:
      self.wspace.factory('ExtendPdf::model(sig,nsig)')
    else:
      self.setup_model_bkg()
      self.wspace.factory('SUM::model(nsig*sig,nbkg*bkg{})'.format(','+','.join(['n'+name+'*'+name for name in self.partialfit.keys()])))
 
      self.var['mean'].setVal(self.params['mean']); self.var['mean'].setConstant(True)
      self.var['width'].setVal(self.params['width']); self.var['width'].setConstant(True)
      self.var['alpha1'].setVal(self.params['alpha1']); self.var['alpha1'].setConstant(True)
      self.var['n1'].setVal(self.params['n1']); self.var['n1'].setConstant(True)
      self.var['alpha2'].setVal(self.params['alpha2']); self.var['alpha2'].setConstant(True)
      self.var['n2'].setVal(self.params['n2']); self.var['n2'].setConstant(True)
      if self.blinded:
        self.var['nsig'].setVal(self.expS); self.var['nsig'].setConstant(True)

    self.pdf['model'] = self.wspace.pdf('model')

    self.mass.setRange("window",self.sig_low,self.sig_up) 
    self.mass.setRange("fom_window",self.fom_low,self.fom_up) 
    self.mass.setRange("SB1",self.fit_low,self.sig_low) 
    self.mass.setRange("SB2",self.sig_up,self.fit_up) 

    ## fit the model to the data.
    print('Fitting data...')
    if (not self.isMC) and self.blinded:
      self.results = self.pdf['model'].fitTo(self.data, RooFit.Extended(True), RooFit.Save(), RooFit.Range("SB1,SB2"), RooFit.SplitRange(True), RooFit.PrintLevel(-1))
    else:
      self.results = self.pdf['model'].fitTo(self.data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(self.fit_low, self.fit_up), RooFit.PrintLevel(-1))

    self.results.Print()
    self.fitConverged = True if self.results.status() == 0 else False
    self.nsig_total = self.var['nsig'].getVal()
    self.nsig_total_err = self.var['nsig'].getError()

    if not self.isMC:
      self.bkg_unc_propagation(outputfile)
    self.plot(outputfile)
    
    if self.isMC:
      return 0.0
    else:
      output = {}
      output['Stot'] = self.nsig_total
      output['StotErr'] = self.nsig_total_err
      output['S'] = self.nsig_interested
      output['SErr'] = self.nsig_interested_err
      output['B'] = self.nbkg_total
      output['BErr'] = self.nbkg_total_err
      output['exp_alpha'] = self.fit_params.find('exp_alpha').getVal()
      output['fitConverged'] = self.fitConverged
      return output


  def fit_kde(self, tree, outputfile):
    if not self.fit_init: return
    msgservice = ROOT.RooMsgService.instance()
    msgservice.setGlobalKillBelow(RooFit.FATAL)
    self.load_tree(tree)
    print('Fitting KDE...')
    self.wspace.factory('KeysPdf::{0}(x,data,MirrorLeft,2.0)'.format(self.pdfname))
    self.pdf['model'] = self.wspace.pdf(self.pdfname)
    self.plot(outputfile)

    wf = ROOT.TFile(outputfile.replace('.pdf','').replace('.root','')+'.root', "RECREATE")
    self.wspace.Write()
    wf.Close()
    print("Created a RooWorkspace file - {}, with a KeysPdf named - {}".format(outputfile.replace('.pdf','').replace('.root','')+'.root', self.pdfname))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Unbinned likelihood fit')
    parser.add_argument("-i", "--inputfile", dest="inputfile", default="", help="Input file")
    parser.add_argument("-o", "--outputfile", dest="outputfile", default="", help="Output file")
    parser.add_argument("-p", "--partial", dest="partial", action="store_true", help="Fit partially reconstructed background")
    parser.add_argument("-n", "--pdfname", dest="pdfname", default="partial", help="PDF name of the Partially reconstructed background")
    args = parser.parse_args()
    
    params = params_jpsi_pf
    partialfit = OrderedDict()
    partialfit['partial'] = {'filename': 'part_workspace_jpsi_pf.root', 'label': 'Partially Reco.', 'color': 40}
    #partialfit['partial'] = {'filename': 'part_workspace_nonresonant_lowq2_pf.root', 'label': 'Partially Reco.', 'color': 40}
    #partialfit['jpsi'] = {'filename': 'jpsi_workspace_lowq2_pf.root', 'label': 'B^{+}#rightarrow K^{+} J/#psi(#rightarrow e^{+}e^{-})', 'color': 46, 'expected_yield': 150}

    tree = ROOT.TChain('tree')
    tree.AddFile(args.inputfile)
    b_fitter = fitter()
    if args.partial:
      b_fitter.init_fit_kde(pdfname=args.pdfname)
      b_fitter.fit_kde(tree, args.outputfile)
    else:
      b_fitter.init_fit_mc()
      #b_fitter.init_fit_data(params=params, partialfit=partialfit, drawSNR=True)
      b_fitter.fit(tree, args.outputfile)


