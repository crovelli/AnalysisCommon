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

import matplotlib as mpl
mpl.use('agg')
import matplotlib.font_manager
from matplotlib import pyplot as plt
from matplotlib import rc
#.Allow for using TeX mode in matplotlib Figures
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]

ratio=5.0/7.0
fig_width_pt = 3*246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = ratio if ratio != 0.0 else (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

params = {'text.usetex' : True,
        'axes.labelsize': 24,
        'font.size': 24,
        'legend.fontsize': 20,
        'xtick.labelsize': 24,
        'ytick.labelsize': 24,
        'font.family' : 'lmodern',
        'text.latex.unicode': True,
        'axes.grid' : True,
        'text.usetex': True,
        'figure.figsize': fig_size}
plt.rcParams.update(params)

nbin_data = 20

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

def define_workspace_bmass_data(wspace_name,mB_branch,tree):
   wspace = ROOT.RooWorkspace(wspace_name)
   fitvars = ROOT.RooArgSet()
   bMass = ROOT.RooRealVar(mB_branch, "m(K^{+}e^{+}e^{-})", 4.7, 5.7, "GeV")
   fitvars.add(bMass)
   dataset = ROOT.RooDataSet('data','data',tree, ROOT.RooArgSet(fitvars))
   theBMassfunc = ROOT.RooFormulaVar("x", "x", "@0", ROOT.RooArgList(bMass) )
   theBMass     = dataset.addColumn(theBMassfunc) ;
   theBMass.setRange(4.7,5.7);
   fitvars.add(theBMass)
   getattr(wspace, "import")(dataset, RooFit.Rename('data'))
   # When a RooWorkspace is set as an attribute of a class, it can trigger a memory error
   # This is solved in the latest ROOT version > 6.14.08
   atexit.register(wspace.Delete)
   return wspace,dataset,bMass,theBMass

def get_visible_yield_error(obs, results, pdf, amplitude):
   intgral_pdf = pdf.createIntegral(obs,obs,"window")
   intgral_pdf_err = intgral_pdf.getPropagatedError(results, obs)  
   visible = amplitude.getVal() * intgral_pdf.getVal()
   if intgral_pdf.getVal()==0:
     visible_err =0
   else:
     visible_err = visible * math.sqrt(pow(amplitude.getError()/amplitude.getVal(), 2) + pow(intgral_pdf_err/intgral_pdf.getVal(), 2)) if amplitude.getVal() != 0.0 else 0.0
   return visible, visible_err

def get_visible_yield(obs, pdf, amplitude):
   intgral_pdf = pdf.createIntegral(obs,obs,"window")
   visible = amplitude * intgral_pdf.getVal()
   return visible


##################### signal fit == Double sided CB ######################
def signal_fit(tree, outputfile, branches):
   print "Signal"
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_signal",branches[0],tree)
   # signal
   wspace.factory('mean[5.272e+00, 5.22e+00, 5.5e+00]')
   wspace.factory('width[4.1858e-02, 1.0e-6, 5.0e-1]')
   wspace.factory('alpha1[1.0, 0.0, 10.0]')
   wspace.factory('n1[1.0, 1.0, 20.0]')
   wspace.factory('alpha2[1.0, 0.0, 10.0]')
   wspace.factory('n2[1.0, 1.0, 20.0]')
   wspace.factory('GenericPdf::sig( "DoubleSidedCB2(x,mean,width,alpha1,n1,alpha2,n2)",{x,mean,width,alpha1,n1,alpha2,n2})')
   
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
   CMS_lumi()
   c1.SaveAs('sgn_eek_'+outputfile+'.png')
   residuals(sgnframe,theBMass,'sgn_eek_'+outputfile)
   params=esig.getParameters(ROOT.RooArgSet(bMass))
   return {"mean":params.getRealValue('mean'),"width": params.getRealValue('width'),"alpha1":params.getRealValue('alpha1'),"n1":params.getRealValue('n1'),"alpha2":params.getRealValue('alpha2'),"n2":params.getRealValue('n2')}


###############################  B->KJpsi fit ##############################
def kjpsi_fit(tree, outputfile, branches):
   print "kjpsi"
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_kjpsi_bkg",branches[0],tree)
   wspace.factory('mean_kjpsi[4.7, 1.0, 5.0]')
   wspace.factory('width_kjpsi[0.1, 0.001, 5.0]')
   wspace.factory("RooGaussian::kjpsi(x,mean_kjpsi,width_kjpsi)")
   wspace.factory('nkjpsi[1000,0,10e+6]')
   wspace.factory('RooExtendPdf::ekjpsi(kjpsi,nkjpsi)')
   nkjpsi=wspace.var('nkjpsi')
   ekjpsi=wspace.pdf('ekjpsi')
   results = ekjpsi.fitTo(dataset,RooFit.Extended(True),RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1))
   results.Print()
   kjpframe=theBMass.frame()
   dataset.plotOn(kjpframe,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   ekjpsi.plotOn(kjpframe,RooFit.Name("ekjpsi"),RooFit.LineColor(30),RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected),RooFit.LineWidth(3))
   c1=canvas_create(kjpframe,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   CMS_lumi()
   c1.SaveAs('bkg_jpsik_'+outputfile+'.png')
   residuals(kjpframe,theBMass,'bkg_jpsik_'+outputfile)
   params=ekjpsi.getParameters(ROOT.RooArgSet(bMass))
   n_param = results.floatParsFinal().getSize()
   print "chi2",kjpframe.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()

   return {"mean_kjpsi":params.getRealValue('mean_kjpsi'),"width_kjpsi":params.getRealValue('width_kjpsi')}


############################ BKG fit #################################
def bkg_fit(tree, outputfile, branches):
   print "combinatorial"
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_comb_bkg",branches[0],tree)
   wspace.factory('exp_alpha[-1.0, -100.0, -1.e-4]')
   wspace.factory('Exponential::bkg(x,exp_alpha)')
   wspace.factory('nbkg[1000,0,10e+6]')
   wspace.factory('RooExtendPdf::ebkg(bkg,nbkg)')
   nbkg=wspace.var('nbkg')
   ebkg=wspace.pdf('ebkg')

   results = ebkg.fitTo(dataset,RooFit.Extended(True),RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1) )
   results.Print()
   bkgframe=theBMass.frame()
   dataset.plotOn(bkgframe,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   ebkg.plotOn(bkgframe,RooFit.Name("ebkg"), RooFit.LineColor(49),RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected),RooFit.LineWidth(3))
   params=ebkg.getParameters(ROOT.RooArgSet(bMass))
   c1=canvas_create(bkgframe,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV]')
   CMS_lumi()
   c1.SaveAs('bkg_comb_'+outputfile+'.png')
   n_param = results.floatParsFinal().getSize()
   print "chi2",bkgframe.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()
   residuals(bkgframe,theBMass,'bkg_comb_'+outputfile)

   return {"exp_alpha":params.getRealValue('exp_alpha')}


############################ KDE fit #######################
def kde_fit(tree, outputfile, branches, pdfname, SavePlot=True):
   print "KDE - {}".format(pdfname)
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_kde",branches[0],tree)
   kde_frame=theBMass.frame()
   wspace.factory('KeysPdf::{0}(x,data,MirrorLeft,2.0)'.format(pdfname))
   kde = wspace.pdf(pdfname)

   dataset.plotOn(kde_frame,RooFit.Binning(nbin_data), RooFit.Name("datas"))
   kde.plotOn(kde_frame, RooFit.LineColor(ROOT.kBlue), RooFit.LineWidth(2) )

   #wf = ROOT.TFile('ws_bkg_kde_'+outputfile+'_{}.root'.format(pdfname), "RECREATE")
   #wspace.Write()
   #wf.Close()

   c1=canvas_create(kde_frame,4.7,5.7,nbin_data,'m (e^{+}e^{-}K) [GeV] ')
   CMS_lumi()
   if SavePlot:
     c1.SaveAs('bkg_kde_'+outputfile+'_{}.png'.format(pdfname))
     residuals(kde_frame,theBMass,'bkg_kde_'+outputfile+'_'+pdfname)
   return kde


############################# total fit ##############################
def total_fit(tree, outputfile, branches, sgn_parameters=None, kjpsi_pdf=None, kstaree_pdf=None, bkg_parameters=None, set_sgn_yield=None,Blind_range={"min":4.7,"max":5.7} , number_of_mctoys=None,mva=None,log='log.csv'):
   wspace,dataset,bMass,theBMass = define_workspace_bmass_data("wspace_total",branches[0],tree)
   print "Total"
   #amplitudes
   wspace.factory('nsig[100.0, 0.0, 1000000.0]' )
   wspace.factory('nbkg[10000.0, 0.0, 100000000.0]')
   wspace.factory('nkjpsi[10000.0, 0.0, 1000000.0]')
   wspace.factory('nkstaree[10000.0, 0.0, 1000000.0]')

   # signal
   wspace.factory('mean[5.272e+00, 5.22e+00, 5.5e+00]')
   wspace.factory('width[4.1858e-02, 1.0e-6, 5.0e-1]')
   wspace.factory('alpha1[1.0, 0.0, 10.0]')
   wspace.factory('n1[1.0, 1.0, 20.0]')
   wspace.factory('alpha2[1.0, 0.0, 10.0]')
   wspace.factory('n2[1.0, 1.0, 20.0]')
   wspace.factory('GenericPdf::sig( "DoubleSidedCB2(x,mean,width,alpha1,n1,alpha2,n2)",{x,mean,width,alpha1,n1,alpha2,n2})')

   # Exponential - bkg
   wspace.factory('exp_alpha[-1.0, -100.0, -1.e-4]')
   alpha = wspace.var('alpha')
   wspace.factory('Exponential::bkg(x,exp_alpha)')
   
   # Gaussian - bkg 
   #wspace.factory('mean_kjpsi[4.7, 1.0, 5.0]')
   #wspace.factory('width_kjpsi[0.1, 0.001, 5.0]')
   #wspace.factory("RooGaussian::kjpsi(x,mean_kjpsi,width_kjpsi)")
   
   # KJpsi - bkg
   getattr(wspace, "import")(kjpsi_pdf, RooFit.Rename('kjpsi'))
   
   # K* ee - bkg
   getattr(wspace, "import")(kstaree_pdf, RooFit.Rename('kstaree'))
   
   #sum
   #wspace.factory('SUM::model(nsig*sig,nbkg*bkg,nkjpsi*kjpsi,nkstaree*kstaree)')
   wspace.factory('SUM::model(nsig*sig,nbkg*bkg,nkjpsi*kjpsi)')
   
   model = wspace.pdf('model');    bkg = wspace.pdf('bkg')
   sig = wspace.pdf('sig');        kjpsi = wspace.pdf('kjpsi');    
   nsig = wspace.var('nsig');      nbkg = wspace.var('nbkg')
   nkjpsi = wspace.var('nkjpsi')
   #mean = wspace.var('mean')   
   nkstaree = wspace.var('nkstaree')
   kstaree = wspace.pdf('kstaree')
   
   if set_sgn_yield!=None:
      nsig.setVal(set_sgn_yield)
      nsig.setConstant(True)
   for par in sgn_parameters.keys():
      (wspace.var(par)).setVal(sgn_parameters[par])
      (wspace.var(par)).setConstant(True)
   #for par in kjpsi_parameters.keys():
   #  (wspace.var(par)).setVal(kjpsi_parameters[par])
   #  (wspace.var(par)).setConstant(True)
   for par in bkg_parameters.keys():
      (wspace.var(par)).setVal(bkg_parameters[par]) 

   results = model.fitTo(dataset, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1))
   print results.Print()
   xframe=theBMass.frame(RooFit.Title(""))


   if Blind_range["min"]>4.7 and Blind_range["max"]<5.7:
      norm = dataset.reduce('(({0} > {1}) & ({0} < {2})) | (({0}> {3}) & ({0} < {4}))'.format(branches[0],"4.7", str(Blind_range["min"]),str(Blind_range["max"]), "5.7")).sumEntries() / dataset.reduce('({0} > {1}) & ({0} < {2})'.format(branches[0],"4.7", "5.7")).sumEntries()
      theBMass.setRange("left",4.7,Blind_range["min"])
      theBMass.setRange("right",Blind_range["max"],5.7)
      norm=1.0
      dataset.plotOn(xframe,RooFit.Binning(nbin_data), RooFit.Name("datas"),RooFit.CutRange("left,right"))
   else:
      norm=1.
      dataset.plotOn(xframe,RooFit.Binning(nbin_data), RooFit.Name("datas")) 
      

   #norm=1.
   #dataset.plotOn(xframe,RooFit.Binning(nbin_data), RooFit.Name("datas"))

   model.plotOn(xframe,RooFit.Name("bkg"),RooFit.Components("bkg"),RooFit.Range("Full"),RooFit.DrawOption("L"),RooFit.VLines(),RooFit.FillColor(49),RooFit.LineColor(49),RooFit.LineStyle(2),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineWidth(3))
   model.plotOn(xframe,RooFit.Name("kjpsi"),RooFit.Components("kjpsi"),RooFit.Range("Full"),RooFit.FillColor(30),RooFit.LineColor(30),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineStyle(2), RooFit.LineWidth(3),RooFit.DrawOption("L"),RooFit.MoveToBack())
   #model.plotOn(xframe,RooFit.Name("kstaree"),RooFit.Components("kstaree"),RooFit.Range("Full"),RooFit.FillColor(30),RooFit.LineColor(12),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineStyle(2), RooFit.LineWidth(3),RooFit.DrawOption("L"),RooFit.MoveToBack())
   model.plotOn(xframe,RooFit.Name("sig"),RooFit.Components("sig"),RooFit.Range("Full"),RooFit.DrawOption("L"),RooFit.LineStyle(2),RooFit.LineColor(ROOT.kBlue),RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected), RooFit.LineWidth(3))
   model.plotOn(xframe, RooFit.Normalization(norm, ROOT.RooAbsReal.RelativeExpected),RooFit.LineColor(ROOT.kRed) )
   
   wspace.defineSet('obs', 'x')
   obs  = wspace.set('obs')  
   theBMass.setRange("window",Blind_range["min"],Blind_range["max"])
   
   #theBMass.setRange("window",5.0,5.4)
   obs2= ROOT.RooRealVar("obs2","obs2",Blind_range["min"],Blind_range["max"])
   nset = ROOT.RooArgSet(obs2)
   print sig.getVal(), sig.getVal(nset)
   print Blind_range
   print nbkg.getVal(),nkjpsi.getVal(),nsig.getVal()
   nbkg_visible, nbkg_visible_err = get_visible_yield_error(obs, results, bkg, nbkg)
   nsig_visible, nsig_visible_err = get_visible_yield_error(obs, results, sig, nsig)
   nkjpsi_visible, nkjpsi_visible_err = get_visible_yield_error(obs, results, kjpsi, nkjpsi)
   print "hereee",nbkg_visible,nsig_visible,nkjpsi_visible
   #   c1=canvas_create(xframe,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]')
   n_param = results.floatParsFinal().getSize()
   print "chi2",xframe.chiSquare(n_param),"ndof",n_param
   print "edm",results.edm(),"log",results.minNll()
   c1=canvas_create(xframe,4.7,5.7,nbin_data,'m(e^{+}e^{-}K) [GeV]')
   
   legend = ROOT.TLegend(0.65,0.65,0.92,0.85)
   legend.AddEntry(xframe.findObject("bkg"),"Combinatorial","l");
   legend.AddEntry(xframe.findObject("kjpsi"),"B -> J/#psiK","l");
   #legend.AddEntry(xframe.findObject("kstaree"),"B -> eeK*","l");
   legend.AddEntry(xframe.findObject("sig"),"B -> eeK","l");
   legend.SetLineColor(ROOT.kWhite)
   legend.SetTextFont(42);
   legend.SetTextSize(0.04);
   legend.AddEntry(xframe.findObject("datas"),"Data","lpe");
   legend.Draw();
   pt=pt_create(mva,nsig_visible,nsig_visible_err,nkjpsi_visible+nbkg_visible)
   pt.Draw()
   CMS_lumi()
   c1.cd()
   c1.Update()
   c1.SaveAs('total_fit_'+outputfile+'.png')
   print nsig_visible, nsig_visible_err, nbkg_visible_err, nkjpsi_visible_err
   residuals(xframe,theBMass,"poutana")
   if number_of_mctoys!=None:
      nsig.setConstant(False)
      mctoys = ROOT.RooMCStudy(model, ROOT.RooArgSet(theBMass),
               RooFit.Binned( ROOT.kTRUE),
               #RooFit.Binned( ROOT.kFALSE),
                  ROOT.RooFit.Silence(),
                  RooFit.Extended(),
                  RooFit.FitOptions(
                     RooFit.Save(), RooFit.Range(4.7,5.7), RooFit.PrintLevel(-1)
                  )
             )
      mctoys.generateAndFit(number_of_mctoys)

      frame1 = mctoys.plotParam(nsig, ROOT.RooFit.Bins(40))
      frame2 = mctoys.plotError(nsig, ROOT.RooFit.Bins(40))
      frame3 = mctoys.plotPull(nsig, ROOT.RooFit.Bins(40), 
                                      ROOT.RooFit.FitGauss(ROOT.kTRUE)
                               )
      # Plot distribution of minimized likelihood
      frame4 = mctoys.plotNLL(ROOT.RooFit.Bins(40))
      cpr=canvas_create(frame1,4.7,5.7,1,'Distribution of the fitted value of N_{sgn}',False)
      cpr.SaveAs("cpr_"+outputfile+".png")
      cerr=canvas_create(frame2,0,1,1,'Distribution of the fitted error of N_{sgn}',False)
      cerr.SaveAs("cerr_"+outputfile+".png")
      cpull=canvas_create(frame3,4.7,5.7,1,' Pull of N_{sgn}',False)
      cpull.SaveAs("cpull_"+outputfile+".png")
      clog=canvas_create(frame4,4.7,5.7,1,'- log (L)')
      clog.SaveAs("clog_"+outputfile+".png")

      postfit_data = mctoys.fitParDataSet()
      postfit_nsig = np.array([postfit_data.get(i).getRealValue("nsig") for i in range(int(postfit_data.sumEntries()))])
      #postfit_nsig = np.array([get_visible_yield(obs, sig, postfit_data.get(i).getRealValue("nsig")) for i in range(int(postfit_data.sumEntries()))])
      #postfit_mu = np.array([s/nsig_visible for s in postfit_nsig])
      postfit_mu = np.array([s/set_sgn_yield for s in postfit_nsig])
      rms_mu = np.std(postfit_mu)
      fig, ax = plt.subplots()
      ax.hist(postfit_mu, bins=50, normed=True, histtype='step', label='MVA={}, RMS={}'.format(mva, rms_mu))
      ax.set_xlabel(r'$\mu$')
      ax.set_ylabel('a.u.')
      ax.legend(loc='best')
      fig.savefig('cmu_{}.png'.format(outputfile), bbox_inches='tight')

   # Save workspace 
   wspace_output = ROOT.RooWorkspace('wspace')
   getattr(wspace_output, "import")(model)
   getattr(wspace_output, "import")(dataset)
   params = model.getParameters(dataset)
   wspace_output.saveSnapshot("nominal_values",params)
   wspace_output.Print("V")
   wspace_output.writeToFile('wspace_'+outputfile+'.root')

   csv_header = ['cut', 'nsig', 'nbkg', 'njpsi', 'snr', 'rms_mu']
   df = {}
   df['cut'] = mva
   df['nsig'] = nsig_visible
   df['nbkg'] = nbkg_visible
   df['njpsi'] = nkjpsi_visible
   df['snr'] = nsig_visible / np.sqrt(nsig_visible + nbkg_visible + nkjpsi_visible)
   df['rms_mu'] = 0.0 if number_of_mctoys == None else rms_mu_
   csv_outputfile = log
   file_exists = os.path.isfile(csv_outputfile)
   with open (csv_outputfile, 'a+') as filedata:                            
     writer = csv.DictWriter(filedata, delimiter=',', fieldnames=csv_header)
     if not file_exists:
       writer.writeheader()
     writer.writerow(df) 


   return (nsig_visible,nbkg_visible,nkjpsi_visible)


def kee_yield_from_kjpsi(kjpsi_yield,tree_kee,tree_kjpsi,total_kee,total_kjpsi,xgb):
    print "\nWARNING!!!!!",total_kee,"will be used as denominsator for kee eff.  and ",total_kjpsi,"for kjpsi. Is this correct ?; code hypothesizes that cross validation is in 8 parts (hardcoded, change if different)\n"
    hkee = ROOT.TH1F("hkee","",50,4.7,5.7)	
    hkjpsi = ROOT.TH1F("hkjpsi","",50,4.7,5.7)  
    #tree_kee.Draw("Bmass>>hkee","1./8.*(Bmass>4.7 && Bmass<5.7 && Mll>1.05 && Mll<2.45 && xgb>"+xgb+")")
    #tree_kjpsi.Draw("Bmass>>hkjpsi","1./8.*(Bmass>4.7 && Bmass<5.7 && Mll>2.8 && Mll<3.25 && xgb>"+xgb+")")
    tree_kee.Draw("Bmass>>hkee","1./8.*(Bmass>4.7 && Bmass<5.7 && Mll>1.05 && Mll<2.45 && ( (Npv<15 && xgb>8) || (Npv>14 && xgb>8.5 ) ))")
    tree_kjpsi.Draw("Bmass>>hkjpsi","1./8.*(Bmass>4.7 && Bmass<5.7 && Mll>2.8 && Mll<3.25 && ( (Npv<15 && xgb>8) || (Npv>14 && xgb>8.5 ) ))")

    eff_kee=float(hkee.Integral()) / total_kee
    eff_kjpsi=float(hkjpsi.Integral()) / total_kjpsi
    result = kjpsi_yield / eff_kjpsi * 4.43*0.01/(1.026*5.93) *eff_kee
    print "expect",result,"eff rare",eff_kee,"eff res",eff_kjpsi,"res data",kjpsi_yield
    print "rare eff num",float(hkee.Integral()),"den",total_kee
    return result

################################### for scan ############################ 
def FitForScan(inputfile, mva, isgnfile, ikspiBkg, ibkg, kjpsi_yield_for_kee,total_nsgn, total_nkjpsi, name ):
    branches=["Bmass","Mll","xgb"]
    cuts = "xgb>"+str(mva)+" && 1.05<Mll && Mll<2.45" 
    outputfile="scan_kee_wp"+str(mva)+"_"+name

    tree_sgn = ROOT.TChain('mytreefit')
    tree_sgn.Add(isgnfile)
    tree_sgn_cut=tree_sgn.CopyTree(cuts)
    signal_parameters = signal_fit(tree_sgn_cut, outputfile+"_sgnMC", branches)
      
    tree_kjpsi = ROOT.TChain('mytreefit')
    tree_kjpsi.Add(ikspiBkg)
    tree_kjpsi_cut=tree_kjpsi.CopyTree(cuts)
    kjpsi_parameters = kjpsi_fit(tree_kjpsi_cut, outputfile+"_kjpsiMC", branches)
    tree_bkg = ROOT.TChain('mytreefit')
    tree_bkg.Add(ibkg)
    cuts_bkg = "xgb>"+str(mva)+" && Mll<5"
    tree_bkg_cut=tree_bkg.CopyTree(cuts_bkg)
    bkg_parameters = bkg_fit(tree_bkg_cut, outputfile+"_SameSign", branches)

    tree_sgn = ROOT.TChain('mytreefit')
    tree_sgn.Add(isgnfile)
    tree_kjpsi = ROOT.TChain('mytreefit')
    tree_kjpsi.Add(ikspiBkg)
    set_expected_sgn = kee_yield_from_kjpsi(kjpsi_yield_for_kee,tree_sgn,tree_kjpsi,total_nsgn,total_nkjpsi,str(mva))
    print "WARNING expected signal set and FIXED to ",set_expected_sgn,"in all q^2"
    tree = ROOT.TChain('mytreefit')
    tree.Add(inputfile)
    tree_cut=tree.CopyTree(cuts)
    nsig, nbkg, nkjpsi =total_fit(tree_cut, outputfile, branches, signal_parameters,  kjpsi_parameters, bkg_parameters, set_expected_sgn, {"min":5.0,"max":5.4},None,str(mva))
    return {"sig":nsig, "bkg":nbkg+nkjpsi, "ntotal_sig":set_expected_sgn }
    



#################################### main ################################
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Unbinned likelihood fit')
    parser.add_argument("-i", "--inputfile", dest="inputfile", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_finalfix_12B_*_data.root", help="Input data file")
    parser.add_argument("-o", "--outputfile", dest="outputfile", default="test", help="Output name")
    parser.add_argument("--mvacut", dest="mva",default=7.2,type=float)
    parser.add_argument("--set_expected_sgn", dest="set_expected_sgn", default=None,type=float,  help="fixes kee yield in the number")
    parser.add_argument("--number_of_mctoys", dest="number_of_mctoys", default=None, type=int,  help="None to deactivate, else genearates as test the given points")
    parser.add_argument("--kjpsi_yield_for_kee", dest="kjpsi_yield_for_kee", default=None,type=float,  help="measured resonant yield in data for estimating kee")
    parser.add_argument("--isgn", "--isgnfile", dest="isgnfile", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_finalfix_12B_*_MC.root", help="Input signal file")
    parser.add_argument("--nsgn", dest="total_nsgn", default=0,type=float, help="Total number of MC evts for denominator of efficiency. DO NOT put N* for N x-validation segments, code divides by 8")
    parser.add_argument("--ikjpsi", "--ikjpsiBkg", dest="ikjpsiBkg", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_finalfix_12B_*_MCres.root", help="Input kjpsi BKG file")
    parser.add_argument("--ikstaree", "--ikstareeBkg", dest="ikstareeBkg", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_finalfix_12B_*_MCres.root", help="Input kstar ee BKG file")
    parser.add_argument("--nkjpsi", dest="total_nkjpsi", default=0,type=float, help="Total number of MC evts for denominator of efficiency. DO NOT put N* for N x-validation segments, code divides by 8")
    parser.add_argument("--ibkg", dest="ibkg", default="../BDT/TrainingCore/forMeas_xgbmodel_kee_finalfix_12B_samesign_model_*_chunk_*_data.root", help="Input combinatorial BKG file")
    parser.add_argument("--fit_primitive", dest="fit_primtv", default=False, action='store_true', help="primitive fits for fixing params or use defaults")
    parser.add_argument("--skip_realfit", dest="skip_realfit", default=False, action='store_true', help="does not perform the final fit. useful for defining parameters, tests on pdfs")
    parser.add_argument("--sel_primitive", dest="sel_primtv", default=None,  help="runs only the selected primitive fits. Options: sgn bkg_comb bkg_kjpsi. they can be combined with ',' in strings. No spaces")
    parser.add_argument("--log", dest="log", default="log.csv", help="log of the fitting results")
    args = parser.parse_args()

    branches=["Bmass","Mll","xgb","KLmassD0"] #+ ["L1pt", "L2pt"]
    cuts = branches[2]+">"+str(args.mva)+" && 1.05<"+branches[1]+" && "+branches[1]+"<2.45"+" && "+branches[3]+">2.0" #+" && {} > 10 && {} > 10".format(branches[4], branches[5])
    cuts_samesign = branches[2]+">"+str(args.mva)+" && 1.05<"+branches[1]+" && "+branches[1]+"<2.45"
    print "cut: ", cuts
    #cuts = " 1.05<"+branches[1]+" && "+branches[1]+"<2.45 && ( (Npv<15 && xgb>8) || (Npv>14 && xgb>8.5 ) )" 
    args.outputfile+="_wp"+str(args.mva)

    bkg_parameters={'exp_alpha':-1.98}

    print "start"
    if args.fit_primtv:
      if args.sel_primtv!= None:
         args.sel_primtv = args.sel_primtv.split(",")
      else:
         args.sel_primtv = ["sgn","bkg_comb","bkg_kjpsi"]
      print "primitive params"
      if "sgn" in args.sel_primtv:
        tree_sgn = ROOT.TChain('mytreefit')
        tree_sgn.Add(args.isgnfile)
        tree_sgn_cut=tree_sgn.CopyTree(cuts)
        signal_parameters = signal_fit(tree_sgn_cut, args.outputfile+"_sgnMC", branches)
        print "parameters SGN", signal_parameters['mean'], signal_parameters['width'],signal_parameters['alpha1'],signal_parameters['n1'],signal_parameters['alpha2'],signal_parameters['n2']
      
      if "bkg_kjpsi" in args.sel_primtv:
        tree_kjpsi = ROOT.TChain('mytreefit')
        tree_kjpsi.Add(args.ikjpsiBkg)
        tree_kjpsi_cut=tree_kjpsi.CopyTree(cuts)
        kjpsi_pdf = kde_fit(tree_kjpsi_cut, args.outputfile+"_kjpsiMC", branches, 'kjpsi')
        print "finished kde fit for K J/psi"
    
      if "bkg_kstaree" in args.sel_primtv:
        tree_kstaree = ROOT.TChain('mytreefit')
        tree_kstaree.Add(args.ikstareeBkg)
        tree_kstaree_cut=tree_kstaree.CopyTree(cuts)
        kstaree_pdf = kde_fit(tree_kstaree_cut, args.outputfile+"_kstareeMC", branches, 'kstaree')
        print "finished kde fit for K* ee"

      if "bkg_comb" in args.sel_primtv:
        #cuts = branches[2]+">8 && "+branches[1]+"<5."
        #cuts = branches[2]+">"+str(args.mva)+" && "+branches[1]+"<5"
        tree_bkg = ROOT.TChain('mytreefit')
        tree_bkg.Add(args.ibkg)
        tree_bkg_cut=tree_bkg.CopyTree(cuts_samesign)
        bkg_parameters = bkg_fit(tree_bkg_cut, args.outputfile+"_SameSign", branches)
        print "parameters Combinatorial BKG", bkg_parameters['exp_alpha']

    else:
      signal_parameters={'mean': 5.272, 'width': 0.057, 'alpha1': 0.652, 'n1': 3.3, 'alpha2': 1.32, 'n2': 2.01}
      kjpsi_parameters={'mean_kjpsi':4.72,'width_kjpsi':1.06}
      bkg_parameters={'exp_alpha':-1.98}

    if args.kjpsi_yield_for_kee!=None:
      tree_sgn = ROOT.TChain('mytreefit')
      tree_sgn.Add(args.isgnfile)
      tree_kjpsi = ROOT.TChain('mytreefit')
      tree_kjpsi.Add(args.ikspiBkg)
      args.set_expected_sgn = kee_yield_from_kjpsi(args.kjpsi_yield_for_kee,tree_sgn,tree_kjpsi,args.total_nsgn,args.total_nkjpsi,str(args.mva))
      print "WARNING expected signal set and FIXED to ",args.set_expected_sgn,"in all q^2"

    #args.set_expected_sgn = None
    if not args.skip_realfit:
      #cuts = branches[2]+">"+str(args.mva)+" && 1.05<"+branches[1]+" && "+branches[1]+"<2.45"+" && "+branches[3]+">2.0 &&"+branches[4]+">2.0"
      tree = ROOT.TChain('mytreefit')
      tree.Add(args.inputfile)
      tree_cut=tree.CopyTree(cuts)
      nsig, nbkg, nkjpsi =total_fit(tree_cut, args.outputfile, branches, signal_parameters, kjpsi_pdf, kstaree_pdf, bkg_parameters, args.set_expected_sgn, {"min":5.0,"max":5.4}, args.number_of_mctoys, str(args.mva), args.log)
      # combinatorial BKG parameters set but not fixed.
      print "sigma",float(nsig)/math.sqrt(nsig+nkjpsi+nbkg),"nsig",float(nsig),"nbkg",float(nbkg),"Kjpsi leak",float(nkjpsi)
