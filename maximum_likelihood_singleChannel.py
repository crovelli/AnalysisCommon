import os,sys,re
from time import sleep
import math
import ROOT as rt
rt.gROOT.SetBatch(True);
rt.gROOT.SetStyle("Plain");
from plotting import *

def importWorkspace(wspace, wsfilename_import, rename):
    wsfile_import = rt.TFile.Open(wsfilename_import)
    wspace_import = wsfile_import.Get('wspace')
    data_import = wspace_import.data('data')
    getattr(wspace, 'import')(data_import, rt.RooCmdArg(RooFit.Rename('data_{}'.format(rename))))
    model_import = wspace_import.pdf('model')
    getattr(wspace, 'import')(model_import, RooFit.RenameAllVariablesExcept(rename, 'x'), RooFit.RenameAllNodes(rename))
    return wspace

def createWorkspace(wsname, wsfilename, wsfilename_tar, wsfilename_norm, eff_tar=0.0, deff_tar=0.0, eff_norm=0.0, deff_norm=0.0):
    wspace = rt.RooWorkspace(wsname)
    wspace = importWorkspace(wspace, wsfilename_tar, 'tar')
    wspace = importWorkspace(wspace, wsfilename_norm, 'norm')

    wspace.factory('cat[tar,norm]')
    x = wspace.var("x")
    data = ROOT.RooDataSet('data', 'data', rt.RooArgSet(x), RooFit.Index(wspace.cat('cat')), RooFit.Import('tar', wspace.data('data_tar')), RooFit.Import('norm', wspace.data('data_norm')))
    getattr(wspace, 'import')(data, rt.RooCmdArg())


    ### Create parameters
    # observations
    params = [
    # efficiency of target channel estimate
              ('eff_hat_tar',   eff_tar,   0,   1),
              ('deff_tar',     deff_tar,   0,   1),
    # efficiency of normalization channel estimate
              ('eff_hat_norm', eff_norm,   0,   1),
              ('deff_norm',   deff_norm,   0,   1),
    # nuisance parameters
              ('n_norm',           9000,   0,  20000),
              ('eff_tar',       eff_tar,   0,  1),
              ('eff_norm',     eff_norm,   0,  1),
    # parameter of interest
              ('R',      0.09,    0.06,   0.15)]

    for t in params:
        cmd = '%s[%f, %f, %f]' % t
        wspace.factory(cmd)        
    wspace.var('R').SetTitle('R')

    # fix all background and signal parameters
    for t in params[0:-4]:
        name = t[0]
        print '=> make %8s = %5.5f constant' % (name,
                                                wspace.var(name).getVal())
        wspace.var(name).setConstant()

    ### Create expressions
    express = ['n_tar("R*n_norm*(eff_tar/eff_norm)", R, n_norm, eff_tar, eff_norm)',
               'nKstarPsi2S_tar_modified("frac_partial_tar*n_tar", frac_partial_tar, n_tar)',
               'nKstarJpsi_norm_modified("frac_partial_norm*n_norm", frac_partial_norm, n_norm)',
              ]
        
    for t in express:
        cmd = 'expr::%s' % t
        wspace.factory(cmd)

    edit = ['model_tar_modified(model_tar, nsignal_tar=n_tar, nKstarPsi2S_tar=nKstarPsi2S_tar_modified)',
            'model_norm_modified(model_norm, nsignal_norm=n_norm, nKstarJpsi_norm=nKstarJpsi_norm_modified)',
           ]

    for t in edit:
        cmd = 'EDIT::%s' % t
        wspace.factory(cmd)

    wspace.Print("V")
    ### Create pdfs
    pdfs = [('Gaussian','peff_tar', '(eff_hat_tar, eff_tar, deff_tar)'),
            ('Gaussian','peff_norm', '(eff_hat_norm, eff_norm, deff_norm)'),
           ]
    
    prodpdf = ''
    for t in pdfs:
        wspace.factory('%s::%s%s' % t)
        name = t[1]
        prodpdf += "%s, " % name
    prodpdf = prodpdf[:-2] # remove last ", "
    

    # multiply the pdfs together. use upper case PROD to
    # do this
    wspace.factory("SIMUL:jointModel(cat,tar=model_tar_modified,norm=model_norm_modified)")
    wspace.factory('PROD::model({},{})'.format(prodpdf, 'jointModel'))

    nuis_excluded = ['x', 'R', 'cat', 'deff_norm', 'deff_tar', 'eff_hat_norm', 'eff_hat_tar']
    nuis = []
    params = wspace.pdf('model').getVariables()
    params_iter = params.createIterator()
    param = params_iter.Next()
    while param :
      if param.GetName() not in nuis_excluded:
        nuis.append(param.GetName())
      param = params_iter.Next()
    nuis = ','.join(nuis)

    sets = [('obs',  'x'),           # observations
            ('poi',  'R'),          # parameter of interest
            ('nuis', nuis)] # nuisance parameters (leave no spaces)
            #('nuis', 'n_norm,eff_tar,eff_norm')] # nuisance parameters (leave no spaces)

    for t in sets:
        name, parlist = t
        wspace.defineSet(name, parlist)
    
      
    #-----------------------------------------------------
    # Create model configuration. This is needed for the
    # statistical analyses
    #-----------------------------------------------------
    cfg = rt.RooStats.ModelConfig('cfg')
    cfg.SetWorkspace(wspace)
    cfg.SetPdf(wspace.pdf('model'))
    cfg.SetParametersOfInterest(wspace.set('poi'))
    cfg.SetNuisanceParameters(wspace.set('nuis'))

    # import model configuration into workspace
    getattr(wspace, 'import')(cfg)

    wspace.Print()
    
    # write out workspace
    wspace.writeToFile(wsfilename)

def analyzeWorkspace(wsname, wsfilename, R_SM=None):

    # Open workspace file
    wsfile = rt.TFile.Open(wsfilename)

    # Get workspace
    wspace = wsfile.Get(wsname) 

    # Get data
    data = wspace.data('data')

    # Get model configuration    
    cfg  = wspace.obj('cfg')

    #-----------------------------------------------------    
    # Fit model to data
    #-----------------------------------------------------
    results = wspace.pdf('model').fitTo(data, RooFit.Extended(True), rt.RooFit.Save(), RooFit.NumCPU(8))
    results.Print()
    
    #-----------------------------------------------------    
    # Compute interval based on profile likelihood
    #-----------------------------------------------------
    # suppress some (apparently) innocuous warnings
    #msgservice = rt.RooMsgService.instance()
    #msgservice.setGlobalKillBelow(rt.RooFit.FATAL)
       
    R_min, R_max = 0.07, 0.12
    NPoints = 20
    xlabel = 'R_{#psi (2S)}^{e}'

    print 'compute interval using profile likelihood'
    plc = rt.RooStats.ProfileLikelihoodCalculator(data, cfg)
    CL  = 0.683
    plc.SetConfidenceLevel(CL)
    plcInterval= plc.GetInterval()
    lowerLimit = plcInterval.LowerLimit(wspace.var('R'))
    upperLimit = plcInterval.UpperLimit(wspace.var('R'))

    print '\tPL %4.1f%s CL interval = [%5.5f, %5.5f]' % \
      (100*CL, '%', lowerLimit, upperLimit)

    plcplot = rt.RooStats.LikelihoodIntervalPlot(plcInterval)      
    plccanvas = rt.TCanvas('fig_PL', 'PL', 800, 600)
    plccanvas.cd()
    pad = setup_pad()
    pad.Draw()
    pad.cd()

    plcplot.SetRange(R_min, R_max)
    plcplot.SetMaximum(10)
    plcplot.SetLineColor(9)
    plcplot.SetNPoints(NPoints)
    plcplot.Draw("tf1")

    # compute an 95% limit on mu by
    CL = 0.95
    plc.SetConfidenceLevel(CL)
    plcInterval = plc.GetInterval()
    lowerLimit = plcInterval.LowerLimit(wspace.var('R'))
    upperLimit = plcInterval.UpperLimit(wspace.var('R'))

    print '\tPL %4.1f%s CL interval = [%5.5f, %5.5f]' % \
      (100*CL, '%', lowerLimit, upperLimit)
      
    pad.cd()
    plcplot2 = rt.RooStats.LikelihoodIntervalPlot(plcInterval)
    plcplot2.SetRange(R_min, R_max)
    plcplot2.SetMaximum(10)
    plcplot2.SetLineColor(8)
    plcplot2.SetNPoints(NPoints)
    plcplot2.Draw("tf1 same")
    
    frame = plcplot.GetPlottedObject()
    frame.GetYaxis().SetTitle('Profile of -log(L/L_{min})')
    frame.GetXaxis().SetTitle(xlabel)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleFont(42)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelFont(42)
    frame.GetXaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetTitleFont(42)
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetLabelFont(42)

    pad.Update()
    uymax = rt.gPad.GetUymax()
    uymin = rt.gPad.GetUymin()
    if R_SM is not None:
      pad.cd()
      l = rt.TLine(R_SM, uymin, R_SM, uymax)
      l.SetLineColor(2)
      l.SetLineWidth(2)
      l.Draw("same")

    pad.cd()
    CMS_lumi(False)
    plccanvas.cd()
    plccanvas.Update()

    # save canvases
    plccanvas.Draw()
    plccanvas.SaveAs('.pdf')
    return plccanvas

if __name__ == "__main__":

    eff_tar, deff_tar = 0.03430, 0.00024
    eff_norm, deff_norm = 0.03680, 0.00023
    R_SM = 0.0812
    
    wsfilename_tar = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_norm = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename = 'wspace_rpsi2s_electron.root'

    createWorkspace('R_psi2s', wsfilename, wsfilename_tar, wsfilename_norm, eff_tar=eff_tar, deff_tar=deff_tar, eff_norm=eff_norm, deff_norm=deff_norm)
    plccanvas = analyzeWorkspace('R_psi2s', wsfilename, R_SM=R_SM)


