import os,sys,re
from time import sleep
import math
import ROOT as rt
rt.gROOT.SetBatch(True);
rt.gROOT.SetStyle("Plain");
from plotting import *

def createWorkspace(wsname, wsfilename, N=0.0, n_norm=0.0, dn_norm=0.0, eff_tar=0.0, deff_tar=0.0, eff_norm=0.0, deff_norm=0.0):
    wspace = rt.RooWorkspace(wsname)

    ### Create parameters
    # observations
    params = [('N',       N,     0,  20000),
    # number of signal in normalization channel estimate        
              ('n_hat_norm',     n_norm,   0,  20000),
              ('dn_norm',       dn_norm,   0,   5000),
    # efficiency of target channel estimate
              ('eff_hat_tar',   eff_tar,   0,   1),
              ('deff_tar',     deff_tar,   0,   1),
    # efficiency of normalization channel estimate
              ('eff_hat_norm', eff_norm,   0,   1),
              ('deff_norm',   deff_norm,   0,   1),
    # nuisance parameters
              ('n_norm',         n_norm,   0,  20000),
              ('eff_tar',       eff_tar,   0,  1),
              ('eff_norm',     eff_norm,   0,  1),
    # parameter of interest
              ('R',      0.08,    0,   1)]

    for t in params:
        cmd = '%s[%f, %f, %f]' % t
        wspace.factory(cmd)        
    wspace.var('R').SetTitle('R')

    # fix all background and signal parameters
    for t in params[1:-4]:
        name = t[0]
        print '=> make %8s = %5.5f constant' % (name,
                                                wspace.var(name).getVal())
        wspace.var(name).setConstant()

    ### Create expressions
    express = ['M_norm("(n_hat_norm/dn_norm)^2", n_hat_norm, dn_norm)',
               'tau_norm("n_hat_norm/dn_norm^2", n_hat_norm, dn_norm)',
               'tau_norm_n("tau_norm*n_norm", tau_norm, n_norm)',
               'n("R*n_norm*(eff_tar/eff_norm)", R, n_norm, eff_tar, eff_norm)']
        
    for t in express:
        cmd = 'expr::%s' % t
        wspace.factory(cmd)

    ### Create pdfs
    pdfs = [('Poisson', 'pN',    '(N, n)'),
            ('Poisson', 'pM_norm', '(M_norm, tau_norm_n, 1)'), 
            ('Gaussian','peff_tar', '(eff_hat_tar, eff_tar, deff_tar)'),
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
    wspace.factory('PROD::model(%s)' % prodpdf)

    sets = [('obs',  'N'),           # observations
            ('poi',  'R'),          # parameter of interest
            ('nuis', 'n_norm,eff_tar,eff_norm')] # nuisance parameters (leave no spaces)
    for t in sets:
        name, parlist = t
        wspace.defineSet(name, parlist)
    
    #-----------------------------------------------------        
    # create a dataset
    #-----------------------------------------------------    
    data = rt.RooDataSet('data', 'data', wspace.set('obs'))
    data.add(wspace.set('obs'))
    # import dataset into workspace
    # need last argument to workaround a PyROOT "feature".
    # the last argument ensures the correct version of
    # the import method is called.
    getattr(wspace, 'import')(data, rt.RooCmdArg())
        
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
    results = wspace.pdf('model').fitTo(data, rt.RooFit.Save())
    results.Print()
    
    #-----------------------------------------------------    
    # Compute interval based on profile likelihood
    #-----------------------------------------------------
    # suppress some (apparently) innocuous warnings
    msgservice = rt.RooMsgService.instance()
    msgservice.setGlobalKillBelow(rt.RooFit.FATAL)
       
    R_min, R_max = 0.07, 0.12
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
    plcplot.Draw()

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
    plcplot2.Draw("same")
    
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

    if R_SM is not None:
      pad.cd()
      l = rt.TLine(R_SM, 0, R_SM, 10)
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

    N = 753.36
    n_norm, dn_norm = 8616.26, 147.97
    eff_tar, deff_tar = 0.03430, 0.00024
    eff_norm, deff_norm = 0.03680, 0.00023
    R_SM = 0.0812
    
    wsfilename = 'wspace_rpsi2s_electron.root'

    createWorkspace('R_psi2s', wsfilename, N=N, n_norm=n_norm, dn_norm=dn_norm, eff_tar=eff_tar, deff_tar=deff_tar, eff_norm=eff_norm, deff_norm=deff_norm)
    plccanvas = analyzeWorkspace('R_psi2s', wsfilename, R_SM=R_SM)


