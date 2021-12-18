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

def createWorkspace(wsname, wsfilename, wsfilename_ele_psi2s, wsfilename_ele_jpsi, wsfilename_mu_psi2s, wsfilename_mu_jpsi,
                    eff_ele_psi2s=0.0, deff_ele_psi2s=0.0, eff_ele_jpsi=0.0, deff_ele_jpsi=0.0,
                    eff_mu_psi2s=0.0, deff_mu_psi2s=0.0, eff_mu_jpsi=0.0, deff_mu_jpsi=0.0,
                    ):

    wspace = rt.RooWorkspace(wsname)
    wspace = importWorkspace(wspace, wsfilename_ele_psi2s, 'ele_psi2s')
    wspace = importWorkspace(wspace, wsfilename_ele_jpsi,  'ele_jpsi')
    wspace = importWorkspace(wspace, wsfilename_mu_psi2s,  'mu_psi2s')
    wspace = importWorkspace(wspace, wsfilename_mu_jpsi,   'mu_jpsi')

    wspace.factory('cat[ele_psi2s,ele_jpsi,mu_psi2s,mu_jpsi]')
    x = wspace.var("x")
    data = ROOT.RooDataSet('data', 'data', rt.RooArgSet(x), RooFit.Index(wspace.cat('cat')), 
                           RooFit.Import('ele_psi2s', wspace.data('data_ele_psi2s')), RooFit.Import('ele_jpsi', wspace.data('data_ele_jpsi')), 
                           RooFit.Import('mu_psi2s', wspace.data('data_mu_psi2s')), RooFit.Import('mu_jpsi', wspace.data('data_mu_jpsi')),
                           )
    getattr(wspace, 'import')(data, rt.RooCmdArg())

    br_jpsi_ee, dbr_jpsi_ee= 0.05971, 0.00032 
    br_psi2s_ee, dbr_psi2s_ee = 7.93e-3, 0.17e-3 
    br_jpsi_mumu, dbr_jpsi_mumu= 0.05961, 0.00033 
    br_psi2s_mumu, dbr_psi2s_mumu = 8.0e-3, 0.6e-3 


    ### Create parameters
    # observations
    params = [
    # efficiency of electron psi(2S) estimate
              ('eff_hat_ele_psi2s',   eff_ele_psi2s,  1.e-6,   1),
              ('deff_ele_psi2s',      deff_ele_psi2s, 1.e-6,   1),
    # efficiency of electron J/psi channel estimate
              ('eff_hat_ele_jpsi',    eff_ele_jpsi,   1.e-6,   1),
              ('deff_ele_jpsi',       deff_ele_jpsi,  1.e-6,   1),
    # branching fraction of J/psi -> ee
              ('br_hat_jpsi_ee',      br_jpsi_ee,     1.e-6,   1),
              ('dbr_jpsi_ee',         dbr_jpsi_ee,    1.e-6,   1),
    # branching fraction of psi(2S) -> ee
              ('br_hat_psi2s_ee',     br_psi2s_ee,    1.e-6,   1),
              ('dbr_psi2s_ee',        dbr_psi2s_ee,   1.e-6,   1),
    # efficiency of muon psi(2S) estimate
              ('eff_hat_mu_psi2s',    eff_mu_psi2s,   1.e-6,   1),
              ('deff_mu_psi2s',       deff_mu_psi2s,  1.e-6,   1),
    # efficiency of muon J/psi channel estimate
              ('eff_hat_mu_jpsi',     eff_mu_jpsi,    1.e-6,   1),
              ('deff_mu_jpsi',        deff_mu_jpsi,   1.e-6,   1),
    # branching fraction of J/psi -> mu mu
              ('br_hat_jpsi_mumu',    br_jpsi_mumu,   1.e-6,   1),
              ('dbr_jpsi_mumu',       dbr_jpsi_mumu,  1.e-6,   1),
    # branching fraction of psi(2S) -> mu mu
              ('br_hat_psi2s_mumu',   br_psi2s_mumu,  1.e-6,   1),
              ('dbr_psi2s_mumu',      dbr_psi2s_mumu, 1.e-6,   1),
    # nuisance parameters
              ('n_ele_jpsi',          9000,           1.e-6,   20000),
              ('eff_ele_psi2s',       eff_ele_psi2s,  1.e-6,   1),
              ('eff_ele_jpsi',        eff_ele_jpsi,   1.e-6,   1),
              ('br_psi2s_ee',         br_psi2s_ee,    1.e-6,   1),
              ('br_jpsi_ee',          br_jpsi_ee,     1.e-6,   1),
              ('n_mu_jpsi',           9000,           1.e-6,   20000),
              ('eff_mu_psi2s',        eff_mu_psi2s,   1.e-6,   1),
              ('eff_mu_jpsi',         eff_mu_jpsi,    1.e-6,   1),
              ('br_psi2s_mumu',       br_psi2s_mumu,  1.e-6,   1),
              ('br_jpsi_mumu',        br_jpsi_mumu,   1.e-6,   1),
    # parameter of interest
              ('r_br',                0.7,            0.5,     1.1),
              ('R',                   1.0,            0.5,     2),
              ]

    for t in params:
        cmd = '%s[%f, %f, %f]' % t
        wspace.factory(cmd)        
    wspace.var('R').SetTitle('R')

    # fix all background and signal parameters
    for t in params[0:-12]:
        name = t[0]
        print '=> make %8s = %5.5f constant' % (name,
                                                wspace.var(name).getVal())
        wspace.var(name).setConstant()

    ### Create expressions
    express = ['n_ele_psi2s("(r_br/R)*n_ele_jpsi*(br_psi2s_ee/br_jpsi_ee)*(eff_ele_psi2s/eff_ele_jpsi)", r_br, R, n_ele_jpsi, br_psi2s_ee, br_jpsi_ee, eff_ele_psi2s, eff_ele_jpsi)',
               'nKstarPsi2S_ele_psi2s_modified("frac_partial_ele_psi2s*n_ele_psi2s", frac_partial_ele_psi2s, n_ele_psi2s)',
               'nKstarJpsi_ele_jpsi_modified("frac_partial_ele_jpsi*n_ele_jpsi", frac_partial_ele_jpsi, n_ele_jpsi)',
               'n_mu_psi2s("r_br*n_mu_jpsi*(br_psi2s_mumu/br_jpsi_mumu)*(eff_mu_psi2s/eff_mu_jpsi)", r_br, n_mu_jpsi, br_psi2s_mumu, br_jpsi_mumu, eff_mu_psi2s, eff_mu_jpsi)',
               'nKstarPsi2S_mu_psi2s_modified("frac_partial_mu_psi2s*n_mu_psi2s", frac_partial_mu_psi2s, n_mu_psi2s)',
               'nKstarJpsi_mu_jpsi_modified("frac_partial_mu_jpsi*n_mu_jpsi", frac_partial_mu_jpsi, n_mu_jpsi)',
              ]
        
    for t in express:
        cmd = 'expr::%s' % t
        wspace.factory(cmd)

    edit = ['model_ele_psi2s_modified(model_ele_psi2s, nsignal_ele_psi2s=n_ele_psi2s, nKstarPsi2S_ele_psi2s=nKstarPsi2S_ele_psi2s_modified)',
            'model_ele_jpsi_modified(model_ele_jpsi, nsignal_ele_jpsi=n_ele_jpsi, nKstarJpsi_ele_jpsi=nKstarJpsi_ele_jpsi_modified)',
            'model_mu_psi2s_modified(model_mu_psi2s, nsignal_mu_psi2s=n_mu_psi2s, nKstarPsi2S_mu_psi2s=nKstarPsi2S_mu_psi2s_modified)',
            'model_mu_jpsi_modified(model_mu_jpsi, nsignal_mu_jpsi=n_mu_jpsi, nKstarJpsi_mu_jpsi=nKstarJpsi_mu_jpsi_modified)',
           ]

    for t in edit:
        cmd = 'EDIT::%s' % t
        wspace.factory(cmd)

    wspace.Print("V")
    ### Create pdfs
    pdfs = [('Gaussian','peff_ele_psi2s', '(eff_hat_ele_psi2s, eff_ele_psi2s, deff_ele_psi2s)'),
            ('Gaussian','peff_ele_jpsi', '(eff_hat_ele_jpsi, eff_ele_jpsi, deff_ele_jpsi)'),
            ('Gaussian','peff_br_psi2s_ee', '(br_hat_psi2s_ee, br_psi2s_ee, dbr_psi2s_ee)'),
            ('Gaussian','peff_br_jpsi_ee', '(br_hat_jpsi_ee, br_jpsi_ee, dbr_jpsi_ee)'),
            ('Gaussian','peff_mu_psi2s', '(eff_hat_mu_psi2s, eff_mu_psi2s, deff_mu_psi2s)'),
            ('Gaussian','peff_mu_jpsi', '(eff_hat_mu_jpsi, eff_mu_jpsi, deff_mu_jpsi)'),
            ('Gaussian','peff_br_psi2s_mumu', '(br_hat_psi2s_mumu, br_psi2s_mumu, dbr_psi2s_mumu)'),
            ('Gaussian','peff_br_jpsi_mumu', '(br_hat_jpsi_mumu, br_jpsi_mumu, dbr_jpsi_mumu)'),
           ]
    
    prodpdf = ''
    for t in pdfs:
        wspace.factory('%s::%s%s' % t)
        name = t[1]
        prodpdf += "%s, " % name
    prodpdf = prodpdf[:-2] # remove last ", "
    

    # multiply the pdfs together. use upper case PROD to
    # do this
    wspace.factory("SIMUL:jointModel(cat,ele_psi2s=model_ele_psi2s_modified,ele_jpsi=model_ele_jpsi_modified,mu_psi2s=model_mu_psi2s_modified,mu_jpsi=model_mu_jpsi_modified)")
    wspace.factory('PROD::model({},{})'.format(prodpdf, 'jointModel'))

    nuis_excluded = ['x', 'R', 'cat', 
                     'deff_ele_jpsi', 'deff_ele_psi2s', 'eff_hat_ele_jpsi', 'eff_hat_ele_psi2s', 
                     'dbr_jpsi_ee', 'dbr_psi2s_ee', 'br_hat_jpsi_ee', 'br_hat_psi2s_ee'
                     'deff_mu_jpsi', 'deff_mu_psi2s', 'eff_hat_mu_jpsi', 'eff_hat_mu_psi2s', 
                     'dbr_jpsi_mumu', 'dbr_psi2s_mumu', 'br_hat_jpsi_mumu', 'br_hat_psi2s_mumu'
                     ]
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
            #('nuis', 'n_ele_jpsi,eff_ele_psi2s,eff_ele_jpsi')] # nuisance parameters (leave no spaces)

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
       
    R_min, R_max = 0.7, 1.1
    NPoints = 20
    xlabel = 'R_{#psi(2S)}'

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

    eff_ele_psi2s, deff_ele_psi2s = 0.03430, 0.00024
    eff_ele_jpsi, deff_ele_jpsi = 0.03680, 0.00023
    eff_mu_psi2s, deff_mu_psi2s = 0.03430, 0.00024
    eff_mu_jpsi, deff_mu_jpsi = 0.03680, 0.00023
    R_SM = 1.0
    
    wsfilename_ele_psi2s = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_ele_jpsi = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename_mu_psi2s = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_mu_jpsi = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename = 'wspace_rpsi2s.root'

    createWorkspace('R_psi2s', wsfilename, wsfilename_ele_psi2s, wsfilename_ele_jpsi, wsfilename_mu_psi2s, wsfilename_mu_jpsi,
                    eff_ele_psi2s=eff_ele_psi2s, deff_ele_psi2s=deff_ele_psi2s, eff_ele_jpsi=eff_ele_jpsi, deff_ele_jpsi=deff_ele_jpsi,
                    eff_mu_psi2s=eff_mu_psi2s, deff_mu_psi2s=deff_mu_psi2s, eff_mu_jpsi=eff_mu_jpsi, deff_mu_jpsi=deff_mu_jpsi,
                    )
    plccanvas = analyzeWorkspace('R_psi2s', wsfilename, R_SM=R_SM)


