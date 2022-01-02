import os,sys,re
import time
from time import sleep
import math
import ROOT as rt
from ROOT import RooFit
rt.gROOT.SetBatch(True);
rt.gROOT.SetStyle("Plain");
from utils_likelihood import analyzeWorkspace, plot_likelihood, importWorkspace 

def createWorkspace(wsname, wsfilename, 
                    wsfilename_pfe_psi2s, wsfilename_pfe_jpsi, 
                    wsfilename_mixe_psi2s, wsfilename_mixe_jpsi, 
                    wsfilename_mu_psi2s, wsfilename_mu_jpsi,
                    eff_pfe_psi2s=0.0, deff_pfe_psi2s=0.0, eff_pfe_jpsi=0.0, deff_pfe_jpsi=0.0,
                    eff_mixe_psi2s=0.0, deff_mixe_psi2s=0.0, eff_mixe_jpsi=0.0, deff_mixe_jpsi=0.0,
                    eff_mu_psi2s=0.0, deff_mu_psi2s=0.0, eff_mu_jpsi=0.0, deff_mu_jpsi=0.0,
                    ):

    wspace = rt.RooWorkspace(wsname)
    wspace = importWorkspace(wspace, wsfilename_pfe_psi2s, 'pfe_psi2s')
    wspace = importWorkspace(wspace, wsfilename_pfe_jpsi,  'pfe_jpsi')
    wspace = importWorkspace(wspace, wsfilename_mixe_psi2s, 'mixe_psi2s')
    wspace = importWorkspace(wspace, wsfilename_mixe_jpsi,  'mixe_jpsi')
    wspace = importWorkspace(wspace, wsfilename_mu_psi2s,  'mu_psi2s')
    wspace = importWorkspace(wspace, wsfilename_mu_jpsi,   'mu_jpsi')

    wspace.factory('cat[pfe_psi2s,pfe_jpsi,mixe_psi2s,mixe_jpsi,mu_psi2s,mu_jpsi]')
    x = wspace.var("x")
    data = rt.RooDataSet('data', 'data', rt.RooArgSet(x), RooFit.Index(wspace.cat('cat')), 
                           RooFit.Import('pfe_psi2s', wspace.data('data_pfe_psi2s')), RooFit.Import('pfe_jpsi', wspace.data('data_pfe_jpsi')), 
                           RooFit.Import('mixe_psi2s', wspace.data('data_mixe_psi2s')), RooFit.Import('mixe_jpsi', wspace.data('data_mixe_jpsi')), 
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
    # efficiency of PF-PF electron psi(2S) estimate
              ('eff_hat_pfe_psi2s',   eff_pfe_psi2s,  1.e-6,   1),
              ('deff_pfe_psi2s',      deff_pfe_psi2s, 1.e-6,   1),
    # efficiency of PF-PF electron J/psi channel estimate
              ('eff_hat_pfe_jpsi',    eff_pfe_jpsi,   1.e-6,   1),
              ('deff_pfe_jpsi',       deff_pfe_jpsi,  1.e-6,   1),
    # efficiency of PF-LP electron psi(2S) estimate
              ('eff_hat_mixe_psi2s',  eff_mixe_psi2s,  1.e-6,   1),
              ('deff_mixe_psi2s',     deff_mixe_psi2s, 1.e-6,   1),
    # efficiency of PF-LP electron J/psi channel estimate
              ('eff_hat_mixe_jpsi',   eff_mixe_jpsi,   1.e-6,   1),
              ('deff_mixe_jpsi',      deff_mixe_jpsi,  1.e-6,   1),
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
              ('n_pfe_jpsi',          8617,           7000,   10000),
              ('eff_pfe_psi2s',       eff_pfe_psi2s,  eff_pfe_psi2s-5.0*deff_pfe_psi2s,     eff_pfe_psi2s+5.0*deff_pfe_psi2s),
              ('eff_pfe_jpsi',        eff_pfe_jpsi,   eff_pfe_jpsi-5.0*deff_pfe_jpsi,       eff_pfe_jpsi+5.0*deff_pfe_jpsi),
              ('n_mixe_jpsi',         8617,           7000,   10000),
              ('eff_mixe_psi2s',      eff_mixe_psi2s, eff_mixe_psi2s-5.0*deff_mixe_psi2s,   eff_mixe_psi2s+5.0*deff_mixe_psi2s),
              ('eff_mixe_jpsi',       eff_mixe_jpsi,  eff_mixe_jpsi-5.0*deff_mixe_jpsi,     eff_mixe_jpsi+5.0*deff_mixe_jpsi),
              ('br_psi2s_ee',         br_psi2s_ee,    br_psi2s_ee-5.0*dbr_psi2s_ee,         br_psi2s_ee+5.0*dbr_psi2s_ee),
              ('br_jpsi_ee',          br_jpsi_ee,     br_jpsi_ee-5.0*dbr_jpsi_ee,           br_jpsi_ee+5.0*dbr_jpsi_ee),
              ('n_mu_jpsi',           8617,           7000,   10000),
              ('eff_mu_psi2s',        eff_mu_psi2s,   eff_mu_psi2s-5.0*deff_mu_psi2s,       eff_mu_psi2s+5.0*deff_mu_psi2s),
              ('eff_mu_jpsi',         eff_mu_jpsi,    eff_mu_jpsi-5.0*deff_mu_jpsi,         eff_mu_jpsi+5.0*deff_mu_jpsi),
              ('br_psi2s_mumu',       br_psi2s_mumu,  br_psi2s_mumu-5.0*dbr_psi2s_mumu,     br_psi2s_mumu+5.0*dbr_psi2s_mumu),
              ('br_jpsi_mumu',        br_jpsi_mumu,   br_jpsi_mumu-5.0*dbr_jpsi_mumu,       br_jpsi_mumu+5.0*dbr_jpsi_mumu),
    # parameter of interest
              ('r_br',                0.7,            0.6,     0.8),
              ('R',                   1.0,            0.5,     1.5),
              ]

    for t in params:
        cmd = '%s[%f, %f, %f]' % t
        wspace.factory(cmd)        
    wspace.var('R').SetTitle('R')

    # fix all background and signal parameters
    for t in params[0:-15]:
        name = t[0]
        print '=> make %8s = %5.5f constant' % (name,
                                                wspace.var(name).getVal())
        wspace.var(name).setConstant()

    ### Create expressions
    express = [
    # PF-PF electron
               'n_pfe_psi2s("(r_br/R)*n_pfe_jpsi*(br_psi2s_ee/br_jpsi_ee)*(eff_pfe_psi2s/eff_pfe_jpsi)", r_br, R, n_pfe_jpsi, br_psi2s_ee, br_jpsi_ee, eff_pfe_psi2s, eff_pfe_jpsi)',
               'nKstarPsi2S_pfe_psi2s_modified("frac_partial_pfe_psi2s*n_pfe_psi2s", frac_partial_pfe_psi2s, n_pfe_psi2s)',
               'nKstarJpsi_pfe_jpsi_modified("frac_partial_pfe_jpsi*n_pfe_jpsi", frac_partial_pfe_jpsi, n_pfe_jpsi)',
    # PF-LP electron
               'n_mixe_psi2s("(r_br/R)*n_mixe_jpsi*(br_psi2s_ee/br_jpsi_ee)*(eff_mixe_psi2s/eff_mixe_jpsi)", r_br, R, n_mixe_jpsi, br_psi2s_ee, br_jpsi_ee, eff_mixe_psi2s, eff_mixe_jpsi)',
               'nKstarPsi2S_mixe_psi2s_modified("frac_partial_mixe_psi2s*n_mixe_psi2s", frac_partial_mixe_psi2s, n_mixe_psi2s)',
               'nKstarJpsi_mixe_jpsi_modified("frac_partial_mixe_jpsi*n_mixe_jpsi", frac_partial_mixe_jpsi, n_mixe_jpsi)',
    # muon
               'n_mu_psi2s("r_br*n_mu_jpsi*(br_psi2s_mumu/br_jpsi_mumu)*(eff_mu_psi2s/eff_mu_jpsi)", r_br, n_mu_jpsi, br_psi2s_mumu, br_jpsi_mumu, eff_mu_psi2s, eff_mu_jpsi)',
               'nKstarPsi2S_mu_psi2s_modified("frac_partial_mu_psi2s*n_mu_psi2s", frac_partial_mu_psi2s, n_mu_psi2s)',
               'nKstarJpsi_mu_jpsi_modified("frac_partial_mu_jpsi*n_mu_jpsi", frac_partial_mu_jpsi, n_mu_jpsi)',
              ]
        
    for t in express:
        cmd = 'expr::%s' % t
        wspace.factory(cmd)

    edit = [
    # PF-PF electron
            'model_pfe_psi2s_modified(model_pfe_psi2s, nsignal_pfe_psi2s=n_pfe_psi2s, nKstarPsi2S_pfe_psi2s=nKstarPsi2S_pfe_psi2s_modified)',
            'model_pfe_jpsi_modified(model_pfe_jpsi, nsignal_pfe_jpsi=n_pfe_jpsi, nKstarJpsi_pfe_jpsi=nKstarJpsi_pfe_jpsi_modified)',
    # PF-LP electron
            'model_mixe_psi2s_modified(model_mixe_psi2s, nsignal_mixe_psi2s=n_mixe_psi2s, nKstarPsi2S_mixe_psi2s=nKstarPsi2S_mixe_psi2s_modified)',
            'model_mixe_jpsi_modified(model_mixe_jpsi, nsignal_mixe_jpsi=n_mixe_jpsi, nKstarJpsi_mixe_jpsi=nKstarJpsi_mixe_jpsi_modified)',
    # muon
            'model_mu_psi2s_modified(model_mu_psi2s, nsignal_mu_psi2s=n_mu_psi2s, nKstarPsi2S_mu_psi2s=nKstarPsi2S_mu_psi2s_modified)',
            'model_mu_jpsi_modified(model_mu_jpsi, nsignal_mu_jpsi=n_mu_jpsi, nKstarJpsi_mu_jpsi=nKstarJpsi_mu_jpsi_modified)',
           ]

    for t in edit:
        cmd = 'EDIT::%s' % t
        wspace.factory(cmd)

    wspace.Print("V")
    ### Create pdfs
    pdfs = [
    # PF-PF electron
            ('Gaussian','peff_pfe_psi2s', '(eff_hat_pfe_psi2s, eff_pfe_psi2s, deff_pfe_psi2s)'),
            ('Gaussian','peff_pfe_jpsi', '(eff_hat_pfe_jpsi, eff_pfe_jpsi, deff_pfe_jpsi)'),
    # PF-LP electron
            ('Gaussian','peff_mixe_psi2s', '(eff_hat_mixe_psi2s, eff_mixe_psi2s, deff_mixe_psi2s)'),
            ('Gaussian','peff_mixe_jpsi', '(eff_hat_mixe_jpsi, eff_mixe_jpsi, deff_mixe_jpsi)'),
            ('Gaussian','peff_br_psi2s_ee', '(br_hat_psi2s_ee, br_psi2s_ee, dbr_psi2s_ee)'),
            ('Gaussian','peff_br_jpsi_ee', '(br_hat_jpsi_ee, br_jpsi_ee, dbr_jpsi_ee)'),
    # muon
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
    wspace.factory("SIMUL:jointModel(cat,pfe_psi2s=model_pfe_psi2s_modified,pfe_jpsi=model_pfe_jpsi_modified,mixe_psi2s=model_mixe_psi2s_modified,mixe_jpsi=model_mixe_jpsi_modified,mu_psi2s=model_mu_psi2s_modified,mu_jpsi=model_mu_jpsi_modified)")
    wspace.factory('PROD::model({},{})'.format(prodpdf, 'jointModel'))

    nuis_excluded = ['x', 'R', 'cat', 
                     'deff_pfe_jpsi', 'deff_pfe_psi2s', 'eff_hat_pfe_jpsi', 'eff_hat_pfe_psi2s', 
                     'deff_mixe_jpsi', 'deff_mixe_psi2s', 'eff_hat_mixe_jpsi', 'eff_hat_mixe_psi2s', 
                     'dbr_jpsi_ee', 'dbr_psi2s_ee', 'br_hat_jpsi_ee', 'br_hat_psi2s_ee'
                     'deff_mu_jpsi', 'deff_mu_psi2s', 'eff_hat_mu_jpsi', 'eff_hat_mu_psi2s', 
                     'dbr_jpsi_mumu', 'dbr_psi2s_mumu', 'br_hat_jpsi_mumu', 'br_hat_psi2s_mumu'
                     ]

    ### Define global observables
    ### they are not being fitted and they are not loaded from a dataset, 
    ### but some knowledge exists that allows to set them to a specific value
    ### Global Observables are generated once per toy
    globs = ['eff_pfe_psi2s', 'eff_pfe_jpsi', 'br_psi2s_ee', 'br_jpsi_ee',
             'eff_mixe_psi2s', 'eff_mixe_jpsi',
             'eff_mu_psi2s', 'eff_mu_jpsi', 'br_psi2s_mumu', 'br_jpsi_mumu',
             ]

    nuis = []
    params = wspace.pdf('model').getVariables()
    params_iter = params.createIterator()
    param = params_iter.Next()
    while param :
      if param.GetName() not in (nuis_excluded+globs):
        nuis.append(param.GetName())
      param = params_iter.Next()
    nuis = ','.join(nuis)
    globs = ','.join(globs)

    sets = [('obs',   'x'),           # observations
            ('poi',   'R'),          # parameter of interest
            ('nuis',  nuis),
            ('globs', globs),
            ] # nuisance parameters (leave no spaces)
            #('nuis', 'n_pfe_jpsi,eff_pfe_psi2s,eff_pfe_jpsi')] # nuisance parameters (leave no spaces)

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
    cfg.SetGlobalObservables(wspace.set('globs'))

    # import model configuration into workspace
    getattr(wspace, 'import')(cfg)

    wspace.Print()
    
    # write out workspace
    wspace.writeToFile(wsfilename)


if __name__ == "__main__":

    eff_pfe_psi2s, deff_pfe_psi2s = 0.03430, 0.00024
    eff_pfe_jpsi, deff_pfe_jpsi = 0.03680, 0.00023
    eff_mixe_psi2s, deff_mixe_psi2s = 0.03430, 0.00024
    eff_mixe_jpsi, deff_mixe_jpsi = 0.03680, 0.00023
    eff_mu_psi2s, deff_mu_psi2s = 0.03430, 0.00024
    eff_mu_jpsi, deff_mu_jpsi = 0.03680, 0.00023
    R_SM = 1.0
    R_min, R_max = 0.8, 1.2
    xlabel = 'R_{#psi(2S)}'

    wsfilename_pfe_psi2s = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_pfe_jpsi = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename_mixe_psi2s = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_mixe_jpsi = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename_mu_psi2s = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_mu_jpsi = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename = 'wspace_rpsi2s_combinedEleCh.root'

    createWorkspace('R_psi2s', wsfilename, wsfilename_pfe_psi2s, wsfilename_pfe_jpsi, wsfilename_mixe_psi2s, wsfilename_mixe_jpsi, wsfilename_mu_psi2s, wsfilename_mu_jpsi,
                    eff_pfe_psi2s=eff_pfe_psi2s, deff_pfe_psi2s=deff_pfe_psi2s, eff_pfe_jpsi=eff_pfe_jpsi, deff_pfe_jpsi=deff_pfe_jpsi,
                    eff_mixe_psi2s=eff_mixe_psi2s, deff_mixe_psi2s=deff_mixe_psi2s, eff_mixe_jpsi=eff_mixe_jpsi, deff_mixe_jpsi=deff_mixe_jpsi,
                    eff_mu_psi2s=eff_mu_psi2s, deff_mu_psi2s=deff_mu_psi2s, eff_mu_jpsi=eff_mu_jpsi, deff_mu_jpsi=deff_mu_jpsi,
                    )
    plccanvas = analyzeWorkspace('R_psi2s', wsfilename, x_SM=R_SM, x_min=R_min, x_max=R_max, xlabel=xlabel, plot=True)


