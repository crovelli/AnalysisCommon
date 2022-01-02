import os,sys,re
from time import sleep
import math
import ROOT as rt
from ROOT import RooFit
rt.gROOT.SetBatch(True);
rt.gROOT.SetStyle("Plain");
from utils_likelihood import analyzeWorkspace, plot_likelihood, importWorkspace 

def createWorkspace(wsname, wsfilename, wsfilename_tar, wsfilename_norm, eff_tar=0.0, deff_tar=0.0, eff_norm=0.0, deff_norm=0.0):
    wspace = rt.RooWorkspace(wsname)
    wspace = importWorkspace(wspace, wsfilename_tar, 'tar')
    wspace = importWorkspace(wspace, wsfilename_norm, 'norm')

    wspace.factory('cat[tar,norm]')
    x = wspace.var("x")
    data = rt.RooDataSet('data', 'data', rt.RooArgSet(x), RooFit.Index(wspace.cat('cat')), RooFit.Import('tar', wspace.data('data_tar')), RooFit.Import('norm', wspace.data('data_norm')))
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
              ('n_norm',           8617,   7000,   10000),
              ('eff_tar',       eff_tar,   eff_tar-5.0*deff_tar,  eff_tar+5.0*deff_tar),
              ('eff_norm',     eff_norm,   eff_norm-5.0*deff_norm,  eff_norm+5.0*deff_norm),
    # parameter of interest
              ('R',      0.09,    0.01,   0.2)]

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

    ### Define global observables
    ### they are not being fitted and they are not loaded from a dataset, 
    ### but some knowledge exists that allows to set them to a specific value
    ### Global Observables are generated once per toy
    globs = ['eff_tar', 'eff_norm']

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
            ('nuis',  nuis), # nuisance parameters (leave no spaces)
            ('globs', globs),
            ]
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
    cfg.SetGlobalObservables(wspace.set('globs'))

    # import model configuration into workspace
    getattr(wspace, 'import')(cfg)

    wspace.Print()
    
    # write out workspace
    wspace.writeToFile(wsfilename)


if __name__ == "__main__":

    eff_tar, deff_tar = 0.03430, 0.00024
    eff_norm, deff_norm = 0.03680, 0.00023
    R_SM = 0.0812
    R_min, R_max = 0.07, 0.12
    xlabel = 'R_{#psi (2S)}^{e}'

    wsfilename_tar = 'wspace_psi2s_fixedPartial_pf_wp5.0.root'
    wsfilename_norm = 'wspace_jpsi_fixedPartial_pf_wp5.0.root'
    wsfilename = 'wspace_rpsi2s_electron.root'

    createWorkspace('R_psi2s', wsfilename, wsfilename_tar, wsfilename_norm, eff_tar=eff_tar, deff_tar=deff_tar, eff_norm=eff_norm, deff_norm=deff_norm)
    plccanvas = analyzeWorkspace('R_psi2s', wsfilename, x_SM=R_SM, x_min=R_min, x_max=R_max, xlabel=xlabel, plot=True)


