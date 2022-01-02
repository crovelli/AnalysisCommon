import os,sys,re
import time
from time import sleep
import math
import ROOT as rt
import root_numpy
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

def analyzeWorkspace(wsname, wsfilename, x_SM=None, x_min=0.8, x_max=1.2, NPoints=20, xlabel='R', plot=False):
    rt.Math.MinimizerOptions.SetDefaultMinimizer("Minuit")
    #rt.Math.MinimizerOptions.SetDefaultPrintLevel(1)

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
    ts_start = time.time()
    results = wspace.pdf('model').fitTo(data, RooFit.Extended(True), rt.RooFit.Save(), RooFit.NumCPU(8))
    results.Print()
    print("Time for fit: {} s".format(time.time() - ts_start))

    #-----------------------------------------------------    
    # Compute interval based on profile likelihood
    #-----------------------------------------------------
    # suppress some (apparently) innocuous warnings
    #msgservice = rt.RooMsgService.instance()
    #msgservice.setGlobalKillBelow(rt.RooFit.FATAL)
       
    print 'compute interval using profile likelihood'

    # compute an 68.3% limit
    ts_start = time.time()
    plc = rt.RooStats.ProfileLikelihoodCalculator(data, cfg)
    CL  = 0.683
    plc.SetConfidenceLevel(CL)
    plcInterval_1sd= plc.GetInterval()
    lowerLimit_1sd = plcInterval_1sd.LowerLimit(wspace.var('R'))
    upperLimit_1sd = plcInterval_1sd.UpperLimit(wspace.var('R'))

    print '\tPL %4.1f%s CL interval = [%5.5f, %5.5f]' % \
      (100*CL, '%', lowerLimit_1sd, upperLimit_1sd)

    print("Time for profile likelihood 1 sigma: {} s".format(time.time() - ts_start))


    # compute an 95% limit
    ts_start = time.time()
    CL = 0.95
    plc.SetConfidenceLevel(CL)
    plcInterval_2sd = plc.GetInterval()
    lowerLimit_2sd = plcInterval_2sd.LowerLimit(wspace.var('R'))
    upperLimit_2sd = plcInterval_2sd.UpperLimit(wspace.var('R'))

    print '\tPL %4.1f%s CL interval = [%5.5f, %5.5f]' % \
      (100*CL, '%', lowerLimit_2sd, upperLimit_2sd)
      
    print("Time for profile likelihood 2 sigma: {} s".format(time.time() - ts_start))
  
    if plot:
      ts_start = time.time()
      Yat_Xmax_1sd = 0.5*rt.Math.chisquared_quantile(plcInterval_1sd.ConfidenceLevel(),1);
      Yat_Xmax_2sd = 0.5*rt.Math.chisquared_quantile(plcInterval_2sd.ConfidenceLevel(),1);
      plccanvas = plot_likelihood(plcInterval_1sd, Yat_Xmax_1sd, lowerLimit_1sd, upperLimit_1sd, Yat_Xmax_2sd, lowerLimit_2sd, upperLimit_2sd, x_min=x_min, x_max=x_max, x_SM=x_SM, NPoints=NPoints, xlabel=xlabel)
      print("Time for plotting: {} s".format(time.time() - ts_start))
      return plccanvas
    else:
      return None
    
def plot_likelihood(plcInterval_1sd, Yat_Xmax_1sd, min_1sd, max_1sd, Yat_Xmax_2sd, min_2sd, max_2sd, x_min=0.8, x_max=1.2, x_SM=1.0, NPoints=20, xlabel='R'):
    plccanvas = rt.TCanvas('fig_PL', 'PL', 800, 600)
    plcplot = rt.RooStats.LikelihoodIntervalPlot(plcInterval_1sd)      
    plccanvas.cd()
    pad = setup_pad()
    pad.Draw()
    pad.cd()

    plcplot.SetRange(x_min, x_max)
    plcplot.SetMaximum(10)
    plcplot.SetLineColor(0)
    plcplot.SetNPoints(NPoints)
    plcplot.Draw("tf1")

    pad.cd()

    '''
    Yline_cutoff = rt.TLine(x_min,Yat_Xmax,x_max,Yat_Xmax);
    Yline_min = rt.TLine(min_2sd,0.,min_2sd,Yat_Xmax);
    Yline_max = rt.TLine(max_2sd,0.,max_2sd,Yat_Xmax);

    Yline_cutoff.SetLineColor(8)
    Yline_min.SetLineColor(8)
    Yline_max.SetLineColor(8)

    Yline_cutoff.Draw("same");
    Yline_min.Draw("same");
    Yline_max.Draw("same");
    '''

    hist = plcplot.GetPlottedObject()
    frame = rt.TGraph(hist)
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
    frame.SetMinimum(0)
    frame.Draw()
    
    pad.Update()

    # plot 95% limit
    hist_2sd = hist.Clone()
    y_2sd, edge_2sd = root_numpy.hist2array(hist_2sd, return_edges=True)
    edge_2sd = np.array(edge_2sd).flatten()
    x_2sd = np.convolve(edge_2sd, np.ones(2), 'valid') / 2.0
    tofill_2sd = [(x, y) for x, y in zip(x_2sd, y_2sd) if min_2sd < x < max_2sd]
    tg_2sd = rt.TGraph()
    root_numpy.fill_graph(tg_2sd, tofill_2sd)
    tg_2sd.SetPoint(tg_2sd.GetN(), min_2sd, Yat_Xmax_2sd)
    tg_2sd.SetPoint(tg_2sd.GetN(), min_2sd, 0)
    tg_2sd.SetPoint(tg_2sd.GetN(), max_2sd, 0)
    tg_2sd.SetPoint(tg_2sd.GetN(), max_2sd, Yat_Xmax_2sd)
    tg_2sd.Sort()
    tg_2sd.SetFillColor(16)
    tg_2sd.SetFillStyle(1001)
    tg_2sd.Draw("F2 SAME")

    pad.Update()

    # plot 68.3% limit
    hist_1sd = hist.Clone()
    y_1sd, edge_1sd = root_numpy.hist2array(hist_1sd, return_edges=True)
    edge_1sd = np.array(edge_1sd).flatten()
    x_1sd = np.convolve(edge_1sd, np.ones(2), 'valid') / 2.0
    tofill_1sd = [(x, y) for x, y in zip(x_1sd, y_1sd) if min_1sd < x < max_1sd]
    tg_1sd = rt.TGraph()
    root_numpy.fill_graph(tg_1sd, tofill_1sd)
    tg_1sd.SetPoint(tg_1sd.GetN(), min_1sd, Yat_Xmax_1sd)
    tg_1sd.SetPoint(tg_1sd.GetN(), min_1sd, 0)
    tg_1sd.SetPoint(tg_1sd.GetN(), max_1sd, 0)
    tg_1sd.SetPoint(tg_1sd.GetN(), max_1sd, Yat_Xmax_1sd)
    tg_1sd.Sort()
    tg_1sd.SetFillColor(14)
    tg_1sd.SetFillStyle(1001)
    tg_1sd.Draw("F2 SAME")

    pad.Update()

    frame.Draw("same")
    pad.Update()

    uymax = rt.gPad.GetUymax()
    uymin = rt.gPad.GetUymin()
    if x_SM is not None:
      pad.cd()
      l = rt.TLine(x_SM, uymin, x_SM, uymax)
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

