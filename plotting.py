import uproot
import pandas as pd
import numpy as np
from collections import OrderedDict
from scipy import interp
from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D
from root_numpy import fill_hist, array2root, array2tree
from root_pandas import to_root
import ROOT
from ROOT import RooFit
import itertools
import PyPDF2
import os, sys, copy
from helper import *

#ROOT.gErrorIgnoreLevel=ROOT.kError
#ROOT.RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

'''
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
'''

def pdf_combine(pdf_list, outputfile):
    merger = PyPDF2.PdfFileMerger()
    for pdf in pdf_list:
        merger.append(pdf)
    merger.write(outputfile)

def pdf_sidebyside(out, inputfile):
    output = PyPDF2.PdfFileWriter()
    reader = [PyPDF2.PdfFileReader(file(in1, "rb")) for in1 in inputfile]
    m = min([in1.getNumPages() for in1 in reader])
    print("common pages",m)
    for i in range(0,m):
        print "adding page common",i
        p = [in1.getPage(i) for in1 in reader]
        nPages = len(p)
        p1 = p[0]
        offset_x = 0.0
        for i, p2 in enumerate(p[1:]):
            offset_y = -(i+1)*p1.cropBox[1] + (i+1)*p1.cropBox[3]
            p1.mergeTranslatedPage(p2, offset_x, offset_y, expand=True)
        bounding_box = copy.deepcopy(p1.cropBox)
        p1.trimBox.lowerLeft = (bounding_box[0], bounding_box[1])
        p1.trimBox.upperRight = (bounding_box[2], bounding_box[1] + nPages*(bounding_box[3] - bounding_box[1]))
        p1.cropBox.lowerLeft = (bounding_box[0], bounding_box[1])
        p1.cropBox.upperRight = (bounding_box[2], bounding_box[1] + nPages*(bounding_box[3] - bounding_box[1]))
        output.addPage(p1)
    outputStream = file(out, "wb")
    output.write(outputStream)
    outputStream.close()

def setup_pad():
    pad = ROOT.TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.12)
    pad.SetLeftMargin(0.11)
    pad.SetRightMargin(0.06)
    return pad

def CMS_lumi(isMC=False):
    mark = ROOT.TLatex()
    mark.SetNDC()
    lumistamp = '2018 (13 TeV)'
    fontScale = 1.0
    cmsTextSize = 0.042 * fontScale * 1.25
    extraOverCmsTextSize  = 0.76
    extraTextSize = extraOverCmsTextSize*cmsTextSize
    mark.SetTextAlign(11)
    mark.SetTextSize(cmsTextSize)
    mark.SetTextFont(61)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), "CMS")
    mark.SetTextSize(0.042 * fontScale)
    mark.SetTextFont(52)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.09, 1 - (ROOT.gPad.GetTopMargin() - 0.017), "Simulation Preliminary" if isMC else "Preliminary")
    mark.SetTextSize(extraTextSize)
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.017), lumistamp)


def draw_hist(histo, histo_name, x_label, y_label, draw_option='E'):
    #histo.SetTitle(histo_name)
    histo.GetYaxis().SetTitle(y_label)
    histo.GetXaxis().SetTitle(x_label)
    histo.SetTitleFont(42)
    histo.SetTitleSize(0.05)
    histo.GetYaxis().SetTitleOffset(0.9)
    histo.GetYaxis().SetTitleFont(42)
    histo.GetYaxis().SetTitleSize(0.04)
    histo.GetYaxis().SetLabelSize(0.04)
    histo.GetYaxis().SetLabelFont(42)
    histo.GetXaxis().SetTitleOffset(0.9)
    histo.GetXaxis().SetTitleFont(42)
    histo.GetXaxis().SetTitleSize(0.04)
    histo.GetXaxis().SetLabelSize(0.04)
    histo.GetXaxis().SetLabelFont(42)
    #histo.Draw(draw_option)
    if (histo.GetEntries() > 0):
      #histo.DrawNormalized(draw_option)
      histo.Draw(draw_option)
    else:
      histo.Draw(draw_option)
    

def plot_hist(selected_branches, name, color, x_label, hist_bins, outputfile, draw_option='E', logScale=False, xlogScale=False, isMC=False, overflow=False, underflow=False, weights=None):
    if not (len(selected_branches) == len(name) == len(color)): 
      print('lengths do not match')
      return

    hist = [Hist(hist_bins['nbins'], hist_bins['xmin'], hist_bins['xmax'], name=n, title='', type='F') for n in name]
    _ = [fill_hist(h, branch) for h, branch in zip(hist, selected_branches)] if weights is None else [fill_hist(h, branch, weights=w) for h, branch, w in zip(hist, selected_branches, weights)]

    canvas_name = "c_{}".format(outputfile)
    ylabel = 'Events' if len(hist) < 2 else 'Fraction'

    c = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    c.cd()
    pad = setup_pad()
    pad.Draw()
    pad.cd()

    if logScale: pad.SetLogy()
    if xlogScale: pad.SetLogx()

    #l1 = ROOT.TLegend(0.6,0.8,0.92,0.9)
    #l1 = ROOT.TLegend(0.7,0.7,0.92,0.9)
    l1 = ROOT.TLegend(0.6,0.7,0.92,0.9)
    #l1 = ROOT.TLegend(0.6,0.2,0.92,0.4)
    l1.SetTextFont(42)
    l1.SetTextSize(0.03)
    
    for h in hist:
      x_min = 1
      x_max = h.GetNbinsX()
      if overflow:
        x_max += 1
      if underflow:
        x_min = 0
      h.GetXaxis().SetRange(x_min, x_max)

    y_max = hist[0].GetMaximum()
    first_hist = True
    for h, n, co in zip(hist, name, color):
      h.SetMarkerStyle(0)
      h.SetFillStyle(0)
      h.SetLineColor(co)
      if first_hist:
        if len(hist) > 1: 
          h.SetMaximum(1.3*y_max)
          #h.SetMaximum(2.0*y_max)
          #h.SetMaximum(1)
        draw_hist(h, n, x_label, ylabel, draw_option=draw_option)
        first_hist = False
      else:
        draw_hist(h, n, x_label, ylabel, draw_option='{} SAME'.format(draw_option))
      l1.AddEntry(h, n)

    if len(hist) > 1: l1.Draw("same")

    pad.cd()
    CMS_lumi(isMC)
    c.cd()
    c.Update()
    c.SaveAs(outputfile)
    c.Close()

def plot_hist_from_hist(hist, name, color, x_label, y_label, outputfile, draw_option='E', logScale=False, xlogScale=False, isMC=False, overflow=False, underflow=False, weights=None):
    if not (len(hist) == len(name) == len(color)): 
      print('lengths do not match')
      return

    canvas_name = "c_{}".format(outputfile)

    c = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    c.cd()
    pad = setup_pad()
    pad.Draw()
    pad.cd()

    if logScale: pad.SetLogy()
    if xlogScale: pad.SetLogx()

    #l1 = ROOT.TLegend(0.6,0.8,0.92,0.9)
    #l1 = ROOT.TLegend(0.7,0.7,0.92,0.9)
    l1 = ROOT.TLegend(0.6,0.7,0.92,0.9)
    #l1 = ROOT.TLegend(0.6,0.2,0.92,0.4)
    l1.SetTextFont(42)
    l1.SetTextSize(0.03)
    
    for h in hist:
      x_min = 1
      x_max = h.GetNbinsX()
      if overflow:
        x_max += 1
      if underflow:
        x_min = 0
      h.GetXaxis().SetRange(x_min, x_max)

    y_max = hist[0].GetMaximum()
    first_hist = True
    for h, n, co in zip(hist, name, color):
      h.SetMarkerStyle(0)
      h.SetFillStyle(0)
      h.SetLineColor(co)
      if first_hist:
        if len(hist) > 1: 
          h.SetMaximum(1.3*y_max)
          #h.SetMaximum(2.0*y_max)
          #h.SetMaximum(1)
        draw_hist(h, n, x_label, y_label, draw_option=draw_option)
        first_hist = False
      else:
        draw_hist(h, n, x_label, y_label, draw_option='{} SAME'.format(draw_option))
      l1.AddEntry(h, n)

    if len(hist) > 1: l1.Draw("same")

    pad.cd()
    CMS_lumi(isMC)
    c.cd()
    c.Update()
    c.SaveAs(outputfile)
    c.Close()

def plot_hist2d(xvar_np, yvar_np, x_label, y_label, hist_bins, outputfile, logScale=False, isMC=False):
    hist =  Hist2D(hist_bins['nbinx'], hist_bins['xmin'], hist_bins['xmax'], hist_bins['nbiny'], hist_bins['ymin'], hist_bins['ymax'], title='', type='F')
    fill_hist(hist, np.vstack((xvar_np[np.isfinite(xvar_np)], yvar_np[np.isfinite(yvar_np)])).T)

    ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)
    c = ROOT.TCanvas('c1', 'c1', 800, 600)
    pad = setup_pad()
    pad.SetLeftMargin(0.08)
    pad.SetRightMargin(0.11)
    pad.Draw()
    pad.cd()

    hist.SetContour(100)
    hist.GetYaxis().SetTitle(y_label)
    hist.GetXaxis().SetTitle(x_label)
    hist.SetTitleFont(42)
    hist.SetTitleSize(0.05)
    hist.GetYaxis().SetTitleOffset(0.9)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetTitleSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetXaxis().SetTitleOffset(0.9)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleSize(0.04)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetXaxis().SetLabelFont(42)
    if logScale:  ROOT.gPad.SetLogz(1)
    hist.Draw("colz")

    pad.cd()
    CMS_lumi(isMC)
    c.cd()
    c.Update()
    c.SaveAs(outputfile)
    c.Close()

