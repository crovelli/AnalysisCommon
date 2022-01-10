import ROOT as rt
from math import sqrt


def CMS_lumi():
    mark = rt.TLatex()
    mark.SetNDC()
    lumistamp = '2018 (13 TeV)'
    fontScale = 1.0
    cmsTextSize = 0.042 * fontScale * 1.25
    extraOverCmsTextSize  = 0.76
    extraTextSize = extraOverCmsTextSize*cmsTextSize

    mark.SetTextAlign(11)
    mark.SetTextSize(cmsTextSize)
    mark.SetTextFont(61)
    mark.DrawLatex(rt.gPad.GetLeftMargin(), 1 - (rt.gPad.GetTopMargin() - 0.017), "CMS")
    mark.SetTextSize(0.042 * fontScale)
    mark.SetTextFont(52)
    mark.DrawLatex(rt.gPad.GetLeftMargin() + 0.09, 1 - (rt.gPad.GetTopMargin() - 0.017),  "Preliminary")
    mark.SetTextSize(extraTextSize)
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - rt.gPad.GetRightMargin(), 1 - ( rt.gPad.GetTopMargin() - 0.017), lumistamp)
    return mark


def canvas_create(xframe,xmin,xmax,nbin,xtitle,boundYto0=True):
    c2 = rt.TCanvas('fig_binnedFit', 'fit', 800, 600)
    c2.SetGrid()
    rt.gPad.SetLeftMargin(0.12)
    rt.gPad.SetRightMargin(0.05)
    rt.gPad.SetBottomMargin(0.15)
    #xframe.GetYaxis().SetTitle("Events / {0:.0f} MeV".format((xmax - xmin)/nbin*1000.))
    xframe.GetXaxis().SetTitle(xtitle)
    xframe.SetStats(0)
    if (boundYto0): xframe.SetMinimum(0)
    xframe.GetYaxis().SetLabelSize(0.045)
    xframe.GetYaxis().SetLabelOffset(0.007)
    xframe.GetYaxis().SetTitleSize(0.06)
    xframe.GetYaxis().SetTitleOffset(0.90)
    xframe.GetXaxis().SetLabelSize(0.045)
    xframe.GetXaxis().SetLabelOffset(0.007)
    xframe.GetXaxis().SetTitleSize(0.06)
    xframe.GetXaxis().SetTitleOffset(0.95)
    xframe.SetTitle("")
    xframe.Draw()
    return c2

def pt_create(mva,nsig,nsigError,nbkg):
  pt = rt.TPaveText(0.72,0.37,0.92,0.63,"brNDC")
  pt.SetFillColor(0)
  pt.SetBorderSize(1)
  pt.SetTextFont(42);
  pt.SetTextSize(0.04);
  pt.SetTextAlign(12)
  pt.AddText("MVA cut: {0}".format(mva))
  pt.AddText("S: {0:.1f}#pm{1:.1f}".format(nsig,nsigError))
  pt.AddText("B: {0:.1f}".format(nbkg))
  pt.AddText("S/#sqrt{{S+B}}: {0:.2f}".format(nsig/sqrt(nsig + nbkg)))
  return pt
