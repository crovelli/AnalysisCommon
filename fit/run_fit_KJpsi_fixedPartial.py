import uproot
import pandas as pd
import numpy as np
from collections import OrderedDict, defaultdict
from scipy import interp
from rootpy.io import root_open
from rootpy.plotting import Hist
from root_numpy import fill_hist, array2root, array2tree
from root_pandas import to_root
import itertools
import PyPDF2
#import makePlot_fitPeak_unbinned as fit_unbinned
import os, sys, copy
import xgboost as xgb

LOWQ2_LOW = 1.05
LOWQ2_UP = 2.45
JPSI_LOW = 2.8
JPSI_UP = 3.25
PSI2S_LOW = 3.45
PSI2S_UP = 3.8
HIGHQ2_LOW = 4.0
HIGHQ2_UP = 4.87

def get_df(root_file_name, tree='mytreefit', branches=['*']):
    print('Opening file {}...'.format(root_file_name))
    f = uproot.open(root_file_name)
    if len(f.allkeys()) == 0:
        return pd.DataFrame()
    print('Not an null file')
    #df = uproot.open(root_file_name)["tree"].pandas.df()
    #df = pd.DataFrame(uproot.open(root_file_name)["tree"].arrays(namedecode="utf-8"))
    df = pd.DataFrame(f[tree].arrays(branches=branches))
    print('Finished opening file {}...'.format(root_file_name))
    return df

if __name__ == "__main__":
  eleType = 'pf'
  log = 'log_jpsi_bparkPU_v7.2_{}_kstarplus.csv'.format(eleType)
  info = defaultdict(dict)

  br_b2jpsi = 1.02e-3
  br_b2kstarjpsi = 1.27e-3 * 2./3.

  nparts = range(8)

  info['pf']['inputfile'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_data_mvaCut0.root'
  info['pf']['jpsi_mc'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_{}_MCres.root'.format('marker')
  info['pf']['partial_mc'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_{}_MC_kstarjpsi_combined.root'.format('marker')
  info['pf']['isgn'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0_MCres.root'
  info['pf']['ibkg'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_samesign_mvaCut0.root'
  info['pf']['iKstarJpsi_BKG'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0_MC_kstarjpsi_combined.root'
  info['pf']['iKstarPlusJpsi_BKG'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0_MC_kstarplusjpsi_kee.root'
  info['pf']['n_mc_jpsi'] = 563421.0
  info['pf']['n_mc_partial'] = 1023330.0

  selection = {}

  selection['jpsi'] = '(Mll > @JPSI_LOW) and (Mll < @JPSI_UP)'
  selection['psi2s'] = '(Mll > @PSI2S_LOW) and (Mll < @PSI2S_UP)'

  mc_branches = ['Bmass', 'Mll', 'xgb']

  jpsi_mc_branches = [get_df(info[eleType]['jpsi_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]
  partial_mc_branches = [get_df(info[eleType]['partial_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]


  if eleType == 'pf':
    mvaCut = np.linspace(3.0, 5.0, 11)
  else:
    mvaCut = np.linspace(8.0, 10.0, 20)

  for cut in mvaCut:
    eff_sig_bdt = np.mean([float(jpsi_mc_branches[i].query(' and '.join([selection['jpsi'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_jpsi'] for i in nparts])
    eff_partial_bdt = np.mean([float(partial_mc_branches[i].query(' and '.join([selection['jpsi'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_partial'] for i in nparts])

    frac_ratio = (eff_partial_bdt / eff_sig_bdt) * (br_b2kstarjpsi / br_b2jpsi)

    com = 'python KJpsi_roofit_plb_modified_kde_fixedPartial.py -i {} -o jpsi_fixedPartial_{} --isgn={} --iKstarJpsi_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarjpsi,bkg_comb" --fit_primitive --partial_ratio={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['ibkg'], frac_ratio, cut, log)
    #com = 'python KJpsi_roofit_plb_modified_kde_kstarPlus.py -i {} -o jpsi_fixedPartial_{} --isgn={} --iKstarJpsi_BKG={} --iKstarPlusJpsi_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarjpsi,bkg_kstarplusjpsi,bkg_comb" --fit_primitive --partial_ratio={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['iKstarPlusJpsi_BKG'], info[eleType]['ibkg'], frac_ratio, cut, log)

    os.system(com)



