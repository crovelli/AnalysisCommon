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
#import PyPDF2
#import makePlot_fitPeak_unbinned as fit_unbinned
import os, sys, copy
import xgboost as xgb

LOWQ2_LOW  = 1.05
LOWQ2_UP   = 2.45
JPSI_LOW   = 2.9
JPSI_UP    = 3.2
PSI2S_LOW  = 3.55
PSI2S_UP   = 3.8
HIGHQ2_LOW = 4.0
HIGHQ2_UP  = 4.8

def get_df(root_file_name, tree='mytreefit', branches=['*']):
    print('Opening file {}...'.format(root_file_name))
    f = uproot.open(root_file_name)
    if len(f.allkeys()) == 0:
        return pd.DataFrame()
    print('Not an null file')
    df = pd.DataFrame(f[tree].arrays(branches=branches))
    print('Finished opening file {}...'.format(root_file_name))
    return df

if __name__ == "__main__":
  eleType = 'mix'        # chiara
  log = 'log_jpsi_bparkPU_v7.3_{}.csv'.format(eleType)
  info = defaultdict(dict)

  # chiara - vedi questi
  br_b2jpsi = 1.02e-3
  br_b2kstarjpsi = 1.27e-3 * 2./3.

  nparts = range(8)

  info['pf']['inputfile'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_data_mvaCut0.root'
  info['pf']['jpsi_mc'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MCres.root'.format('marker')
  info['pf']['partial_mc'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MC_kstarjpsi.root'.format('marker')
  info['pf']['isgn'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCres.root'
  #####info['pf']['ibkg'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_samesign_mvaCut0.root'
  info['pf']['ibkg'] = 'ottoPFPF/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_samesign_mvaCut0.root'        ### chiara - passare a non regressed
  info['pf']['iKstarJpsi_BKG'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_kstarjpsi.root'
  info['pf']['n_mc_jpsi'] = 563421.0
  info['pf']['n_mc_partial'] = 373882.0

  info['mix']['inputfile'] = 'ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0.root'
  info['mix']['jpsi_mc'] = 'ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MCres.root'.format('marker')
  info['mix']['partial_mc'] = 'ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MC_kstarjpsi.root'.format('marker')
  info['mix']['isgn'] = 'ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MCres.root'
  info['mix']['iKstarJpsi_BKG'] = 'ottoPFLPnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_kstarjpsi.root'
  info['mix']['n_mc_jpsi'] = 563421.0
  info['mix']['n_mc_partial'] = 373882.0

  selection = {}

  selection['jpsi']  = '(Mll > 2.9) and (Mll < 3.2)'
  selection['psi2s'] = '(Mll > 3.55) and (Mll < 3.8)'
  selection['Dmass'] = '(KLmassD0 > 2.0)'

  mc_branches = ['Bmass', 'Mll', 'xgb', 'KLmassD0']

  jpsi_mc_branches = [get_df(info[eleType]['jpsi_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]
  partial_mc_branches = [get_df(info[eleType]['partial_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]


  if eleType == 'pf':
    mvaCut = np.linspace(5.0, 5.0, 1)
  else:
    mvaCut = np.linspace(5.0, 5.0, 1)

  for cut in mvaCut:
    eff_sig_bdt = np.mean([float(jpsi_mc_branches[i].query(' and '.join([selection['jpsi'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_jpsi'] for i in nparts])
    eff_partial_bdt = np.mean([float(partial_mc_branches[i].query(' and '.join([selection['jpsi'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_partial'] for i in nparts])

    frac_ratio = (eff_partial_bdt / eff_sig_bdt) * (br_b2kstarjpsi / br_b2jpsi)

    if eleType == 'pf':
      com = 'python KJpsi_roofit_plb_modified_kde_fixedPartial.py -i {} -o jpsi_fixedPartial_{} --isgn={} --iKstarJpsi_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarjpsi,bkg_comb" --fit_primitive --partial_ratio={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['ibkg'], frac_ratio, cut, log)

    else:
      com = 'python KJpsi_roofit_plb_modified_kde_fixedPartial.py -i {} -o jpsi_fixedPartial_{} --isgn={} --iKstarJpsi_BKG={} --sel_primitive="sgn,bkg_kstarjpsi" --fit_primitive --partial_ratio={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iKstarJpsi_BKG'], frac_ratio, cut, log)

    os.system(com)



