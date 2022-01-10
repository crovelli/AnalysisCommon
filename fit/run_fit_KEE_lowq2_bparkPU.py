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
  log = 'log_kee_bparkPU_v7.2_{}.csv'.format(eleType)
  info = defaultdict(dict)

  nparts = range(8)

  info['pf']['inputfile'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_data_mvaCut0.root'
  info['pf']['nonresonant_mc'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_{}_MC.root'.format('marker')
  info['pf']['jpsi_mc'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_{}_MCres.root'.format('marker')
  info['pf']['isgn'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0_MC.root'
  info['pf']['ikjpsi'] = 'data_PFe_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.2_0_MCres.root'
  info['pf']['jpsi_mva_wp'] = 4.4
  info['pf']['n_data_jpsi'] = 9424.52633491936
  info['pf']['n_mc_jpsi'] = 563421.0
  info['pf']['n_mc_lowq2'] = 617615.0

  #info['pf']['inputfile'] = 'data_PFe_v7.2_run3_elePt5/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_elePt5_v7.2_data.root'
  #info['pf']['isgn'] = 'data_PFe_v7.2_run3_elePt5/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_elePt5_v7.2_0_MC.root'
  #info['pf']['ikjpsi'] = 'data_PFe_v7.2_run3_elePt5/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_elePt5_v7.2_0_MCres.root'

  #info['pf']['inputfile'] = 'data_PFe_v7.2_run3_elePt10/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_elePt10_v7.2_data.root'
  #info['pf']['isgn'] = 'data_PFe_v7.2_run3_elePt10/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_elePt10_v7.2_0_MC.root'
  #info['pf']['ikjpsi'] = 'data_PFe_v7.2_run3_elePt10/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_elePt10_v7.2_0_MCres.root'

  info['mix']['inputfile'] = 'data_LowPtPF_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_data_mvaCut3.root'
  info['mix']['nonresonant_mc'] = 'data_LowPtPF_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_{}_MC.root'.format('marker')
  info['mix']['jpsi_mc'] = 'data_LowPtPF_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_{}_MCres.root'.format('marker')
  info['mix']['isgn'] = 'data_LowPtPF_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_0_MC.root'
  info['mix']['ikjpsi'] = 'data_LowPtPF_v7.2/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.2_0_MCres.root'
  info['mix']['jpsi_mva_wp'] = 5.2
  info['mix']['n_data_jpsi'] = 6710.791959557816
  info['mix']['n_mc_jpsi'] = 563421.0
  info['mix']['n_mc_lowq2'] = 617615.0

  num_mctoys = None #5000

  rk_electron = 7.20e-3

  jpsi_mva_wp = info[eleType]['jpsi_mva_wp']
  n_data_jpsi = info[eleType]['n_data_jpsi']

  selection = {}

  selection['jpsi'] = '(Mll > @JPSI_LOW) and (Mll < @JPSI_UP)'
  selection['psi2s'] = '(Mll > @PSI2S_LOW) and (Mll < @PSI2S_UP)'
  selection['lowq2'] = '(Mll > @LOWQ2_LOW) and (Mll < @LOWQ2_UP)'
  selection['highq2'] = '(Mll > @PSI2S_LOW) and (Mll < @PSI2S_UP)'

  selection['Dmass'] = '(KLmassD0 > 2.0)'
  #selection['elePt'] = '(L1pt > 10) and (L2pt > 10)'

  mc_branches = ['Bmass', 'Mll', 'KLmassD0', 'xgb'] #+ ['L1pt', 'L2pt'] 

  jpsi_mc_branches = [get_df(info[eleType]['jpsi_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]
  nonresonant_mc_branches = [get_df(info[eleType]['nonresonant_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]

  eff = {}
  # efficieny of preselection and di-electron type (after nanoAOD preselection)
  eff['jpsi_presel'] = float(jpsi_mc_branches[0].shape[0]) / info[eleType]['n_mc_jpsi']
  eff['lowq2_presel'] = float(nonresonant_mc_branches[0].shape[0]) / info[eleType]['n_mc_lowq2']
  # efficiency of q2
  eff['jpsi_q2'] = float(jpsi_mc_branches[0].query(selection['jpsi']).shape[0]) / float(jpsi_mc_branches[0].shape[0])
  eff['lowq2_q2'] = float(nonresonant_mc_branches[0].query(selection['lowq2']).shape[0]) / float(nonresonant_mc_branches[0].shape[0])
  # efficiency of bdt of jpsi
  eff['jpsi_bdt'] = np.mean([float(jpsi_mc_branches[i].query(' and '.join([selection['jpsi'], '(xgb > {})'.format(jpsi_mva_wp)])).shape[0]) / float(jpsi_mc_branches[i].query(selection['jpsi']).shape[0]) for i in nparts])

  if eleType == 'pf':
    mvaCut = np.linspace(7.0, 9.0, 20)
  else:
    mvaCut = np.linspace(8.0, 10.0, 20)

  for cut in mvaCut:
    eff_lowq2_bdt = np.mean([float(nonresonant_mc_branches[i].query(' and '.join([selection['lowq2'], '(xgb > @cut)', selection['Dmass']])).shape[0]) / float(nonresonant_mc_branches[i].query(selection['lowq2']).shape[0]) for i in nparts])

    expected_signal = rk_electron * n_data_jpsi * eff['lowq2_presel'] * eff['lowq2_q2'] * eff_lowq2_bdt / (eff['jpsi_presel'] * eff['jpsi_q2'] * eff['jpsi_bdt'])

    if num_mctoys is not None:
      com = 'python KEE_lowq2_roofit_plb_pull_modified.py -i {0} -o lowq2_{1} --isgn={2} --ikjpsi={3} --sel_primitive="sgn,bkg_kjpsi" --fit_primitive --mvacut={4} --set_expected_sgn={5} --number_of_mctoys={6} --log={7}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['ikjpsi'], cut, expected_signal, num_mctoys, log)
    else:
      com = 'python KEE_lowq2_roofit_plb_pull_modified.py -i {0} -o lowq2_{1} --isgn={2} --ikjpsi={3} --sel_primitive="sgn,bkg_kjpsi" --fit_primitive --mvacut={4} --set_expected_sgn={5} --log={7}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['ikjpsi'], cut, expected_signal, num_mctoys, log)

    os.system(com)
    print('cut: {} \n \t expected signal: {} \n \t lowq2 bdt eff: {}'.format(cut, expected_signal, eff_lowq2_bdt))

  for key, value in eff.items():
    print('{}: {}'.format(key, value))




