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
  eleType = 'pf'           # chiara
  log = 'log_kee_bparkPU_v7.3_{}.csv'.format(eleType)
  info = defaultdict(dict)

  nparts = range(8)

  info['pf']['inputfile'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_data_mvaCut0.root'
  info['pf']['nonresonant_mc'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MC.root'.format('marker')
  info['pf']['jpsi_mc'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MCres.root'.format('marker')
  info['pf']['isgn'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC.root'
  info['pf']['ikjpsi'] = 'ottoPFPFnonReg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCres.root'
  info['pf']['ikstaree'] = 'ottoPFPF/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_0_MC_kstaree.root'         ## chiara, passare a non regressed
  info['pf']['jpsi_mva_wp'] = 6.0
  info['pf']['n_data_jpsi'] = 7069.5
  info['pf']['n_mc_jpsi'] = 563421.0
  info['pf']['n_mc_lowq2'] = 617615.0

  info['mix']['inputfile'] = '../data/data_LowPtPF_v7.3/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_data_mvaCut0.root'
  info['mix']['nonresonant_mc'] = '../data/data_LowPtPF_v7.3/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_{}_MC.root'.format('marker')
  info['mix']['jpsi_mc'] = '../data/data_LowPtPF_v7.3/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_{}_MCres.root'.format('marker')
  info['mix']['isgn'] = '../data/data_LowPtPF_v7.3/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_0_MC.root'
  info['mix']['ikjpsi'] = '../data/data_LowPtPF_v7.3/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_0_MCres.root'
  info['mix']['ikstaree'] = '../data/data_LowPtPF_v7.3/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_0_MC_kstaree.root'
  info['mix']['jpsi_mva_wp'] = 5.0
  info['mix']['n_data_jpsi'] = 5533.45
  info['mix']['n_mc_jpsi'] = 563421.0
  info['mix']['n_mc_lowq2'] = 617615.0

  num_mctoys = None #5000

  rk_electron = 7.20e-3

  jpsi_mva_wp = info[eleType]['jpsi_mva_wp']
  n_data_jpsi = info[eleType]['n_data_jpsi']

  selection = {}

  selection['jpsi'] = '(Mll > 2.9) and (Mll < 3.2)'
  selection['psi2s'] = '(Mll > 3.55) and (Mll < 3.8)'
  selection['lowq2'] = '(Mll > 1.05) and (Mll < 2.45)'
  selection['highq2'] = '(Mll > @PSI2S_LOW) and (Mll < @PSI2S_UP)'

  selection['Dmass'] = '(KLmassD0 > 2.0)'

  mc_branches = ['Bmass', 'Mll', 'KLmassD0', 'xgb'] 

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
  eff['jpsi_bdt'] = np.mean([float(jpsi_mc_branches[i].query(' and '.join([selection['jpsi'], selection['Dmass'], '(xgb > {})'.format(jpsi_mva_wp)])).shape[0]) / float(jpsi_mc_branches[i].query(selection['jpsi']).shape[0]) for i in nparts])

  if eleType == 'pf':
    mvaCut = np.linspace(8.44, 8.44, 1)
  else:
    mvaCut = np.linspace(8.44, 8.44, 1)

  for cut in mvaCut:
    eff_lowq2_bdt = np.mean([float(nonresonant_mc_branches[i].query(' and '.join([selection['lowq2'], '(xgb > @cut)', selection['Dmass']])).shape[0]) / float(nonresonant_mc_branches[i].query(selection['lowq2']).shape[0]) for i in nparts])

    expected_signal = rk_electron * n_data_jpsi * eff['lowq2_presel'] * eff['lowq2_q2'] * eff_lowq2_bdt / (eff['jpsi_presel'] * eff['jpsi_q2'] * eff['jpsi_bdt'])

    if num_mctoys is not None:
      com = 'python KEE_lowq2_roofit_plb_pull_modified.py -i {0} -o lowq2_{1} --isgn={2} --ikjpsi={3} --sel_primitive="sgn,bkg_kjpsi" --fit_primitive --mvacut={4} --set_expected_sgn={5} --number_of_mctoys={6} --log={7}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['ikjpsi'], cut, expected_signal, num_mctoys, log)
    else:
      com = 'python KEE_lowq2_roofit_plb_pull_modified_kde.py -i {0} -o lowq2_{1} --isgn={2} --ikjpsi={3} --ikstaree={4} --sel_primitive="sgn,bkg_kjpsi,bkg_kstaree" --fit_primitive --mvacut={5} --set_expected_sgn={6} --log={8}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['ikjpsi'], info[eleType]['ikstaree'], cut, expected_signal, num_mctoys, log)

    os.system(com)
    print('cut: {} \n \t expected signal: {} \n \t lowq2 bdt eff: {}'.format(cut, expected_signal, eff_lowq2_bdt))

  for key, value in eff.items():
    print('{}: {}'.format(key, value))




