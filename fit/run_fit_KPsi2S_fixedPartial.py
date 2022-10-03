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
    df = pd.DataFrame(f[tree].arrays(branches=branches))
    print('Finished opening file {}...'.format(root_file_name))
    return df

if __name__ == "__main__":
  eleType = 'pf'
  log = 'log_psi2s_bparkPU_v7.3_{}_nonreg_tightWP_otherB.csv'.format(eleType)
  info = defaultdict(dict)

  br_b2kpsi2s = 6.24e-4
  br_b2kstarpsi2s = 5.9e-4 * 2./3.
  br_b2kstarjpsi = 1.27e-3 * 2./3.
  br_psi2s2ee = 7.93e-3 
  br_jpsi2ee = 0.05971

  nparts = range(8)

  # PF-PF nonreg
  info['pf']['inputfile'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_data_mvaCut0.root'
  info['pf']['psi2s_mc'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MCPsi2S.root'.format('marker')
  info['pf']['partial_mc'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MC_kstarpsi2s.root'.format('marker')
  info['pf']['isgn'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCPsi2S.root'
  info['pf']['iKstarPsi2S_BKG'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_kstarpsi2s.root'
  info['pf']['iKstarJpsi_BKG'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_kstarjpsi.root'
  info['pf']['iKJpsiee_BKG'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCres.root'
  info['pf']['n_mc_psi2s'] = 180840.0
  info['pf']['n_mc_partial'] = 443318.0
  info['pf']['ibkg'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_samesign_mvaCut0.root'
  info['pf']['iotherB'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_otherb.root'

  # PF-LP nonreg
  info['mix']['inputfile'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0.root'
  info['mix']['psi2s_mc'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MCPsi2S.root'.format('marker')
  info['mix']['partial_mc'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MC_kstarpsi2s.root'.format('marker')
  info['mix']['isgn'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MCPsi2S.root'
  info['mix']['iKstarPsi2S_BKG'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_kstarpsi2s.root'
  info['mix']['iKstarJpsi_BKG'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_kstarjpsi.root'
  info['mix']['iKJpsiee_BKG'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MCres.root'
  info['mix']['n_mc_psi2s'] = 180840.0
  info['mix']['n_mc_partial'] = 443318.0
  info['mix']['ibkg'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_samesign_mvaCut0.root'
  info['mix']['iotherB'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_otherb.root'

  selection = {}

  selection['jpsi'] = '(Mll > 2.9) and (Mll < 3.2)'
  selection['psi2s'] = '(Mll > 3.55) and (Mll < 3.8)'
  selection['Dmass'] = '(KLmassD0 > 2.0)'

  mc_branches = ['Bmass', 'Mll', 'xgb', 'KLmassD0']

  psi2s_mc_branches = [get_df(info[eleType]['psi2s_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]
  partial_mc_branches = [get_df(info[eleType]['partial_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]

  if eleType == 'pf':
    mvaCut = np.array([8.3,])      # 8.3
  else: 
    mvaCut = np.array([6.0,])      # 8.6

  for cut in mvaCut:
    eff_sig_bdt = np.mean([float(psi2s_mc_branches[i].query(' and '.join([selection['psi2s'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_psi2s'] for i in nparts])
    eff_partial_bdt = np.mean([float(partial_mc_branches[i].query(' and '.join([selection['psi2s'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_partial'] for i in nparts])

    frac_ratio = (eff_partial_bdt / eff_sig_bdt) * (br_b2kstarpsi2s / br_b2kpsi2s)
    print(frac_ratio, eff_partial_bdt, eff_sig_bdt)

    if eleType == 'pf':
      com = 'python KPsi2S_roofit_plb_modified_kde_fixedPartial.py -i {} -o psi2s_fixedPartial_{} --isgn={} --iotherB={} --iKstarPsi2S_BKG={} --iKstarJpsi_BKG={} --iKJpsiee_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarpsi2s,bkg_kstarjpsi,bkg_kjpsi_ee,bkg_comb,bkg_otherb" --fit_primitive --partial_ratio={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iotherB'], info[eleType]['iKstarPsi2S_BKG'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['iKJpsiee_BKG'], info[eleType]['ibkg'], frac_ratio, cut, log)

    else:
      com = 'python KPsi2S_roofit_plb_modified_kde_fixedPartial.py -i {} -o psi2s_fixedPartial_{} --isgn={} --iotherB={} --iKstarPsi2S_BKG={} --iKstarJpsi_BKG={} --iKJpsiee_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarpsi2s,bkg_kstarjpsi,bkg_kjpsi_ee,bkg_comb,bkg_otherb" --fit_primitive --partial_ratio={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iotherB'], info[eleType]['iKstarPsi2S_BKG'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['iKJpsiee_BKG'], info[eleType]['ibkg'], frac_ratio, cut, log)



    os.system(com)



