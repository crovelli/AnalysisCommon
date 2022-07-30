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
  eleType = 'pf'        # chiara
  log = 'log_jpsi_bparkPU_v7.3_{}.csv'.format(eleType)
  info = defaultdict(dict)

  ## 
  br_b2kjpsi = 1.02e-3
  br_b2kstarjpsi = 1.27e-3 * 2./3.
  br_b2kstarplusjpsi = 1.43e-3

  nparts = range(8)

  # PF-PF
  #info['pf']['inputfile'] = '../splots/outputFiles_PFPF_forEffOfBDTonly/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_data_mvaCut0____withWeights.root'
  info['pf']['inputfile'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_data_mvaCut0.root'
  info['pf']['jpsi_mc'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MCres.root'.format('marker')
  info['pf']['partial_mc'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MC_kstarjpsi.root'.format('marker')
  info['pf']['kstarplus_mc'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_{}_MC_kstarplusjpsi.root'.format('marker')
  info['pf']['isgn'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MCres.root'
  info['pf']['iKstarJpsi_BKG'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_kstarjpsi.root'
  info['pf']['iKstarPlusJpsi_BKG'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_kstarplusjpsi.root'
  info['pf']['iotherB'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_0_MC_otherb.root'
  # chiara: questi sono i numeri di eventi nel dataset, prima della preselezione dei nanoAOD
  # secondo me ci andrebbero questi divisi per 8 blocchi
  #info['pf']['n_mc_jpsi'] = 649505.38          ## 5196043./8.
  #info['pf']['n_mc_partial'] = 1171642.8       ## 9373142./8.                
  #info['pf']['n_mc_kstarplus'] = 1163462.6     ## 9307701./8.
  # chiara: questi sono i numeri che usava otto, da capire cosa fossero 
  info['pf']['n_mc_jpsi'] = 211912.0                   
  info['pf']['n_mc_partial'] = 379086.0                
  info['pf']['n_mc_kstarplus'] = 314672.0
  info['pf']['ibkg'] = 'otto_PFPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_PFe_v7.3_nonreg_ottoCut_samesign_mvaCut0.root'

  # PF-LP
  info['mix']['inputfile'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0.root'
  #info['mix']['inputfile'] = '../splots/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withEle1PtWeight.root'
  #info['mix']['inputfile'] = '../splots/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withEle2PtWeight.root'
  #info['mix']['inputfile'] = '../splots/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_data_mvaCut0_withElesDrWeight.root'

  info['mix']['jpsi_mc'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MCres.root'.format('marker')
  info['mix']['partial_mc'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MC_kstarjpsi.root'.format('marker')
  info['mix']['kstarplus_mc'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_{}_MC_kstarplusjpsi.root'.format('marker')
  info['mix']['isgn'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MCres.root'
  info['mix']['iKstarJpsi_BKG'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_kstarjpsi.root'
  info['mix']['iKstarPlusJpsi_BKG'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_kstarplusjpsi.root'
  info['mix']['iotherB'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_0_MC_otherb.root'
  info['mix']['n_mc_jpsi'] = 211912.0
  info['mix']['n_mc_partial'] = 379086.0
  info['mix']['n_mc_kstarplus'] = 314672.0
  info['mix']['ibkg'] = 'otto_LowPtPF_v7.3_nonreg/forMeas_xgbmodel_kee_12B_kee_correct_pu_Depth17_LowPtPF_v7.3_nonreg_samesign_mvaCut0.root'

  selection = {}

  selection['jpsi']  = '(Mll > 2.9) and (Mll < 3.2)'
  selection['psi2s'] = '(Mll > 3.55) and (Mll < 3.8)'
  selection['Dmass'] = '(KLmassD0 > 2.0)'                ## chiara
  ##selection['Dmass'] = '(KLmassD0 > -10000000000)'     # per calcolo efficienze

  mc_branches = ['Bmass', 'Mll', 'xgb', 'KLmassD0']

  jpsi_mc_branches = [get_df(info[eleType]['jpsi_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]
  partial_mc_branches = [get_df(info[eleType]['partial_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]
  kstarplus_mc_branches = [get_df(info[eleType]['kstarplus_mc'].replace('marker', str(i)), branches=mc_branches) for i in nparts]

  if eleType == 'pf':
    mvaCut = np.linspace(8.3, 8.3, 1)               # otto: 8.3
    #mvaCut = np.linspace(0., 0., 1)               # otto: 8.3
  else:
    mvaCut = np.linspace(8.6, 8.6, 1)               # otto: 8.6
    #mvaCut = np.linspace(0., 0., 1)               # otto: 8.6

  for cut in mvaCut:
    eff_sig_bdt = np.mean([float(jpsi_mc_branches[i].query(' and '.join([selection['jpsi'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_jpsi'] for i in nparts])
    eff_partial_bdt = np.mean([float(partial_mc_branches[i].query(' and '.join([selection['jpsi'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_partial'] for i in nparts])
    eff_kstarplus_bdt = np.mean([float(kstarplus_mc_branches[i].query(' and '.join([selection['jpsi'], selection['Dmass'], '(xgb > @cut)'])).shape[0]) / info[eleType]['n_mc_kstarplus'] for i in nparts])

    frac_ratio = (eff_partial_bdt / eff_sig_bdt) * (br_b2kstarjpsi / br_b2kjpsi)
    frac_ratio_kstarplus = (eff_kstarplus_bdt / eff_partial_bdt) * (br_b2kstarplusjpsi / br_b2kstarjpsi)
    total_ratio = ((eff_partial_bdt + eff_kstarplus_bdt) / eff_sig_bdt) * ((br_b2kstarjpsi + br_b2kstarplusjpsi) / br_b2kjpsi)

    print(frac_ratio, " ", frac_ratio_kstarplus, " ", total_ratio)    

    if eleType == 'pf':
      com = 'python KJpsi_roofit_plb_modified_kde_fixedPartial.py -i {} -o jpsi_fixedPartial_{} --isgn={} --iotherB={} --iKstarJpsi_BKG={} --iKstarPlusJpsi_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarjpsi,bkg_kstarplusjpsi,bkg_otherb,bkg_comb" --fit_primitive --partial_ratio={} --partial_ratio_kstarplus={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iotherB'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['iKstarPlusJpsi_BKG'], info[eleType]['ibkg'], frac_ratio, frac_ratio_kstarplus, cut, log)
    else:
      com = 'python KJpsi_roofit_plb_modified_kde_fixedPartial.py -i {} -o jpsi_fixedPartial_{} --isgn={} --iotherB={} --iKstarJpsi_BKG={} --iKstarPlusJpsi_BKG={} --ibkg={} --sel_primitive="sgn,bkg_kstarjpsi,bkg_kstarplusjpsi,bkg_otherb,bkg_comb" --fit_primitive --partial_ratio={} --partial_ratio_kstarplus={} --mvacut={} --log={}'.format(info[eleType]['inputfile'], eleType, info[eleType]['isgn'], info[eleType]['iotherB'], info[eleType]['iKstarJpsi_BKG'], info[eleType]['iKstarPlusJpsi_BKG'], info[eleType]['ibkg'], frac_ratio, frac_ratio_kstarplus, cut, log)

    os.system(com)



