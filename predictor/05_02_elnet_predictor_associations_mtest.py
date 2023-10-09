#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena
# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH
# source ~/anaconda3/etc/profile.d/conda.sh

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import colorsys
import datatable as dt
import math
from scipy import stats
import statsmodels.stats.multitest as smm
import scipy.stats


### Import data
#########################################################

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/associations/"
lbc_target = read.table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", header = T)

# Covars/lifestyle
assocs_21 = pd.read_table(output_dir + "trait_assocs_episcore_lbc1921.tsv")
assocs_36 = pd.read_table(output_dir + "trait_assocs_episcore_lbc1936.tsv")
assocs_21_ma = pd.read_table(output_dir + "trait_assocs_measuredalc_lbc1921.tsv")
assocs_36_ma = pd.read_table(output_dir + "trait_assocs_measuredalc_lbc1936.tsv")

# Disease
assocs_21_disease = pd.read_table(output_dir + "disease_assocs_episcore_lbc1921.tsv")
assocs_36_disease = pd.read_table(output_dir + "disease_assocs_episcore_lbc1936.tsv")
assocs_21_ma_disease = pd.read_table(output_dir + "disease_assocs_measuredalc_lbc1921.tsv")
assocs_36_ma_disease = pd.read_table(output_dir + "disease_assocs_measuredalc_lbc1936.tsv")

# Biomarkers
assocs_21_bm = pd.read_table(output_dir + "biomarker_assocs_episcore_lbc1921.tsv")
assocs_36_bm = pd.read_table(output_dir + "biomarker_assocs_episcore_lbc1936.tsv")
assocs_21_ma_bm = pd.read_table(output_dir + "biomarker_assocs_measuredalc_lbc1921.tsv")
assocs_36_ma_bm = pd.read_table(output_dir + "biomarker_assocs_measuredalc_lbc1936.tsv")

# MRI
assocs_36_mri = pd.read_table(output_dir + "brainmri_assocs_lbc1936.tsv")

# Survival
assocs_21_surv = pd.read_table(output_dir + "survival_assocs_episcore_lbc1921.tsv")
assocs_36_surv = pd.read_table(output_dir + "survival_assocs_episcore_lbc1936.tsv")
assocs_21_ma_surv = pd.read_table(output_dir + "survival_assocs_measuredalc_lbc1921.tsv")
assocs_36_ma_surv = pd.read_table(output_dir + "survival_assocs_measuredalc_lbc1936.tsv")


### Fuse data
#########################################################

# Prep MRI data
mri_epi = assocs_36_mri.loc[assocs_36_mri["Predictor"] == "acpred_en_ud_sc_t", ["Outcome", "Beta", "SE", "P"]].set_index("Outcome")
mri_ma = assocs_36_mri.loc[assocs_36_mri["Predictor"] == "alcunitsupw_sc_t", ["Outcome", "Beta", "SE", "P"]].set_index("Outcome")
mri_epi.columns = ["Beta", "SE", "p"]
mri_ma.columns = ["Beta", "SE", "p"]

# Prep survival data
assocs_21_surv_ii = assocs_21_surv.loc[assocs_21_surv["term"] == "scale(acpred_en_ud)", ["term", "estimate", "std.error", "p.value"]].set_index("term")
assocs_21_surv_ii.columns = ["Beta", "SE", "p"]
assocs_21_surv_ii.index = ["survival"]
assocs_21_ma_surv_ii = assocs_21_ma_surv.loc[assocs_21_ma_surv["term"] == "scale(alcunitsupw)", ["term", "estimate", "std.error", "p.value"]].set_index("term")
assocs_21_ma_surv_ii.columns = ["Beta", "SE", "p"]
assocs_21_ma_surv_ii.index = ["survival"]

assocs_36_surv_ii = assocs_36_surv.loc[assocs_36_surv["term"] == "scale(acpred_en_ud)", ["term", "estimate", "std.error", "p.value"]].set_index("term")
assocs_36_surv_ii.columns = ["Beta", "SE", "p"]
assocs_36_surv_ii.index = ["survival"]
assocs_36_ma_surv_ii = assocs_36_ma_surv.loc[assocs_36_ma_surv["term"] == "scale(alcunitsupw)", ["term", "estimate", "std.error", "p.value"]].set_index("term")
assocs_36_ma_surv_ii.columns = ["Beta", "SE", "p"]
assocs_36_ma_surv_ii.index = ["survival"]

# Fuse it all
assocs_21_all = pd.concat([assocs_21, assocs_21_disease, assocs_21_bm, assocs_21_surv_ii], axis = 0)
assocs_36_all = pd.concat([assocs_36, assocs_36_disease, assocs_36_bm, assocs_36_surv_ii, mri_epi], axis = 0)
assocs_21_all = assocs_21_all.drop(["demhist", "parkinhist"])
assocs_36_all = assocs_36_all.drop(["dement_w1", "parkin_w1"])

assocs_21_all_ma = pd.concat([assocs_21_ma, assocs_21_ma_disease, assocs_21_ma_bm, assocs_21_ma_surv_ii], axis = 0)
assocs_36_all_ma = pd.concat([assocs_36_ma, assocs_36_ma_disease, assocs_36_ma_bm, assocs_36_ma_surv_ii, mri_ma], axis = 0)
assocs_21_all_ma = assocs_21_all_ma.drop(["demhist", "parkinhist"])
assocs_36_all_ma = assocs_36_all_ma.drop(["dement_w1", "parkin_w1"])


### Multiple testing for LBC1921 and LCB1936 separately
#########################################################

assocs_21_all["p_FDR"] = smm.multipletests(assocs_21_all["p"], method='fdr_bh')[1]
assocs_36_all["p_FDR"] = smm.multipletests(assocs_36_all["p"], method='fdr_bh')[1]
assocs_21_all_ma["p_FDR"] = smm.multipletests(assocs_21_all_ma["p"], method='fdr_bh')[1]
assocs_36_all_ma["p_FDR"] = smm.multipletests(assocs_36_all_ma["p"], method='fdr_bh')[1]
assocs_21_all.columns = ["%s_EpiScore" % i for i in assocs_21_all.columns]
assocs_21_all_ma.columns = ["%s_Measured" % i for i in assocs_21_all_ma.columns]
assocs_36_all.columns = ["%s_EpiScore" % i for i in assocs_36_all.columns]
assocs_36_all_ma.columns = ["%s_Measured" % i for i in assocs_36_all_ma.columns]

assocs_21_final = pd.concat([assocs_21_all, assocs_21_all_ma], axis = 1)
assocs_36_final = pd.concat([assocs_36_all, assocs_36_all_ma], axis = 1)


### Export
#########################################################

assocs_21_final.to_csv(output_dir + "all_assocs_lbc1921.tsv", sep = "\t", index_label = "trait")
assocs_36_final.to_csv(output_dir + "all_assocs_lbc1936.tsv", sep = "\t", index_label = "trait")




