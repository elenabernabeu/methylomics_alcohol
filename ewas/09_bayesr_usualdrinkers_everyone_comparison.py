#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena
# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate py37

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
from statsmodels.graphics.tsaplots import plot_pacf
from statsmodels.graphics.tsaplots import plot_acf
import pymc3 as pm3
import arviz as az
sys.path.append('/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/scripts')      
from plot_funcs import *
from scipy.stats import gaussian_kde
import scipy.stats
from sklearn import metrics

# Color palette
palette = "custom2"

def flatten(xss):
    return [x for xs in xss for x in xs]


# Import data
###########################################################

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/explore_final/"

every = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_16717_logalcohol_residualized_everyone_meanbeta_pip.tsv", index_col = 1) 
usual = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_10506_logalcohol_residualized_usualdrinkers_onelesssample_meanbeta_pip.tsv", index_col = 1) 
every = every[every["PIP"] > 0.8]
usual = usual[usual["PIP"] > 0.8]
every["Model"] = "Full cohort"
usual["Model"] = "Usual drinkers"

# Union of both
###########################################################

union = set(every.index.values.tolist() + usual.index.values.tolist())
intersection = set(every.index.values.tolist()) & set(usual.index.values.tolist())
union_df = pd.concat([every.loc[[i for i in union if i in every.index],], usual.loc[[i for i in union if i in usual.index],]])
union_df = union_df.sort_values(by = "PIP", ascending = False)
union_df.loc[intersection, "Model"] = "Both"
union_df = union_df[~union_df.index.duplicated(keep='first')]
union_df.to_csv(output_dir + "union_bothmodels.tsv", sep = "\t")
union_genes = [i.split(";")[0] for i in set(union_df["UCSC_RefGene_Name"].dropna())]
pd.Series(union_genes).to_csv(output_dir + "union_genes_0.8.tsv", sep = "\t", index = False, header = False)

# Compare betas
###########################################################

r = stats.pearsonr(every["Mean_Beta"], usual["Mean_Beta"]) # (0.5479692895898396, 0.0) for all


# Compare PIPs
###########################################################

r = stats.pearsonr(every["PIP"], usual["PIP"]) # (0.4390568886300766, 0.0)


# Plot
###########################################################


fig, axes = plt.subplots(1, 2, figsize = (9, 3))

# Beta
r = stats.pearsonr(every["Mean_Beta"], usual["Mean_Beta"]) # (0.5479692895898396, 0.0)
m, b = np.polyfit(every["Mean_Beta"], usual["Mean_Beta"], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([every["Mean_Beta"], usual["Mean_Beta"]])
z = gaussian_kde(xy)(xy)
axes[0].scatter(y = every["Mean_Beta"], x = usual["Mean_Beta"], c = z, s = 15, cmap = cmap, rasterized = True)
axes[0].plot(usual["Mean_Beta"], m*usual["Mean_Beta"]+b, color = "black", lw = 1)
axes[0].set_ylabel("Everyone")
axes[0].set_xlabel("Usual Drinkers")
axes[0].set_title("Beta \n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

# PIP
r = stats.pearsonr(every["PIP"], usual["PIP"]) # (0.4390568886300766, 0.0)
m, b = np.polyfit(every["PIP"], usual["PIP"], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([every["PIP"], usual["PIP"]])
z = gaussian_kde(xy)(xy)
axes[1].scatter(y = every["PIP"], x = usual["PIP"], c = z, s = 15, cmap = cmap, rasterized = True)
axes[1].plot(usual["PIP"], m*usual["PIP"]+b, color = "black", lw = 1)
axes[1].set_ylabel("Everyone")
axes[1].set_xlabel("Usual Drinkers")
axes[1].set_title("PIP \n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.savefig(output_dir + "usual_vs_everyone.pdf", dpi = 300)
plt.close(fig)
