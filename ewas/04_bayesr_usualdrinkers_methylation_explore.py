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

# Color palette
palette = "custom2"

def flatten(xss):
    return [x for xs in xss for x in xs]


# Import data
###########################################################

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/explore_final/"
output_dir_ii = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/gene_enrichment/"

usual = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_10506_logalcohol_residualized_usualdrinkers_onelesssample_meanbeta_pip.tsv", index_col = 1) 
sigma_usual_file = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/sigma/allindividuals_10506_logalcohol_residualized_usualdrinkers_onelesssample_processed.csv"
sigma_usual = pd.read_table(sigma_usual_file, sep = ",", header = None)
sigma_usual["varexp"] = sigma_usual.iloc[:,1]/(sigma_usual.iloc[:,0]+sigma_usual.iloc[:,1])
varexplained_usual = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_10506_logalcohol_residualized_usualdrinkers_onelesssample_varianceexplained.tsv", index_col = 0)
prop_varexplained_usual = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_10506_logalcohol_residualized_usualdrinkers_onelesssample_varianceexplained_periteration.tsv")

# Fix chromosome column
usual["chr"] = [int(i.replace("chr", "")) for i in usual["chr"]]


# Convergence
###########################################################

## Parameter values over iterations
i = 0
i_dic = {0: "usualdrinkers"}
seed = 1
s = 0

for df in [sigma_usual_file]:
    # Convergence for each seed (sum of sigmas)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3), sharex = True)
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
    seed_df["ite"] = range(1,1001)
    axes.plot(seed_df["ite"], seed_df["sum"], linestyle = "solid", c = col[s], zorder = 2, linewidth = 1)
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Sigma Sum")
    plt.title("Sigma Convergence")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s.pdf" % i_dic[i], dpi = 300)
    plt.close()
    # Convergence for each seed (proportion variance explained by G)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3), sharex = True)
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["prop"] = seed_df[list(seed_df)[1]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
    seed_df["ite"] = range(1,1001)
    axes.plot(seed_df["ite"], seed_df["prop"], linestyle = "solid", c = col[s], zorder = 2, linewidth = 1)
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Sigma G Prop")
    plt.title("Sigma G/(Sigma G + Sigma E) Convergence")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_proportion.pdf" % i_dic[i], dpi = 300)
    plt.close()
    # Convergence for each seed (individual sigmas)
    fig, axes = plt.subplots(2, 1, figsize = (7, 5), sharex = True)
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["ite"] = range(1,1001)
    axes[0].plot(seed_df["ite"], seed_df[0], linestyle = "solid", c = col[s], zorder = 2, linewidth = 1)
    axes[1].plot(seed_df["ite"], seed_df[1], linestyle = "solid", c = col[s], zorder = 2, linewidth = 1)
    sns.despine(offset=10, trim=True);
    axes[1].set_xlabel("Iteration")
    axes[0].set_ylabel("Sigma E")
    axes[1].set_ylabel("Sigma G")
    plt.suptitle("Sigma Convergence")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_perparameter.pdf" % i_dic[i], dpi = 300)
    plt.close()
    # Convergence rolling median plot (sum of sigmas)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3))
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
    seed_df["ite"] = range(1,1001)
    for j in range(0, len(seed_df.index)):
        if j == 0:
            subset = seed_df.loc[seed_df.index[0]]
        else:
            subset = seed_df.loc[seed_df.index[0:j]]
        seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median"] = median(subset["sum"])
    axes.plot(seed_df["ite"], seed_df["rolling_median"], linestyle = "solid", c = col[s], zorder = 2, linewidth = 2)
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Median Sigma G + Sigma E")
    plt.title("Rolling Median Sigma G + Sigma E")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_rollingmedian.pdf" % i_dic[i], dpi = 300)
    plt.close()
    # Convergence rolling median plot (proportion)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3))
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["prop"] = seed_df[list(seed_df)[1]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
    seed_df["ite"] = range(1,1001)
    for j in range(0, len(seed_df.index)):
        if j == 0:
            subset = seed_df.loc[seed_df.index[0]]
        else:
            subset = seed_df.loc[seed_df.index[0:j]]
        seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median"] = median(subset["prop"])
    axes.plot(seed_df["ite"], seed_df["rolling_median"], linestyle = "solid", c = col[s], zorder = 2, linewidth = 2)
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Median Sigma G Prop")
    plt.title("Rolling Median Sigma G/(Sigma G + Sigma E)")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_rollingmedian_proportion.pdf" % i_dic[i], dpi = 300)
    plt.close()
    # Convergence rolling median plot (per sigma)
    fig, axes = plt.subplots(2, 1, figsize = (7, 5), sharex = True)
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["ite"] = range(1,1001)
    for j in range(0, len(seed_df.index)):
        if j == 0:
            subset = seed_df.loc[seed_df.index[0]]
        else:
            subset = seed_df.loc[seed_df.index[0:j]]
        seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median_sigmaG"] = median(subset[1])
        seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median_sigmaE"] = median(subset[0])
    axes[0].plot(seed_df["ite"], seed_df["rolling_median_sigmaE"], linestyle = "solid", c = col[s], zorder = 2, linewidth = 2)
    axes[1].plot(seed_df["ite"], seed_df["rolling_median_sigmaG"], linestyle = "solid", c = col[s], zorder = 2, linewidth = 2)
    sns.despine(offset=10, trim=True);
    axes[1].set_xlabel("Iteration")
    axes[0].set_ylabel("Median Sigma E")
    axes[1].set_ylabel("Median Sigma G")
    plt.suptitle("Rolling Median Sigma G and Sigma E")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_rollingmedian_perparameter.pdf" % i_dic[i], dpi = 300)
    plt.close()
    i += 1

       
## Autocorrelation
i = 0
for df in [sigma_usual_file]:
    # Convergence for each seed (sum of sigmas)
    fig, axes = plt.subplots(1, 1, figsize = (3, 3), sharex = True, sharey = True)
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(df, sep = ",", header = None)
    seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
    seed_df["ite"] = range(1,1001)
    plot_acf(seed_df["sum"], ax = axes, color = col[s], lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
    axes.set_title("Sigma Sum Autocorrelation")
    axes.set_ylabel("Autocorrelation")
    axes.set_xlabel("Lag")
    axes.set_ylim([-0.1,1.1])
    sns.despine(offset=10, trim=True);
    plt.tight_layout()
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_autocorrelation.pdf" % i_dic[i], dpi = 300)
    plt.close()
    # Convergence for each seed (individual sigmas)
    fig, axes = plt.subplots(2, 1, figsize = (5, 5), sharex = True, sharey = True)
    col = set_colors(4, palette)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
    seed_df["ite"] = range(1,1001)
    plot_acf(seed_df[1], ax = axes[1], color = col[s], lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
    plot_acf(seed_df[0], ax = axes[0], color = col[s], lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
    axes[0].set(xlabel = None)
    #axes[1,s].set_xlabel("Lag")
    axes[1].set(title = None)
    axes[1].set_ylabel("Sigma G")
    axes[0].set_ylabel("Sigma E")
    sns.despine(offset=10, trim=True);
    plt.suptitle("Sigma Autocorrelation")
    fig.supylabel("Autocorrelation")
    fig.supxlabel("Lag")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir + "convergence_alcoholconsumption_%s_autocorrelation_perparameter.pdf" % i_dic[i], dpi = 300)
    plt.close()
    i += 1
    

## Geweke diagnostics
i = 0
fig, axes = plt.subplots(1, 1, figsize = (5, 2))
col = set_colors(4, palette)
for df in [sigma_usual_file]:
    print(i)
    f = df
    seed_df = pd.read_table(f, sep = ",", header = None)
    geweke_sigmaG = pm3.geweke(seed_df[1], intervals = 1, first=0.1, last=0.5)
    geweke_sigmaE = pm3.geweke(seed_df[0], intervals = 1, first=0.1, last=0.5)
    print(geweke_sigmaG[0][1])
    print(geweke_sigmaE[0][1])
    axes.plot(geweke_sigmaG[0][1], 0, marker = 'o', label = "Chain %s" % seed, color = col[s], zorder = 2)
    axes.plot(geweke_sigmaE[0][1], 1, marker = 'o', label = "Chain %s" % seed, color = col[s], zorder = 2)
    axes.set_xlim([-5, 5])
    axes.set_title("Geweke Z Diagnostic")
    axes.set_yticks([0, 1])
    axes.set_yticklabels(["Sigma G", "Sigma E"])
    axes.set_ylim([-0.1, 1.1])
    axes.axvline(-2, color = 'grey', alpha = 0.7, linestyle = "dashed", linewidth = 1.5, zorder = 1)
    axes.axvline(2, color = 'grey', alpha = 0.7, linestyle = "dashed", linewidth = 1.5, zorder = 1)
    i += 1

axes.set_xlabel("Z")
#fig.suptitle("Geweke Z Diagnostic")
plt.tight_layout()
sns.despine(offset=10, trim=True);
fig.savefig(output_dir + "convergence_alcoholconsumption_usualdrinkers_geweke_perparameter.pdf", dpi = 300)
plt.close()


## Effective sample size
i = 0
ess_df = pd.DataFrame(index = flatten([list(repeat(i_dic[i], 4)) for i in range(0,1)]), columns = ["Parameter"] + ["Chain %s" % i for i in range(1,2)])
ess_df["Parameter"] = ["Sigma G", "Sigma E", "Sigma Sum", "Sigma Prop"]
for df in [sigma_usual_file]:
    print(i)
    for s in range(0,1):
        print(s)
        seed = s + 1
        f = df
        seed_df = pd.read_table(f, sep = ",", header = None)
        seed_df["sum"] = seed_df[1] + seed_df[0]
        seed_df["prop"] = seed_df[list(seed_df)[1]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
        ess_G = az.ess(np.array(seed_df[0]))
        ess_E = az.ess(np.array(seed_df[1]))
        ess_sum = az.ess(np.array(seed_df["sum"]))
        ess_prop = az.ess(np.array(seed_df["prop"]))
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma G"), "Chain %s" % seed] = ess_G
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma E"), "Chain %s" % seed] = ess_E
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma Sum"), "Chain %s" % seed] = ess_sum
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma Prop"), "Chain %s" % seed] = ess_prop
    i += 1

ess_df.to_csv(output_dir + "ess_acrosschains_usualdrinkers.tsv", sep = "\t", index_label = "Model", na_rep = "NA")


# Explore metrics
###########################################################

## Manhattan plots PIP
fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = usual["PIP"], meta = usual, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "PIP", colors = "custom2")
axes.set_title("Alcohol Consumption")
axes.axhline(y=0.95, color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir + "manhattan_alcoholconsumption_usual.pdf", dpi=300)
plt.close(fig)


# Lists of genes for gene set enrichment analysis (FUMA)
###########################################################

# From each EWAS
pip_thresh = 0.6
usual_genes = list(set(flatten([i.split(";") for i in usual.loc[(usual["PIP"] > pip_thresh) & (usual["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))

# Export lists
new_file = open(output_dir + "alcoholconsumption_usual_genes_0.6.txt", "w")
new_file.write("\n".join(usual_genes))
new_file.close()
