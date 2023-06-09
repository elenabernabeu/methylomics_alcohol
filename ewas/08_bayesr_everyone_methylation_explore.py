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

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/sexcomp/"
output_dir_ii = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/explore/"


sexagnos = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/output/summary/gs20k_alcoholconsumption_16689_meanbeta_pip.tsv", index_col = 1) 
sigma_sexagnos_file = "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/output/sigma/gs20k_alcoholconsumption_16689_processed.csv"
sigma_sexagnos = pd.read_table(sigma_sexagnos_file, sep = ",")
sigma_sexagnos["varexp"] = sigma_sexagnos.iloc[:,1]/(sigma_sexagnos.iloc[:,0]+sigma_sexagnos.iloc[:,1])
varexplained_sexagnos = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/output/summary/gs20k_alcoholconsumption_16689_varianceexplained.tsv", index_col = 0)
prop_varexplained_sexagnos = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/output/summary/gs20k_alcoholconsumption_16689_varianceexplained_periteration.tsv")

# Fix chromosome column
sexagnos["chr"] = [int(i.replace("chr", "")) for i in sexagnos["chr"]]


# Convergence
###########################################################

## Parameter values over iterations
i = 0
i_dic = {0 : "sexagnos", 1 : "females", 2 : "males"}
i_dic_ii = {0 : "", 1 : " - Females", 2 : " - Males"}
i_dic_iii = {0 : "Sex-agnostic", 1 : "Females", 2 : "Males"}
for df in [sigma_sexagnos_file, sigma_females_file, sigma_males_file]:
    # Convergence for each seed (sum of sigmas)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3), sharex = True)
    col = set_colors(4, palette)
    axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    for s in range(4):
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
        seed_df["ite"] = range(1,1001)
        axes.plot(seed_df["ite"], seed_df["sum"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 1)
    
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Sigma Sum")
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.title("Sigma Convergence%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s.pdf" % i_dic[i], dpi = 300)
    plt.close()
    
    # Convergence for each seed (proportion variance explained by G)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3), sharex = True)
    col = set_colors(4, palette)
    axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    for s in range(4):
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
        seed_df["ite"] = range(1,1001)
        axes.plot(seed_df["ite"], seed_df["prop"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 1)
    
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Sigma G Prop")
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.title("Sigma G/(Sigma G + Sigma E) Convergence%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_proportion.pdf" % i_dic[i], dpi = 300)
    plt.close()
    
    # Convergence for each seed (individual sigmas)
    fig, axes = plt.subplots(2, 1, figsize = (7, 5), sharex = True)
    col = set_colors(4, palette)
    axes[0].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    axes[1].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    for s in range(4):
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["ite"] = range(1,1001)
        axes[0].plot(seed_df["ite"], seed_df["sigmaE"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 1)
        axes[1].plot(seed_df["ite"], seed_df["sigmaG[1]"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 1)
    
    sns.despine(offset=10, trim=True);
    axes[1].set_xlabel("Iteration")
    axes[0].set_ylabel("Sigma E")
    axes[1].set_ylabel("Sigma G")
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.suptitle("Sigma Convergence%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_perparameter.pdf" % i_dic[i], dpi = 300)
    plt.close()
    
    # Convergence rolling median plot (sum of sigmas)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3))
    col = set_colors(4, palette)
    axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    for s in range(4):
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
        seed_df["ite"] = range(1,1001)
        for j in range(0, len(seed_df.index)):
            if j == 0:
                subset = seed_df.loc[seed_df.index[0]]
            else:
                subset = seed_df.loc[seed_df.index[0:j]]
            seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median"] = median(subset["sum"])
        axes.plot(seed_df["ite"], seed_df["rolling_median"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 2)
    
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Median Sigma G + Sigma E")
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.title("Rolling Median Sigma G + Sigma E%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_rollingmedian.pdf" % i_dic[i], dpi = 300)
    plt.close()
    
    # Convergence rolling median plot (proportion)
    fig, axes = plt.subplots(1, 1, figsize = (7, 3))
    col = set_colors(4, palette)
    axes.axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    for s in range(4):
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
        seed_df["ite"] = range(1,1001)
        for j in range(0, len(seed_df.index)):
            if j == 0:
                subset = seed_df.loc[seed_df.index[0]]
            else:
                subset = seed_df.loc[seed_df.index[0:j]]
            seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median"] = median(subset["prop"])
        axes.plot(seed_df["ite"], seed_df["rolling_median"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 2)
    
    sns.despine(offset=10, trim=True);
    axes.set_xlabel("Iteration")
    axes.set_ylabel("Median Sigma G Prop")
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.title("Rolling Median Sigma G/(Sigma G + Sigma E)%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_rollingmedian_proportion.pdf" % i_dic[i], dpi = 300)
    plt.close()
    
    # Convergence rolling median plot (per sigma)
    fig, axes = plt.subplots(2, 1, figsize = (7, 5), sharex = True)
    col = set_colors(4, palette)
    axes[0].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    axes[1].axvspan(749, 1005, color = 'grey', alpha = 0.15, linewidth = 0, zorder = 1)
    for s in range(4):
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["ite"] = range(1,1001)
        for j in range(0, len(seed_df.index)):
            if j == 0:
                subset = seed_df.loc[seed_df.index[0]]
            else:
                subset = seed_df.loc[seed_df.index[0:j]]
            seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median_sigmaG"] = median(subset["sigmaG[1]"])
            seed_df.loc[seed_df.index.values.tolist()[j], "rolling_median_sigmaE"] = median(subset["sigmaE"])
        axes[0].plot(seed_df["ite"], seed_df["rolling_median_sigmaE"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 2)
        axes[1].plot(seed_df["ite"], seed_df["rolling_median_sigmaG"], linestyle = "solid", c = col[s], zorder = 2, label = "Chain %s" % seed, linewidth = 2)
    
    sns.despine(offset=10, trim=True);
    axes[1].set_xlabel("Iteration")
    axes[0].set_ylabel("Median Sigma E")
    axes[1].set_ylabel("Median Sigma G")
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.suptitle("Rolling Median Sigma G and Sigma E%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_rollingmedian_perparameter.pdf" % i_dic[i], dpi = 300)
    plt.close()
    i += 1

## Histogram per parameter
i = 0
for df in [sigma_sexagnos_file, sigma_females_file, sigma_males_file]:
    fig, axes = plt.subplots(1, 2, figsize = (6, 3), sharex = True, sharey = True)
    col = set_colors(4, palette)
    seed = 1
    f = sigma_sexagnos_file
    seed_df = pd.read_table(f, sep = ",")
    #seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
    #seed_df["ite"] = range(1,1001)
    sns.histplot(seed_df["sigmaE"], ax = axes[0], color = col[0], label = "Chain %s" % seed, alpha = None, linewidth = 0, edgecolor = "white")
    sns.histplot(seed_df["sigmaG[1]"], ax = axes[1], color = col[0], label = "Chain %s" % seed, alpha = None, linewidth = 0, edgecolor = "white")
    axes[0].set_xlabel("Sigma E")
    axes[1].set_xlabel("Sigma G")
    axes[0].set_title("Sigma E")
    axes[1].set_title("Sigma G")
    axes[1].set_ylabel("Count")
    axes[0].set_ylabel("Count")
    #axes[0].set_xlim([0,1])
    #axes[1].set_xlim([0,1])
    sns.despine(offset=10, trim=True);
    #plt.suptitle("Sigma Autocorrelation")
    plt.tight_layout()
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_autocorrelation.pdf" % i_dic[i], dpi = 300)
    plt.close()
    i += 1
       
## Autocorrelation
i = 0
for df in [sigma_sexagnos_file, sigma_females_file, sigma_males_file]:
    # Convergence for each seed (sum of sigmas)
    fig, axes = plt.subplots(1, 4, figsize = (10, 3), sharex = True, sharey = True)
    col = set_colors(4, palette)
    for s in range(4):
        print(s)
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
        seed_df["ite"] = range(1,1001)
        plot_acf(seed_df["sum"], ax = axes[s], color = col[s], label = "Chain %s" % seed, lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
        axes[s].set_title("Chain %s" % seed)
        axes[s].set_ylabel("Autocorrelation")
        axes[s].set_xlabel("Lag")
        axes[s].set_ylim([-0.1,1.1])

    sns.despine(offset=10, trim=True);
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.suptitle("Sigma Sum Autocorrelation%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_autocorrelation.pdf" % i_dic[i], dpi = 300)
    plt.close()

    # Convergence for each seed (proportion)
    fig, axes = plt.subplots(1, 4, figsize = (10, 3), sharex = True, sharey = True)
    col = set_colors(4, palette)
    for s in range(4):
        print(s)
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
        seed_df["ite"] = range(1,1001)
        plot_acf(seed_df["prop"], ax = axes[s], color = col[s], label = "Chain %s" % seed, lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
        axes[s].set_title("Chain %s" % seed)
        axes[s].set_ylabel("Autocorrelation")
        axes[s].set_xlabel("Lag")
        axes[s].set_ylim([-0.1,1.1])

    sns.despine(offset=10, trim=True);
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.suptitle("Sigma G/(Sigma G + Sigma E) Autocorrelation%s" % i_dic_ii[i])
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_autocorrelation_proportion.pdf" % i_dic[i], dpi = 300)
    plt.close()

    # Convergence for each seed (individual sigmas)
    fig, axes = plt.subplots(2, 4, figsize = (10, 5), sharex = True, sharey = True)
    col = set_colors(4, palette)
    for s in range(4):
        print(s)
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["sum"] = seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]]
        seed_df["ite"] = range(1,1001)
        plot_acf(seed_df["sigmaE"], ax = axes[0,s], color = col[s], label = "Chain %s" % seed, lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
        plot_acf(seed_df["sigmaG[1]"], ax = axes[1,s], color = col[s], label = "Chain %s" % seed, lags = 50, alpha = None, vlines_kwargs={"colors": col[s]})
        #axes[0,s].acorr(seed_df["sigmaE"], color = col[s], label = "Chain %s" % seed, maxlags = 50)
        #axes[1,s].acorr(seed_df["sigmaG[1]"], color = col[s], label = "Chain %s" % seed, maxlags = 50)
        axes[0,s].set_title("Chain %s" % seed)
        axes[0,s].set(xlabel = None)
        #axes[1,s].set_xlabel("Lag")
        axes[1,s].set(title = None)

    axes[1,0].set_ylabel("Sigma G")
    axes[0,0].set_ylabel("Sigma E")
    sns.despine(offset=10, trim=True);
    elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
    fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
    plt.suptitle("Sigma Autocorrelation%s" % i_dic_ii[i])
    fig.supylabel("Autocorrelation")
    fig.supxlabel("Lag")
    plt.tight_layout()
    fig.subplots_adjust(right=0.87) 
    fig.savefig(output_dir_ii + "convergence_alcoholconsumption_%s_autocorrelation_perparameter.pdf" % i_dic[i], dpi = 300)
    plt.close()
    i += 1
    
## Geweke diagnostics
i = 0
fig, axes = plt.subplots(3, 1, figsize = (6, 5))
col = set_colors(4, palette)
for df in [sigma_sexagnos_file, sigma_females_file, sigma_males_file]:
    print(i)
    for s in range(4):
        print(s)
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        geweke_sigmaG = pm3.geweke(seed_df["sigmaG[1]"], intervals = 1, first=0.1, last=0.5)
        geweke_sigmaE = pm3.geweke(seed_df["sigmaE"], intervals = 1, first=0.1, last=0.5)
        print(geweke_sigmaG[0][1])
        print(geweke_sigmaE[0][1])
        axes[i].plot(geweke_sigmaG[0][1], 0, marker = 'o', label = "Chain %s" % seed, color = col[s], zorder = 2)
        axes[i].plot(geweke_sigmaE[0][1], 1, marker = 'o', label = "Chain %s" % seed, color = col[s], zorder = 2)
    axes[i].set_xlim([-5, 5])
    axes[i].set_title(i_dic_iii[i])
    axes[i].set_yticks([0, 1])
    axes[i].set_yticklabels(["Sigma G", "Sigma E"])
    axes[i].set_ylim([-0.1, 1.1])
    axes[i].axvline(-2, color = 'grey', alpha = 0.7, linestyle = "dashed", linewidth = 1.5, zorder = 1)
    axes[i].axvline(2, color = 'grey', alpha = 0.7, linestyle = "dashed", linewidth = 1.5, zorder = 1)
    i += 1
axes[2].set_xlabel("Z")
fig.suptitle("Geweke Z Diagnostic")
plt.tight_layout()
sns.despine(offset=10, trim=True);
elements = [mpatches.Patch(edgecolor = col[i], facecolor = col[i], label = "Chain %s" % (i+1)) for i in range(0, 4)]
fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
fig.subplots_adjust(right=0.87) 
fig.savefig(output_dir_ii + "convergence_alcoholconsumption_geweke_perparameter.pdf", dpi = 300)
plt.close()

## Effective sample size
i = 0
ess_df = pd.DataFrame(index = flatten([list(repeat(i_dic[i], 4)) for i in range(0,3)]), columns = ["Parameter"] + ["Chain %s" % i for i in range(1,5)])
ess_df["Parameter"] = ["Sigma G", "Sigma E", "Sigma Sum", "Sigma Prop"]*3
for df in [sigma_sexagnos_file, sigma_females_file, sigma_males_file]:
    print(i)
    for s in range(0,4):
        print(s)
        seed = s + 1
        f = df.replace("_processed", "_seed%s" % seed)
        seed_df = pd.read_table(f, sep = ",")
        seed_df["sum"] = seed_df["sigmaE"] + seed_df["sigmaG[1]"]
        seed_df["prop"] = seed_df[list(seed_df)[0]]/(seed_df[list(seed_df)[0]] + seed_df[list(seed_df)[1]])
        ess_G = az.ess(np.array(seed_df["sigmaG[1]"]))
        ess_E = az.ess(np.array(seed_df["sigmaE"]))
        ess_sum = az.ess(np.array(seed_df["sum"]))
        ess_prop = az.ess(np.array(seed_df["prop"]))
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma G"), "Chain %s" % seed] = ess_G
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma E"), "Chain %s" % seed] = ess_E
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma Sum"), "Chain %s" % seed] = ess_sum
        ess_df.loc[(ess_df.index == i_dic[i]) & (ess_df["Parameter"] == "Sigma Prop"), "Chain %s" % seed] = ess_prop
    i += 1

ess_df.to_csv(output_dir_ii + "ess_acrosschains.tsv", sep = "\t", index_label = "Model", na_rep = "NA")


# Explore metrics
###########################################################

## Manhattan plots PIP
fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = sexagnos["PIP"], meta = sexagnos, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "PIP", colors = "custom2")
axes.set_title("Alcohol Consumption")
axes.axhline(y=0.95, color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir_ii + "manhattan_alcoholconsumption_sexagnos.pdf", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = females["PIP"], meta = females, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "PIP", colors = "custom2")
axes.set_title("Alcohol Consumption - Females")
axes.axhline(y=0.95, color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir_ii + "manhattan_alcoholconsumption_females.pdf", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = males["PIP"], meta = males, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "PIP", colors = "custom2")
axes.set_title("Alcohol Consumption - Males")
axes.axhline(y=0.95, color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir_ii + "manhattan_alcoholconsumption_males.pdf", dpi=300)
plt.close(fig)

## Beta & PIP distributions
for var in ["Mean_Beta", "PIP"]:
    # Sex-agnos
    fig, axes = plt.subplots(1, 2, figsize = (6,3.3))
    plt.suptitle(var)
    histogram(sexagnos, var, ax = axes[0], col = palette)
    densityplot(sexagnos, var, ax = axes[1], col = palette)
    plt.tight_layout()
    plt.savefig(output_dir_ii + "%s_distribution_plots_sexagnostic.pdf" % var)
    plt.close(fig)
    
    # Females
    fig, axes = plt.subplots(1, 2, figsize = (6,3.3))
    plt.suptitle("%s - Females" % var)
    histogram(females, var, ax = axes[0], col = palette)
    densityplot(females, var, ax = axes[1], col = palette)
    plt.tight_layout()
    plt.savefig(output_dir_ii + "%s_distribution_plots_females.pdf" % var)
    plt.close(fig)
    
    # Males
    fig, axes = plt.subplots(1, 2, figsize = (6,3.3))
    plt.suptitle("%s - Males" % var)
    histogram(males, var, ax = axes[0], col = palette)
    densityplot(males, var, ax = axes[1], col = palette)
    plt.tight_layout()
    plt.savefig(output_dir_ii + "%s_distribution_plots_males.pdf" % var)
    plt.close(fig)


# Sig dif in variance explained?
###########################################################

col = set_colors(5, palette)
sns.set_palette(sns.color_palette(col))

## Violin plots
df = pd.DataFrame({"Sex-agnostic\n%s [%s, %s]" % (round(varexplained_sexagnos.iloc[0,0],2), round(varexplained_sexagnos.iloc[0,1],2), round(varexplained_sexagnos.iloc[0,2],2)) : sigma_sexagnos[list(sigma_sexagnos)[2]], "Females\n%s [%s, %s]" % (round(varexplained_females.iloc[0,0],2), round(varexplained_females.iloc[0,1],2), round(varexplained_females.iloc[0,2],2)) : sigma_females[list(sigma_females)[2]], "Males\n%s [%s, %s]" % (round(varexplained_males.iloc[0,0],2), round(varexplained_males.iloc[0,1],2), round(varexplained_males.iloc[0,2],2)) : sigma_males[list(sigma_males)[2]]})
fig, axes = plt.subplots(1, 1, figsize = (5,4))
sns.violinplot(data = df, ax = axes, linewidth = 1, saturation=1)
sns.despine(offset=10, trim=True);
plt.ylabel("Variance Explained")
plt.title("Variance Explained - Methylome")
plt.tight_layout()
plt.savefig(output_dir_ii + "varianceexplained_distribution_plots.pdf")
plt.close(fig)

## Sig difference between sexes?
sigma_m_F = mean(df.iloc[:,1])
sigma_m_M = mean(df.iloc[:,2])
sigma_se_F = stats.sem(df.iloc[:,1])
sigma_se_M = stats.sem(df.iloc[:,2])
t = (sigma_m_F - sigma_m_M)/np.sqrt(sigma_se_F**2 + sigma_se_M**2) # 28.758885998180837
p = 2*stats.norm.sf(abs(t)) # 7.013110999745043e-182


# Variance explained by small, medium and large effects
###########################################################

## Violin plots
fig, axes = plt.subplots(3, 1, figsize = (5,7))
sns.violinplot(data = prop_varexplained_sexagnos[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]], ax = axes[0], linewidth = 1, saturation=1)
sns.violinplot(data = prop_varexplained_females[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]], ax = axes[1], linewidth = 1, saturation=1)
sns.violinplot(data = prop_varexplained_males[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]], ax = axes[2], linewidth = 1, saturation=1)
axes[0].set_title("Sex agnostic")
axes[1].set_title("Females")
axes[2].set_title("Males")
axes[0].set_ylabel("PVE")
axes[1].set_ylabel("PVE")
axes[2].set_ylabel("PVE")
axes[0].set_ylim([-0.1,1])
axes[1].set_ylim([-0.1,1])
axes[2].set_ylim([-0.1,1])
axes[0].set_xticklabels(["0.0001\n(Mean = %s)" % round(varexplained_sexagnos.iloc[0,-3],2), "0.001\n(Mean = %s)" % round(varexplained_sexagnos.iloc[0,-2],2), "0.01\n(Mean = %s)" % round(varexplained_sexagnos.iloc[0,-1],2)])
axes[1].set_xticklabels(["0.0001\n(Mean = %s)" % round(varexplained_females.iloc[0,-3],2), "0.001\n(Mean = %s)" % round(varexplained_females.iloc[0,-2],2), "0.01\n(Mean = %s)" % round(varexplained_females.iloc[0,-1],2)])
axes[2].set_xticklabels(["0.0001\n(Mean = %s)" % round(varexplained_males.iloc[0,-3],2), "0.001\n(Mean = %s)" % round(varexplained_males.iloc[0,-2],2), "0.01\n(Mean = %s)" % round(varexplained_males.iloc[0,-1],2)])
axes[2].set_xlabel("Variance of mixture")
sns.despine(offset=10, trim=True);
#sns.despine(bottom=True, ax = axes[0]);
#sns.despine(bottom=True, ax = axes[1]);
#axes[0].set_xticks([])
#axes[1].set_xticks([])
plt.suptitle("Proportion of Variance Explained by Methylome")
plt.tight_layout()
plt.savefig(output_dir_ii + "varianceexplained_effectsize_distribution_plots.pdf")
plt.close(fig)

## Violin plots for all together
values = pd.concat([prop_varexplained_sexagnos["Proportion_Small_Effects"], prop_varexplained_sexagnos["Proportion_Medium_Effects"], prop_varexplained_sexagnos["Proportion_Large_Effects"], prop_varexplained_females["Proportion_Small_Effects"], prop_varexplained_females["Proportion_Medium_Effects"], prop_varexplained_females["Proportion_Large_Effects"], prop_varexplained_males["Proportion_Small_Effects"], prop_varexplained_males["Proportion_Medium_Effects"], prop_varexplained_males["Proportion_Large_Effects"]])
prop_varexplained_all = pd.DataFrame(data = {"Iteration" : list(range(1,1001))*9, "Size" : list(["0.0001"]*1000 + ["0.001"]*1000 + ["0.01"]*1000)*3, "Model" : ["Sex-Agnostic"]*3000 + ["Females"]*3000 + ["Males"]*3000, "Value" : values})
#prop_varexplained_all = pd.concat([prop_varexplained_sexagnos[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]], prop_varexplained_females[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]], prop_varexplained_males[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]]], axis = 0)
#prop_varexplained_all["group"] = ["Sex-Agnostic"] * 1000 + ["Females"] * 1000 + ["Males"] * 1000
fig, axes = plt.subplots(1, 1, figsize = (7,4))
#sns.violinplot(x = prop_varexplained_all[["Proportion_Small_Effects", "Proportion_Medium_Effects", "Proportion_Large_Effects"]], ax = axes, hue = prop_varexplained_all["group"], linewidth = 1, saturation=1)
sns.violinplot(data = prop_varexplained_all, x = "Size", y = "Value", ax = axes, hue = "Model", linewidth = 1, saturation=1)
axes.set_title("Proportion of Variance Explained by Methylome")
axes.set_ylabel("PVE")
axes.set_ylim([-0.1,1])
#axes.set_xticklabels(["0.0001\n(Mean = %s)" % round(varexplained_sexagnos.iloc[0,-3],2), "0.001\n(Mean = %s)" % round(varexplained_sexagnos.iloc[0,-2],2), "0.01\n(Mean = %s)" % round(varexplained_sexagnos.iloc[0,-1],2)])
#axes[1].set_xticklabels(["0.0001\n(Mean = %s)" % round(varexplained_females.iloc[0,-3],2), "0.001\n(Mean = %s)" % round(varexplained_females.iloc[0,-2],2), "0.01\n(Mean = %s)" % round(varexplained_females.iloc[0,-1],2)])
#axes[2].set_xticklabels(["0.0001\n(Mean = %s)" % round(varexplained_males.iloc[0,-3],2), "0.001\n(Mean = %s)" % round(varexplained_males.iloc[0,-2],2), "0.01\n(Mean = %s)" % round(varexplained_males.iloc[0,-1],2)])
axes.set_xlabel("Variance of mixture")
sns.despine(offset=10, trim=True);
axes.get_legend().remove()
#sns.despine(bottom=True, ax = axes[0]);
#sns.despine(bottom=True, ax = axes[1]);
#axes[0].set_xticks([])
#axes[1].set_xticks([])
elements = [mpatches.Patch(edgecolor = col[0], facecolor = col[0], label = "Sex-Agnostic"), mpatches.Patch(edgecolor = col[1], facecolor = col[1], label = "Females"), mpatches.Patch(edgecolor = col[2], facecolor = col[2], label = "Males")]
fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.tight_layout()
fig.subplots_adjust(right=0.87) 
plt.savefig(output_dir_ii + "varianceexplained_effectsize_distribution_plots_together.pdf")
plt.close(fig)

## Stacked barplots
stacked_df = pd.DataFrame(index = ["Sex Agnostic", "Females", "Males"], columns = ["Small Effect (0.0001)", "Medium Effect (0.001)", "Large Effect (0.01)"])
stacked_df.loc["Sex Agnostic",] = [varexplained_sexagnos.loc["gs20k_alcoholconsumption_16689", "Proportion_%s_Effects" % i] * varexplained_sexagnos.loc["gs20k_alcoholconsumption_16689", "Epigenetic_Mean_Variance_Explained"] for i in ["Small", "Medium", "Large"]]
stacked_df.loc["Males",] = [varexplained_males.loc["gs20k_alcoholconsumption_6954_males", "Proportion_%s_Effects" % i] * varexplained_males.loc["gs20k_alcoholconsumption_6954_males", "Epigenetic_Mean_Variance_Explained"] for i in ["Small", "Medium", "Large"]]
stacked_df.loc["Females",] = [varexplained_females.loc["gs20k_alcoholconsumption_9735_females", "Proportion_%s_Effects" % i] * varexplained_females.loc["gs20k_alcoholconsumption_9735_females", "Epigenetic_Mean_Variance_Explained"] for i in ["Small", "Medium", "Large"]]
# stacked_df.index = ["Sex-agnostic\n%s [%s, %s]" % (round(varexplained_sexagnos.iloc[0,0],2), round(varexplained_sexagnos.iloc[0,1],2), round(varexplained_sexagnos.iloc[0,2],2)), "Females\n%s [%s, %s]" % (round(varexplained_females.iloc[0,0],2), round(varexplained_females.iloc[0,1],2), round(varexplained_females.iloc[0,2],2)), "Males\n%s [%s, %s]" % (round(varexplained_males.iloc[0,0],2), round(varexplained_males.iloc[0,1],2), round(varexplained_males.iloc[0,2],2))]
fig, axes = plt.subplots(1, 1, figsize = (6,4))
stacked_df.plot(kind='bar', stacked=True, ax = axes)
axes.set_title("Proportion of Variance Explained by Methylome")
axes.set_ylabel("Phenotypic Variance Explained")
axes.set_xlabel("Variance of mixture")
#for c in axes.containers:
#    # Optional: if the segment is small or 0, customize the labels
#    labels = [round(v.get_height(),2) if round(v.get_height(),2) > 0 else '' for v in c]
#    # remove the labels parameter if it's not needed for customized labels
#    axes.bar_label(c, labels=labels, label_type='center')
axes.vlines(x = 0, ymin = round(varexplained_sexagnos.iloc[0,1],2), ymax = round(varexplained_sexagnos.iloc[0,2],2), colors="black")
axes.vlines(x = 1, ymin = round(varexplained_females.iloc[0,1],2), ymax = round(varexplained_females.iloc[0,2],2), colors="black")
axes.vlines(x = 2, ymin = round(varexplained_males.iloc[0,1],2), ymax = round(varexplained_males.iloc[0,2],2), colors="black")
sns.despine(offset=10, trim=True);
axes.get_legend().remove()
elements = [mpatches.Patch(edgecolor = col[0], facecolor = col[0], label = "Small Effect\n(0.0001)"), mpatches.Patch(edgecolor = col[1], facecolor = col[1], label = "Medium Effect\n(0.001)"), mpatches.Patch(edgecolor = col[2], facecolor = col[2], label = "Large Effect\n(0.01)")]
fig.legend(handles=elements, loc = "upper right", ncol = 1, fancybox = False, frameon = False, bbox_to_anchor=(1, 0.93))
plt.tight_layout()
fig.subplots_adjust(right=0.77) 
plt.savefig(output_dir_ii + "varianceexplained_effectsize_distribution_plots_stackedbarplot.pdf")
plt.close(fig)


# Beta comp F v M
###########################################################

resultdf = females[list(females)[0:4]]
m_colname = "Mean_Beta"
se_colname = "SE"
m_F = females[m_colname]
m_M = males[m_colname]
se_F = females[se_colname]
se_M = males[se_colname]

# Simple comp formula (no third term in denominator like in NG paper)
resultdf["beta_F"] = m_F
resultdf["beta_M"] = m_M
resultdf["beta_dif"] = m_F - m_M
resultdf["se_F"] = se_F
resultdf["se_M"] = se_M
resultdf["PIP_F"] = females["PIP"]
resultdf["PIP_M"] = males["PIP"]
resultdf["sexcomp_t"] = (m_M - m_F)/np.sqrt(se_M**2 + se_F**2)
resultdf["sexcomp_p"] = 2*stats.norm.sf(resultdf["sexcomp_t"].abs()) # t to norm when huge degrees of freedom
resultdf.loc[resultdf["sexcomp_p"].isna(),"sexcomp_p"] = 1 # There are some missing values
resultdf["sexcomp_p_fdr"] = smm.multipletests(resultdf["sexcomp_p"], method='fdr_bh')[1]
resultdf = resultdf.sort_values(by = "sexcomp_p", ascending = True)
resultdf.to_csv(output_dir + "gs20k_alcoholconsumption_methylation_sexcomp.tsv", sep = "\t")

# Sex diffs in beta?
thresh = 3.6*10**(-8)
resultdf_sig = resultdf.loc[resultdf["sexcomp_p"] < thresh,] # 1409
resultdf_sig_fdr5 = resultdf.loc[resultdf["sexcomp_p_fdr"] < 0.05,] # 5671
resultdf_sig.to_csv(output_dir + "gs20k_alcoholconsumption_methylation_sexcomp_sig.tsv", sep = "\t")
resultdf_sig_fdr5.to_csv(output_dir + "gs20k_alcoholconsumption_methylation_sexcomp_sig_fdr0.05.tsv", sep = "\t")

# Manhattan plot
fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = -np.log10(resultdf["sexcomp_p"]), meta = resultdf, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "-log10(P-Value)", colors = "custom2")
axes.set_title("Alcohol Consumption Sex Comparison")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir + "manhattan_alcoholconsumption_sexcomp.pdf", dpi=300)
plt.close(fig)

# Manhattan plot beta diff
fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = resultdf["beta_dif"], meta = resultdf, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "Beta F - Beta M", colors = "custom2")
axes.set_title("Alcohol Consumption Sex Comparison")
#axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
#sns.despine(top = True, right = True, left = True)
sns.despine(offset=10, trim=True);
plt.tight_layout()
fig.savefig(output_dir + "manhattan_alcoholconsumption_sexcomp_betadif.pdf", dpi=300)
plt.close(fig)

# Divide into female and male dominant
resultdf_sig_Fdom = resultdf_sig.loc[resultdf_sig["beta_F"].abs() > resultdf_sig["beta_M"].abs(),] # 969
resultdf_sig_Mdom = resultdf_sig.loc[resultdf_sig["beta_F"].abs() < resultdf_sig["beta_M"].abs(),] # 440
resultdf_sig_Fdom_fdr = resultdf_sig_fdr5.loc[resultdf_sig_fdr5["beta_F"].abs() > resultdf_sig_fdr5["beta_M"].abs(),] # 4036
resultdf_sig_Mdom_fdr = resultdf_sig_fdr5.loc[resultdf_sig_fdr5["beta_F"].abs() < resultdf_sig_fdr5["beta_M"].abs(),] # 1635

# Sex differences but limiting to just PIP above certain threshold
threshs = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
thresh_sigdf = pd.DataFrame(index = threshs, columns = ["SigDif_Nominal", "SigDif_FDR"])
for thresh in threshs:
    resultdf_sig_thresh = resultdf_sig.loc[(resultdf_sig["PIP_F"] > thresh) | (resultdf_sig["PIP_M"] > thresh)]
    resultdf_sig_thresh_fdr = resultdf_sig_fdr5.loc[(resultdf_sig_fdr5["PIP_F"] > thresh) | (resultdf_sig_fdr5["PIP_M"] > thresh)]
    thresh_sigdf.loc[thresh, "SigDif_Nominal"] = len(resultdf_sig_thresh.index)
    thresh_sigdf.loc[thresh, "SigDif_FDR"] = len(resultdf_sig_thresh_fdr.index)
    # Also output the table
    resultdf_sig_thresh.to_csv(output_dir + "gs20k_alcoholconsumption_methylation_sexcomp_sig_pip%s.tsv" % thresh, sep = "\t")
    resultdf_sig_thresh_fdr.to_csv(output_dir + "gs20k_alcoholconsumption_methylation_sexcomp_sig_pip%s_fdr.tsv" % thresh, sep = "\t")

thresh_sigdf.to_csv(output_dir + "gs20k_alcoholconsumption_methylation_sexcomp_pip_counts.tsv", sep = "\t", index_label = "PIP>")


# Lists of genes for gene set enrichment analysis (FUMA)
###########################################################

# From each EWAS
pip_thresh = 0.8
sexagnos_genes = list(set(flatten([i.split(";") for i in sexagnos.loc[(sexagnos["PIP"] > pip_thresh) & (sexagnos["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))
females_genes = list(set(flatten([i.split(";") for i in females.loc[(females["PIP"] > pip_thresh) & (females["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))
males_genes = list(set(flatten([i.split(";") for i in males.loc[(males["PIP"] > pip_thresh) & (males["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))

# For sex diff
sexdif_genes = list(set(flatten([i.split(";") for i in resultdf_sig.loc[(resultdf_sig["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))
sexdif_genes_femaledom = list(set(flatten([i.split(";") for i in resultdf_sig_Fdom.loc[(resultdf_sig_Fdom["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))
sexdif_genes_maledom = list(set(flatten([i.split(";") for i in resultdf_sig_Mdom.loc[(resultdf_sig_Mdom["UCSC_RefGene_Name"].notna()), "UCSC_RefGene_Name"]])))

# Export lists
new_file = open(output_dir_iii + "alcoholconsumption_sexagnostic_genes.txt", "w")
new_file.write("\n".join(sexagnos_genes))
new_file.close()

new_file = open(output_dir_iii + "alcoholconsumption_females_genes.txt", "w")
new_file.write("\n".join(females_genes))
new_file.close()

new_file = open(output_dir_iii + "alcoholconsumption_males_genes.txt", "w")
new_file.write("\n".join(males_genes))
new_file.close()

new_file = open(output_dir_iii + "alcoholconsumption_sexdif_genes.txt", "w")
new_file.write("\n".join(sexdif_genes))
new_file.close()

new_file = open(output_dir_iii + "alcoholconsumption_sexdif_femaledom_genes.txt", "w")
new_file.write("\n".join(sexdif_genes_femaledom))
new_file.close()

new_file = open(output_dir_iii + "alcoholconsumption_sexdif_maledom_genes.txt", "w")
new_file.write("\n".join(sexdif_genes_maledom))
new_file.close()


# Masking?
###########################################################

mask_threshs = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95]
pip_df = pd.concat([sexagnos[["chr", "pos", "UCSC_RefGene_Name", "strand", "PIP"]], females["PIP"], males["PIP"]], axis = 1)
pip_df.columns = ["Chrom", "Position", "Gene", "Strand", "SexAg_PIP", "F_PIP", "M_PIP"]
mask_df = pd.DataFrame(index = mask_threshs, columns = ["F", "M", "All", "FnotAll", "MnotAll", "ForMnotAll", "FandMnotAll"])

for i in mask_threshs:
    # Get values
    F = pip_df.loc[pip_df["F_PIP"] > i]
    M = pip_df.loc[pip_df["M_PIP"] > i]
    All = pip_df.loc[pip_df["SexAg_PIP"] > i]
    FnotAll = pip_df.loc[(pip_df["SexAg_PIP"] < i) & (pip_df["F_PIP"] > i),]
    MnotAll = pip_df.loc[(pip_df["SexAg_PIP"] < i) & (pip_df["M_PIP"] > i),]
    ForMnotAll = pip_df.loc[(pip_df["SexAg_PIP"] < i) & ((pip_df["F_PIP"] > i) | (pip_df["M_PIP"] > i)),]
    FandMnotAll = pip_df.loc[(pip_df["SexAg_PIP"] < i) & ((pip_df["F_PIP"] > i) & (pip_df["M_PIP"] > i)),]
    # Fill table
    mask_df.loc[i, "F"] = len(F.index)
    mask_df.loc[i, "M"] = len(M.index)
    mask_df.loc[i, "All"] = len(All.index)
    mask_df.loc[i, "FnotAll"] = len(FnotAll.index)
    mask_df.loc[i, "MnotAll"] = len(MnotAll.index)
    mask_df.loc[i, "ForMnotAll"] = len(ForMnotAll.index)
    mask_df.loc[i, "FandMnotAll"] = len(FandMnotAll.index)
    # Export F or M not All table
    ForMnotAll.to_csv(output_dir_iv + "alcoholconsumption_masked_PIP%s.tsv" % i, sep = "\t")

mask_df.to_csv(output_dir_iv + "alcoholconsumption_masking.tsv", sep = "\t") 