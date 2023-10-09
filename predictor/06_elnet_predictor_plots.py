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

# Font config
font = "Galvji"
path = '/Cluster_Filespace/Marioni_Group/Elena/other/fonts/Galvji/Galvji-01.ttf'
fe = matplotlib.font_manager.FontEntry(
    fname= path,
    name=font)
matplotlib.font_manager.fontManager.ttflist.insert(0, fe)
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = font

def set_colors(N, name = "default"):
    colors = []
    if name == "default":
        colors = sns.color_palette("Set2", 6)
    elif name == "custom1":
        colors = ["#DD1155", "#7D82B8", "#FCBA04", "#63C7B2", "#45364B"] #Hotpink, lilacky-blue, teal, yellow, dark purple
    elif name == "custom2":
        colors = ["#47C3F4", "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B", "#07B1D8", "#C5E3EA", "#326B64", "#FFD15F", "#EA4C15"] #Lightblue, orange, teal, yellow, red, dark blue, midtone blue, very light blue, 
    colors = (colors*int(math.ceil(N/len(colors))))[:N]
    return colors

def lighten_color(color, amount=0.5):  
    # --------------------- SOURCE: @IanHincks ---------------------
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def rmse_(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def mae_(predictions, targets):
    return(np.median(abs(predictions - targets)))


#### USUAL DRINKERS VS EVERYONE PREDICTOR -- FILTERED VS NON-FILTERED LBC
#####################################################################################
#####################################################################################

lbc_target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", index_col = 0)

lbc_target_36 = lbc_target.loc[(lbc_target["WAVE"] == 1) & (lbc_target["cohort"] == "LBC36"),] # 895
lbc_target_21 = lbc_target.loc[(lbc_target["WAVE"] == 1) & (lbc_target["cohort"] == "LBC21"),] # 436

## Elnet
lbc21_en_usual = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth.tsv", index_col = 0)
lbc36_en_usual = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth.tsv", index_col = 0)
lbc21_en_every = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth.tsv", index_col = 0)
lbc36_en_every = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth.tsv", index_col = 0)

## Elnet filtered CpGs
lbc21_en_usual_filt = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)
lbc36_en_usual_filt = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)
lbc21_en_every_filt = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)
lbc36_en_every_filt = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)

## Add to target
lbc_target_21["acpred_en_ud"] = lbc21_en_usual["ac_pred"]
lbc_target_21["acpred_en_ev"] = lbc21_en_every["ac_pred"]
lbc_target_21["acpred_en_ud_filt"] = lbc21_en_usual_filt["ac_pred"]
lbc_target_21["acpred_en_ev_filt"] = lbc21_en_every_filt["ac_pred"]

lbc_target_36["acpred_en_ud"] = lbc36_en_usual["ac_pred"]
lbc_target_36["acpred_en_ev"] = lbc36_en_every["ac_pred"]
lbc_target_36["acpred_en_ud_filt"] = lbc36_en_usual_filt["ac_pred"]
lbc_target_36["acpred_en_ev_filt"] = lbc36_en_every_filt["ac_pred"]

## All together
df = pd.concat([lbc_target_21, lbc_target_36], axis = 0)

## Pheno
df["alcunitsupw_log"] = np.log(df["alcunitsupw"]+1)
lbc_target_21["alcunitsupw_log"] = np.log(lbc_target_21["alcunitsupw"]+1)
lbc_target_36["alcunitsupw_log"] = np.log(lbc_target_36["alcunitsupw"]+1)


## Scatterplot: Usual drinkers + everyone drinker (elnet + filtering, all cohorts) in LBC
##############################################################################

fig, axes = plt.subplots(2, 2, figsize = (8, 6), sharex = True, sharey = True)

# Elnet usual
true_col = "alcunitsupw_log"
pred_col = "acpred_en_ud"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[0,0].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,0].set_title("Usual Drinkers \n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

# Elnet everyone
true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[0,1].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,1].set_title("Everyone\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

# Elnet usual - Filtered CpGs
true_col = "alcunitsupw_log"
pred_col = "acpred_en_ud_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[1,0].set_ylabel("AC EpiScore")
axes[1,0].set_xlabel("Alcohol Consumption")
axes[1,0].set_title("Usual Drinkers | Filtered CpGs\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

# Elnet everyone - Filtered CpGs
true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[1,1].set_ylabel("AC EpiScore")
axes[1,1].set_xlabel("Alcohol Consumption")
axes[1,1].set_title("Everyone | Filtered CpGs\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc.pdf", dpi = 300)
plt.close(fig)


## Scatterplot: Usual drinkers + everyone drinker (elnet + bayesR, LBC21)
##############################################################################

df = lbc_target_21
fig, axes = plt.subplots(2, 2, figsize = (8, 6), sharex = True, sharey = True)

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ud"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"]) # (0.4230203658674733, 2.348358481612314e-20)
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[0,0].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,0].set_title("Usual Drinkers\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[0,1].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,1].set_title("Everyone\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ud_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[1,0].set_ylabel("AC EpiScore")
axes[1,0].set_xlabel("Alcohol Consumption")
axes[1,0].set_title("Usual Drinkers | Filtered CpGs\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[1,1].set_ylabel("AC EpiScore")
axes[1,1].set_xlabel("Alcohol Consumption")
axes[1,1].set_title("Everyone | Filtered CpGs\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc_1921.pdf", dpi = 300)
plt.close(fig)


## Scatterplot: Usual drinkers + everyone drinker (elnet + bayesR, LBC36)
##############################################################################

df = lbc_target_36
fig, axes = plt.subplots(2, 2, figsize = (8, 6), sharex = True, sharey = True)

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ud"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"]) # (0.4230203658674733, 2.348358481612314e-20)
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[0,0].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,0].set_title("Usual Drinkers\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[0,1].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,1].set_title("Everyone\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ud_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[1,0].set_ylabel("AC EpiScore")
axes[1,0].set_xlabel("Alcohol Consumption")
axes[1,0].set_title("Usual Drinkers | Filtered CpGs\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[1,1].set_ylabel("AC EpiScore")
axes[1,1].set_xlabel("Alcohol Consumption")
axes[1,1].set_title("Everyone | Filtered CpGs\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc_1936.pdf", dpi = 300)
plt.close(fig)


## Scatterplot: Final model
##############################################################################

fig, axes = plt.subplots(1, 2, figsize = (9, 3), sharex = True, sharey = True)

df = lbc_target_21
true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"]) # (0.4230203658674733, 2.348358481612314e-20)
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[0].set_ylabel("AC EpiScore")
axes[0].set_xlabel("Alcohol Consumption")
axes[0].set_title("LBC 1921 \n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

df = lbc_target_36
true_col = "alcunitsupw_log"
pred_col = "acpred_en_ev_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
axes[1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[0,1].set_ylabel("AC EpiScore")
axes[1].set_xlabel("Alcohol Consumption")
axes[1].set_title("LBC 1936\n(r = %s)" % round(r[0],2))
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc_finalmodel.pdf", dpi = 300)
plt.close(fig)


#### GS WAVE 4 - EVERYONE vs USUAL, FILTERED vs NON-FILTERED
#####################################################################################
#####################################################################################

gs_target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/alcohol_16717.tsv", index_col = 1)
gs_target = gs_target.loc[gs_target["Set"] == "W4",] # 8033

gsw4_en = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_usualdrinkerpredictor_noadjustments_scaledmeth.tsv", index_col = 0)
gsw4_en_all = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_noadjustments_scaledmeth.tsv", index_col = 0)
gsw4_en_filt = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_usualdrinkerspredictor_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)
gsw4_en_all_filt = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)
gs_target["ac_pred"] = gsw4_en["ac_pred"]
gs_target["ac_pred_filt"] = gsw4_en_filt["ac_pred"]
gs_target["ac_pred_all"] = gsw4_en_all["ac_pred"]
gs_target["ac_pred_all_filt"] = gsw4_en_all_filt["ac_pred"]

# Plots
df = gs_target
fig, axes = plt.subplots(2, 2, figsize = (8, 6), sharex = True, sharey = True)

true_col = "alcohol_units_log"
pred_col = "ac_pred"
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units_log"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap, rasterized = True)
axes[0,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[0,0].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,0].set_title("Usual Drinkers \n(r = %s)" % round(r[0],3))
sns.despine(offset=10, trim=True);

true_col = "alcohol_units_log"
pred_col = "ac_pred_all"
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units_log"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[0,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap, rasterized = True)
axes[0,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[0,1].set_ylabel("AC EpiScore")
#axes[0].set_xlabel("Alcohol Consumption")
axes[0,1].set_title("Everyone\n(r = %s)" % round(r[0],3))
sns.despine(offset=10, trim=True);

true_col = "alcohol_units_log"
pred_col = "ac_pred_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units_log"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap, rasterized = True)
axes[1,0].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
axes[1,0].set_ylabel("AC EpiScore")
axes[1,0].set_xlabel("Alcohol Consumption")
axes[1,0].set_title("Usual Drinkers | Filtered CpGs\n(r = %s)" % round(r[0],3))
sns.despine(offset=10, trim=True);

true_col = "alcohol_units_log"
pred_col = "ac_pred_all_filt"
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units_log"])
m, b = np.polyfit(df[true_col], df[pred_col], 1)
cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
xy = np.vstack([df[true_col], df[pred_col]])
z = gaussian_kde(xy)(xy)
axes[1,1].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap, rasterized = True)
axes[1,1].plot(df[true_col], m*df[true_col]+b, color = "black", lw = 1)
#axes[1,1].set_ylabel("AC EpiScore")
axes[1,1].set_xlabel("Alcohol Consumption")
axes[1,1].set_title("Everyone | Filtered CpGs\n(r = %s)" % round(r[0],3))
sns.despine(offset=10, trim=True);

plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4.pdf", dpi = 300)
plt.close(fig)


#### CATEGORIES DIFFERENCES - ELNET
#####################################################################################
#####################################################################################

# Separate
gs_target_cats = gs_target.dropna(subset = ["alcohol_usual"]) # 7642/8033
gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2, "Category"] = "Usual drinking units" # 4888/7642
gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1, "Category"] = "Drank more than usual" # 1920/7642
gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3, "Category"] = "Drank less than usual" # 834/7642

## Scatterplot for each category (usual, frequent, less), tested in w4
col = set_colors(3, palette)
sns.set_palette(sns.color_palette(col))
true_col = "alcohol_units_log"
pred_col = "ac_pred_all_filt"
df = gs_target_cats
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
fig, axes = plt.subplots(1, 1, figsize = (8, 4))
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Usual (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[0], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "More (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[1], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Less (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[2], 1.3), levels = 4)
sns.despine(offset=10, trim=True);
plt.legend(frameon = False)
plt.ylabel("AC EpiScore")
plt.xlabel("Alcohol Consumption")
plt.title("EpiScore across categories")
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories.pdf", dpi = 300)
plt.close(fig)

## Scatterplot for each category (usual, frequent, less), tested in w4, SEPARATE PLOTS
fig, axes = plt.subplots(1, 3, figsize = (9, 3), sharey = True)
# Usual
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = lighten_color(col[0], 1.3))
axes[0].set_ylabel("AC EpiScore")
axes[0].set_xlabel("Alcohol Consumption")
axes[0].set_title("Usual drinking\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
# More
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = lighten_color(col[1], 1.3))
axes[1].set_xlabel("Alcohol Consumption")
axes[1].set_title("Drank more than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[1].set_ylabel(None)
# Less
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = lighten_color(col[2], 1.3))
axes[2].set_xlabel("Alcohol Consumption")
axes[2].set_title("Drank less than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[2].set_ylabel(None)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_separate.pdf", dpi = 300)
plt.close(fig)

#### CATEGORIES DIFFERENCES - FOR OTHER MODELS
#####################################################################################
#####################################################################################

## USUAL DRINKERS, NOT FILTERED
col = set_colors(3, palette)
sns.set_palette(sns.color_palette(col))
true_col = "alcohol_units_log"
pred_col = "ac_pred"
df = gs_target_cats
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
fig, axes = plt.subplots(1, 1, figsize = (8, 4))
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Usual (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[0], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "More (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[1], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Less (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[2], 1.3), levels = 4)
sns.despine(offset=10, trim=True);
plt.legend(frameon = False)
plt.ylabel("AC EpiScore")
plt.xlabel("Alcohol Consumption")
plt.title("EpiScore across categories")
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_usualdrinkers_notfiltered.pdf", dpi = 300)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize = (9, 3), sharey = True)
# Usual
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = lighten_color(col[0], 1.3))
axes[0].set_ylabel("AC EpiScore")
axes[0].set_xlabel("Alcohol Consumption")
axes[0].set_title("Usual drinking\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
# More
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = lighten_color(col[1], 1.3))
axes[1].set_xlabel("Alcohol Consumption")
axes[1].set_title("Drank more than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[1].set_ylabel(None)
# Less
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = lighten_color(col[2], 1.3))
axes[2].set_xlabel("Alcohol Consumption")
axes[2].set_title("Drank less than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[2].set_ylabel(None)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_separate_usualdrinkers_notfiltered.pdf", dpi = 300)
plt.close(fig)

## USUAL DRINKERS, FILTERED
col = set_colors(3, palette)
sns.set_palette(sns.color_palette(col))
true_col = "alcohol_units_log"
pred_col = "ac_pred_filt"
df = gs_target_cats
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
fig, axes = plt.subplots(1, 1, figsize = (8, 4))
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Usual (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[0], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "More (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[1], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Less (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[2], 1.3), levels = 4)
sns.despine(offset=10, trim=True);
plt.legend(frameon = False)
plt.ylabel("AC EpiScore")
plt.xlabel("Alcohol Consumption")
plt.title("EpiScore across categories")
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_usualdrinkers.pdf", dpi = 300)
plt.close(fig)


fig, axes = plt.subplots(1, 3, figsize = (9, 3), sharey = True)
# Usual
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = lighten_color(col[0], 1.3))
axes[0].set_ylabel("AC EpiScore")
axes[0].set_xlabel("Alcohol Consumption")
axes[0].set_title("Usual drinking\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
# More
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = lighten_color(col[1], 1.3))
axes[1].set_xlabel("Alcohol Consumption")
axes[1].set_title("Drank more than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[1].set_ylabel(None)
# Less
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = lighten_color(col[2], 1.3))
axes[2].set_xlabel("Alcohol Consumption")
axes[2].set_title("Drank less than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[2].set_ylabel(None)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_separate_usualdrinkers.pdf", dpi = 300)
plt.close(fig)

## EVERYONE, NOT FILTERED
col = set_colors(3, palette)
sns.set_palette(sns.color_palette(col))
true_col = "alcohol_units_log"
pred_col = "ac_pred_all"
df = gs_target_cats
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"]) 
fig, axes = plt.subplots(1, 1, figsize = (8, 4))
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Usual (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[0], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "More (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[1], 1.3), levels = 4)
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Less (r = %s)" % round(r[0], 3))
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes, color = lighten_color(col[2], 1.3), levels = 4)
sns.despine(offset=10, trim=True);
plt.legend(frameon = False)
plt.ylabel("AC EpiScore")
plt.xlabel("Alcohol Consumption")
plt.title("EpiScore across categories")
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_everyone_notfiltered.pdf", dpi = 300)
plt.close(fig)


fig, axes = plt.subplots(1, 3, figsize = (9, 3), sharey = True)
# Usual
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 2,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = lighten_color(col[0], 1.3))
axes[0].set_ylabel("AC EpiScore")
axes[0].set_xlabel("Alcohol Consumption")
axes[0].set_title("Usual drinking\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
# More
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 1,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = lighten_color(col[1], 1.3))
axes[1].set_xlabel("Alcohol Consumption")
axes[1].set_title("Drank more than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[1].set_ylabel(None)
# Less
df = gs_target_cats.loc[gs_target_cats["alcohol_usual"] == 3,]
r = scipy.stats.pearsonr(df[pred_col], df["alcohol_units"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.3, 'linewidth':0, 'edgecolor':None})
sns.kdeplot(y = pred_col, x = true_col, data = df, ax = axes[2], color = lighten_color(col[2], 1.3))
axes[2].set_xlabel("Alcohol Consumption")
axes[2].set_title("Drank less than normal\nN = %s | r = %s" % (len(df.index), round(r[0],3)))
axes[2].set_ylabel(None)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_gsw4_categories_separate_everyone_notfiltered.pdf", dpi = 300)
plt.close(fig)


#### SEX SPECIFIC MODELS - ELNET
#####################################################################################

pred_21 = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_sexspecific_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", index_col = 0)
pred_36 = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_sexspecific_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", index_col = 0)
pred_all = pd.concat([pred_21, pred_36], axis = 0)

## Scatterplot for each category (sex agnos, sex opp, same sex), tested in LBC
df = pred_all
true_col = "alcunitsupw_log"
col = set_colors(3, palette)
sns.set_palette(sns.color_palette(col))
fig, axes = plt.subplots(1, 1, figsize = (8, 4))
pred_col = "ac_pred_sexagnos"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Sex-agnostic (r = %s)" % round(r[0], 3))
pred_col = "ac_pred_samesex"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Same (r = %s)" % round(r[0], 3))
pred_col = "ac_pred_opposex"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes, color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Opposite (r = %s)" % round(r[0], 3))
sns.despine(offset=10, trim=True);
plt.legend(frameon = False)
plt.ylabel("AC EpiScore")
plt.xlabel("Alcohol Consumption")
plt.title("Sex-specific EpiScore")
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc_sexstratified.pdf", dpi = 300)
plt.close(fig)

## Scatterplot for each category (sex agnos, sex opp, same sex), tested in LBC, per cohort
true_col = "alcunitsupw_log"
col = set_colors(3, palette)
sns.set_palette(sns.color_palette(col))
fig, axes = plt.subplots(1, 2, figsize = (9, 3), sharey = True)
df = pred_21
pred_col = "ac_pred_sexagnos"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Sex-agnostic (r = %s)" % round(r[0], 3))
pred_col = "ac_pred_samesex"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Same (r = %s)" % round(r[0], 3))
pred_col = "ac_pred_opposex"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[0], color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Opposite (r = %s)" % round(r[0], 3))
axes[0].set_title("LBC1921")
axes[0].legend(frameon = False)
axes[0].set_ylabel("AC EpiScore")
axes[0].set_xlabel("Alcohol Consumption")
df = pred_36
pred_col = "ac_pred_sexagnos"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[0], line_kws={"color": lighten_color(col[0], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Sex-agnostic (r = %s)" % round(r[0], 3))
pred_col = "ac_pred_samesex"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[1], line_kws={"color": lighten_color(col[1], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Same (r = %s)" % round(r[0], 3))
pred_col = "ac_pred_opposex"
r = scipy.stats.pearsonr(df[pred_col], df["alcunitsupw"])
sns.regplot(y = pred_col, x = true_col, data = df, ax = axes[1], color = col[2], line_kws={"color": lighten_color(col[2], 1.3)}, scatter_kws={'alpha':0.5, 'linewidth':0, 'edgecolor':None, 'rasterized':True}, label = "Opposite (r = %s)" % round(r[0], 3))
axes[1].set_title("LBC1936")
axes[1].legend(frameon = False)
axes[1].set_xlabel("Alcohol Consumption")
axes[1].set_ylabel(None)
sns.despine(offset=10, trim=True);
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc_sexstratified_separatecohorts.pdf", dpi = 300)
plt.close(fig)


#### AUC HEAVY DRINKERS - ELNET
#####################################################################################

lbc_target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", index_col = 0)

# Binarize alcohol consumption, 0 = no drinkers/light drinkers, 1 = heavy drinkers
lbc_target["alcohol_cat_bi"] = [1 if i == "moderate-heavy_drinker" else 0 for i in lbc_target["alcohol_cat"]]

# Import predictions
lbc_predictions_21 = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)
lbc_predictions_36 = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", index_col = 0)

# Separate true values into 1936 and 1921
lbc_target_36 = lbc_target.loc[lbc_target["cohort"] == "LBC36"] # 2797   24
lbc_target_21 = lbc_target.loc[lbc_target["cohort"] == "LBC21"] # 692  24

# Match rownames
lbc_target_36 = lbc_target_36.loc[lbc_predictions_36.index]
lbc_target_21 = lbc_target_21.loc[lbc_predictions_21.index]

# ROC curves + AUC calc
fpr_21, tpr_21, threshold_21 = metrics.roc_curve(lbc_target_21["alcohol_cat_bi"], lbc_predictions_21["ac_pred"])
roc_auc_21 = metrics.auc(fpr_21, tpr_21) # 0.8992405444018348
fpr_36, tpr_36, threshold_36 = metrics.roc_curve(lbc_target_36["alcohol_cat_bi"], lbc_predictions_36["ac_pred"])
roc_auc_36 = metrics.auc(fpr_36, tpr_36) # 0.7777449664429531

col = set_colors(3, palette)
fig, axes = plt.subplots(1, 1, figsize = (4, 4))
axes.plot(fpr_21, tpr_21, col[0], label = 'LBC1921 AUC = %0.2f' % roc_auc_21)
axes.plot(fpr_36, tpr_36, col[1], label = 'LBC1936 AUC = %0.2f' % roc_auc_36)
sns.despine(offset=10, trim=True);
axes.set_ylabel('True Positive Rate')
axes.set_xlabel('False Positive Rate')
axes.legend(loc = 'lower right', frameon = False)
axes.set_xlim([0, 1.05])
axes.set_ylim([0, 1.05])
axes.set_ylabel('True Positive Rate')
axes.set_xlabel('False Positive Rate')
axes.set_title("Heavy Drinker Classification")
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/final_plots/prediction_lbc_heavydrinkers_auc.pdf", dpi = 300)
plt.close(fig)




