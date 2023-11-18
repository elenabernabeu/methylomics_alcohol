#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")
library("glmnet")
library("pROC")
library("PRROC")
library("precrec")

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

datadir <- "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/"
localdir <- "/Local_Data/methylation/GS_20k/Chromosomes/" # p17


## Performance across categories (non-drinkers/moderate/heavy drinkers) AUC in LBC
###################################################################################

lbc_target <- read.table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T, row.names = 1)

# Binarize alcohol consumption, 0 = no drinkers/light drinkers, 1 = heavy drinkers
lbc_target$alcohol_cat_bi <- ifelse(lbc_target$alcohol_cat == "moderate-heavy_drinker", 1, 0)

# Import predictions
lbc_predictions_36 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)
lbc_predictions_21 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)
lbc_predictions_36_sex <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_sexspecific_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", header = T)
lbc_predictions_21_sex <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_sexspecific_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", header = T)

# Import Daniel's predictions
lbc_prediction_36_daniel <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_danielpredictor_noadjustments.tsv", header = T)
lbc_prediction_21_daniel <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_danielpredictor_noadjustments.tsv", header = T)

# Separate true values into 1936 and 1921
lbc_target_36 <- lbc_target[lbc_target$cohort == "LBC36",] # 2797   24
lbc_target_21 <- lbc_target[lbc_target$cohort == "LBC21",] # 692  24

# Match rownames
lbc_target_36 <- lbc_target_36[lbc_predictions_36$basename,]
lbc_target_21 <- lbc_target_21[lbc_predictions_21$basename,]

# ROC curves - All
roc_36 <- roc(response = lbc_target_36$alcohol_cat_bi, predictor = lbc_predictions_36$ac_pred)
roc_21 <- roc(response = lbc_target_21$alcohol_cat_bi, predictor = lbc_predictions_21$ac_pred)
auc_36 <- auc(roc_36) # 0.7777
auc_21 <- auc(roc_21) # 0.8992
ci_auc_36 <- ci.auc(roc_36) # 0.7353-0.8202
ci_auc_21 <- ci.auc(roc_21) # 0.8458-0.9527

# ROC curves - Sexagnos
roc_36_sexagnos <- roc(response = lbc_target_36$alcohol_cat_bi, predictor = lbc_predictions_36_sex$ac_pred_sexagnos)
roc_21_sexagnos <- roc(response = lbc_target_21$alcohol_cat_bi, predictor = lbc_predictions_21_sex$ac_pred_sexagnos)
auc_36_sexagnos <- auc(roc_36_sexagnos) # 0.7623
auc_21_sexagnos <- auc(roc_21_sexagnos) # 0.8754
ci_auc_36_sexagnos <- ci.auc(roc_36_sexagnos) # 0.7189-0.8057
ci_auc_21_sexagnos <- ci.auc(roc_21_sexagnos) # 0.8121-0.9387

# ROC curves - Oppo sex
roc_36_opposex <- roc(response = lbc_target_36$alcohol_cat_bi, predictor = lbc_predictions_36_sex$ac_pred_opposex)
roc_21_opposex <- roc(response = lbc_target_21$alcohol_cat_bi, predictor = lbc_predictions_21_sex$ac_pred_opposex)
auc_36_opposex <- auc(roc_36_opposex) # 0.7622
auc_21_opposex <- auc(roc_21_opposex) # 0.8665
ci_auc_36_opposex <- ci.auc(roc_36_opposex) # 0.7202-0.8043
ci_auc_21_opposex <- ci.auc(roc_21_opposex) # 0.7969-0.9362

# ROC curves - Same sex
roc_36_samesex <- roc(response = lbc_target_36$alcohol_cat_bi, predictor = lbc_predictions_36_sex$ac_pred_samesex)
roc_21_samesex <- roc(response = lbc_target_21$alcohol_cat_bi, predictor = lbc_predictions_21_sex$ac_pred_samesex)
auc_36_samesex <- auc(roc_36_samesex) # 0.7768
auc_21_samesex <- auc(roc_21_samesex) # 0.8683
ci_auc_36_samesex <- ci.auc(roc_36_samesex) # 0.7344-0.8192
ci_auc_21_samesex <- ci.auc(roc_21_samesex) # 0.8073-0.9293

# ROC curves Daniel
roc_36_dan <- roc(response = lbc_target_36$alcohol_cat_bi, predictor = lbc_prediction_36_daniel$ac_pred)
roc_21_dan <- roc(response = lbc_target_21$alcohol_cat_bi, predictor = lbc_prediction_21_daniel$ac_pred)
auc_36_dan <- auc(roc_36_dan) # 0.7262
auc_21_dan <- auc(roc_21_dan) # 0.7692
ci_auc_36_dan <- ci.auc(roc_36_dan) # 0.6798-0.7727
ci_auc_21_dan <- ci.auc(roc_21_dan) # 0.6807-0.8576