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

datadir <- "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/"
localdir <- "/Local_Data/methylation/GS_20k/Chromosomes/" # p17

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)


## Import data
##########################################################################

#df <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")
df <- readRDS("/Local_Data/methylation/GS_20k/mvals.rds")
df <- t(df)
alcohol <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/alcohol_10506_usualdrinkers.tsv")
alcohol_all <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/alcohol_16717.tsv")
rownames(alcohol) <- alcohol$ID
rownames(alcohol_all) <- alcohol_all$ID
probes <- read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)$V1

# Filter meth
df <- df[alcohol_all$ID, probes]
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]
df <- df[,colnames(df) %in% rownames(common_anno)] # 386399 CpGs left
gc()

df <- m2beta(df)
gc()

df <- scale(df)
gc()


## Adjust phenos
##########################################################################

alcohol["alcohol_units_log_resid"] <- resid(lm(alcohol$alcohol_units_log ~ alcohol$age + factor(alcohol$sex), na.action = na.exclude)) 
alcohol_all["alcohol_units_log_resid"] <- resid(lm(alcohol_all$alcohol_units_log ~ alcohol_all$age + factor(alcohol_all$sex), na.action = na.exclude)) 
match_alcohol <- sample(rownames(alcohol_all), nrow(alcohol))
alcohol_all_match <- alcohol_all[match_alcohol,]


## Divide into training and testing
##########################################################################

alcohol_w1w3 <- alcohol[alcohol$Set %in% c("W1", "W3"),]
alcohol_w4 <- alcohol[alcohol$Set == "W4",]
alcohol_all_w1w3 <- alcohol_all[alcohol_all$Set %in% c("W1", "W3"),]
alcohol_all_w4 <- alcohol_all[alcohol_all$Set == "W4",]

# Match sample sizes everyone to test
match_alcohol_w1w3 <- sample(rownames(alcohol_all_w1w3), nrow(alcohol_w1w3))
match_alcohol_w4 <- sample(rownames(alcohol_all_w4), nrow(alcohol_w4))
alcohol_all_w1w3_match <- alcohol_all_w1w3[match_alcohol_w1w3,]
alcohol_all_w4_match <- alcohol_all[match_alcohol_w4,]
gc()


## Elnet
##########################################################################
seed <- 1234

## Everyone, w1w3
x <- df[rownames(alcohol_all_w1w3),]
y <- alcohol_all_w1w3$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3_noadjustments_scaledmeth.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Usual drinkers, w1w3
x <- df[rownames(alcohol_w1w3),]
y <- alcohol_w1w3$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3_noadjustments_scaledmeth.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Everyone, w1w3 (match sample sizes)
x <- df[rownames(alcohol_all_w1w3_match),]
y <- alcohol_all_w1w3_match$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_matchsamplesize_w1w3_noadjustments_scaledmeth.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Everyone, w1w3w4
x <- df[rownames(alcohol_all),]
y <- alcohol_all$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3w4_noadjustments_scaledmeth.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Usual drinkers, w1w3w4
x <- df[rownames(alcohol),]
y <- alcohol$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3w4_noadjustments_scaledmeth.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Everyone, w1w3w4
x <- df[rownames(alcohol_all_match),]
y <- alcohol_all_match$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_matchsamplesize_w1w3w4_noadjustments_scaledmeth.tsv", row.names = F, sep = "\t", quote = F)
gc()


## Elnet - FILTERED CpGs UNION
##########################################################################

ewas_cpgs_carreras <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/cpgs_ewas_carrerasgallo.txt", header = F)$V1 # 2569
ewas_cpgs_dugue <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/cpgs_ewas_dugue.txt", header = F)$V1 # 1414
ewas_cpgs_liu <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/cpgs_ewas_liu.txt", header = F)$V1 # 363

# Union
ewas_cpgs <- union(ewas_cpgs_carreras, union(ewas_cpgs_dugue, ewas_cpgs_liu)) # 3999

# Intersection with 450K
ewas_cpgs <- ewas_cpgs[ewas_cpgs %in% colnames(df)]
df_filt <- df[,ewas_cpgs]

## Everyone, w1w3
x <- df_filt[rownames(alcohol_all_w1w3),]
y <- alcohol_all_w1w3$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Usual drinkers, w1w3
x <- df_filt[rownames(alcohol_w1w3),]
y <- alcohol_w1w3$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Everyone, w1w3, matching sample sizes
x <- df_filt[rownames(alcohol_all_w1w3_match),]
y <- alcohol_all_w1w3_match$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_matchsamplesize_w1w3_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Everyone, w1w3w4
x <- df_filt[rownames(alcohol_all),]
y <- alcohol_all$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Usual drinkers, w1w3w4
x <- df_filt[rownames(alcohol),]
y <- alcohol$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
gc()

## Everyone, w1w3w4, matching sample sizes
x <- df_filt[rownames(alcohol_all_match),]
y <- alcohol_all_match$alcohol_units_log_resid
folds <- 10
cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)
gc()
coefs <- coef(fit) # Extract coeficients 
coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
names(coefs)[1] <- "Coefficient" # Tidy naming 
coefs$CpG_Site <- rownames(coefs) # Create cpg column
coefs <- coefs[c(2,1)] # order 
write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_matchsamplesize_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
gc()


## Test in w4
##########################################################################

coefs_ev <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3_noadjustments_scaledmeth.tsv", header = T)
coefs_ud <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3_noadjustments_scaledmeth.tsv", header = T)
coefs_ev_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)
coefs_ud_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_usualdrinkers_w1w3_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)
coefs_ev_filt_match <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_matchsamplesize_w1w3_noadjustments_scaledmeth_filteredcpgs.tsv", header = T)

rownames(coefs_ev) <- coefs_ev$CpG_Site
rownames(coefs_ud) <- coefs_ud$CpG_Site
rownames(coefs_ev_filt) <- coefs_ev_filt$CpG_Site
rownames(coefs_ud_filt) <- coefs_ud_filt$CpG_Site
rownames(coefs_ev_filt_match) <- coefs_ev_filt_match$CpG_Site

coefs_ni_ev <- coefs_ev[2:nrow(coefs_ev),] # 551
coefs_ni_ud <- coefs_ud[2:nrow(coefs_ud),] # 460
coefs_ni_ev_filt <- coefs_ev_filt[2:nrow(coefs_ev_filt),] # 343
coefs_ni_ud_filt <- coefs_ud_filt[2:nrow(coefs_ud_filt),] # 360
coefs_ni_ev_filt_match <- coefs_ev_filt_match[2:nrow(coefs_ev_filt_match),] # 275

# Meth
df_w4 <- df[rownames(alcohol_all_w4),]

# Predictions
pred_w4_ev <- df_w4[,rownames(coefs_ni_ev)] %*% coefs_ni_ev$Coefficient
pred_w4_ev <- as.data.frame(pred_w4_ev)
names(pred_w4_ev) <- c("ac_pred")

pred_w4_ud <- df_w4[,rownames(coefs_ni_ud)] %*% coefs_ni_ud$Coefficient
pred_w4_ud <- as.data.frame(pred_w4_ud)
names(pred_w4_ud) <- c("ac_pred")

pred_w4_ev_filt <- df_w4[,rownames(coefs_ni_ev_filt)] %*% coefs_ni_ev_filt$Coefficient
pred_w4_ev_filt <- as.data.frame(pred_w4_ev_filt)
names(pred_w4_ev_filt) <- c("ac_pred")

pred_w4_ud_filt <- df_w4[,rownames(coefs_ni_ud_filt)] %*% coefs_ni_ud_filt$Coefficient
pred_w4_ud_filt <- as.data.frame(pred_w4_ud_filt)
names(pred_w4_ud_filt) <- c("ac_pred")

pred_w4_ev_filt_match <- df_w4[,rownames(coefs_ni_ev_filt_match)] %*% coefs_ni_ev_filt_match$Coefficient
pred_w4_ev_filt_match <- as.data.frame(pred_w4_ev_filt_match)
names(pred_w4_ev_filt_match) <- c("ac_pred")

# Export
write.table(data.frame(basename = rownames(pred_w4_ev), pred_w4_ev), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_w4_ud), pred_w4_ud), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_usualdrinkerspredictor_noadjustments_scaledmeth.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_w4_ev_filt), pred_w4_ev_filt), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_w4_ud_filt), pred_w4_ud_filt), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_usualdrinkerspredictor_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_w4_ev_filt_match), pred_w4_ev_filt_match), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_matchsamplesize_noadjustments_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)

# Re-import
pred_w4_ev <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_noadjustments_scaledmeth.tsv", row.names = 1, header = T)
pred_w4_ud <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_usualdrinkerspredictor_noadjustments_scaledmeth.tsv", row.names = 1, header = T)
pred_w4_ev_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = 1, header = T)
pred_w4_ud_filt <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_usualdrinkerspredictor_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = 1, header = T)
pred_w4_ev_filt_match <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_gs_w4_everyonepredictor_matchsamplesize_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = 1, header = T)

alcohol_all_w4 <- alcohol_all_w4[rownames(pred_w4_ev),]
alcohol_all_w4$ac_pred_ev <- pred_w4_ev$ac_pred
alcohol_all_w4$ac_pred_ud <- pred_w4_ud$ac_pred
alcohol_all_w4$ac_pred_ev_filt <- pred_w4_ev_filt$ac_pred
alcohol_all_w4$ac_pred_ud_filt <- pred_w4_ud_filt$ac_pred
alcohol_all_w4$ac_pred_ev_filt_match <- pred_w4_ev_filt_match$ac_pred

# Correlations non-log
r_w4_ev <- cor(alcohol_all_w4$alcohol_units, alcohol_all_w4$ac_pred_ev, use="pairwise.complete.obs") # 0.4203013
r_w4_ud <- cor(alcohol_all_w4$alcohol_units, alcohol_all_w4$ac_pred_ud, use="pairwise.complete.obs") # 0.4036746
r_w4_ev_filt <- cor(alcohol_all_w4$alcohol_units, alcohol_all_w4$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.4447544 Union, 0.4529906 Carreras
r_w4_ud_filt <- cor(alcohol_all_w4$alcohol_units, alcohol_all_w4$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.4337525 Union, 0.4483648 Carreras
r_w4_ev_filt_match <- cor(alcohol_all_w4$alcohol_units, alcohol_all_w4$ac_pred_ev_filt_match, use="pairwise.complete.obs") # 0.442754 Union

# Correlations log
r_w4_ev <- cor(alcohol_all_w4$alcohol_units_log, alcohol_all_w4$ac_pred_ev, use="pairwise.complete.obs") # 0.4174333
r_w4_ud <- cor(alcohol_all_w4$alcohol_units_log, alcohol_all_w4$ac_pred_ud, use="pairwise.complete.obs") # 0.3963741
r_w4_ev_filt <- cor(alcohol_all_w4$alcohol_units_log, alcohol_all_w4$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.4412519
r_w4_ud_filt <- cor(alcohol_all_w4$alcohol_units_log, alcohol_all_w4$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.4198474
r_w4_ev_filt_match <- cor(alcohol_all_w4$alcohol_units_log, alcohol_all_w4$ac_pred_ev_filt_match, use="pairwise.complete.obs") # 0.4328199


# Incremental R^2?
null_w4_ev <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ev <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ev - null_w4_ev), 3) # 16.17

null_w4_ud <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ud <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ud - null_w4_ud), 3) # 15.757

null_w4_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ev_filt - null_w4_ev_filt), 3) # 18.282 Union, 18.191 Carreras

null_w4_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ud_filt - null_w4_ud_filt), 3) # 17.821 Union, 18.381 Carreras

null_w4_ev_filt_match <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ev_filt_match <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev_filt_match, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ev_filt_match - null_w4_ev_filt_match), 3) # 17.456 Union

# Incremental R^2 (log)?
null_w4_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ev - null_w4_ev), 3) # 16.213

null_w4_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ud - null_w4_ud), 3) # 15.285

null_w4_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ev_filt - null_w4_ev_filt), 3) # 18.154

null_w4_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ud_filt - null_w4_ud_filt), 3) # 16.807

null_w4_ev_filt_match <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_all_w4))$r.squared
full_w4_ev_filt_match <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev_filt_match, data=alcohol_all_w4))$r.squared
round(100*(full_w4_ev_filt_match - null_w4_ev_filt_match), 3) # 16.863


## Performance across categories (usual/frequent/less)
##########################################################################

alcohol <- alcohol_all_w4
alcohol_cats <- alcohol[!is.na(alcohol$alcohol_usual),] # 7642/8033
alcohol_usual <- alcohol_cats[alcohol_cats$alcohol_usual == 2,] # 4888/7642
alcohol_more <- alcohol_cats[alcohol_cats$alcohol_usual == 1,] # 1920/7642
alcohol_less <- alcohol_cats[alcohol_cats$alcohol_usual == 3,] # 834/7642

# Non-Log
r_usual_ud <- cor(alcohol_usual$alcohol_units, alcohol_usual$ac_pred_ud, use="pairwise.complete.obs") # 0.4606127
r_more_ud <- cor(alcohol_more$alcohol_units, alcohol_more$ac_pred_ud, use="pairwise.complete.obs") # 0.3456038
r_less_ud <- cor(alcohol_less$alcohol_units, alcohol_less$ac_pred_ud, use="pairwise.complete.obs") # 0.3213474

r_usual_ev <- cor(alcohol_usual$alcohol_units, alcohol_usual$ac_pred_ev, use="pairwise.complete.obs") # 0.47604
r_more_ev <- cor(alcohol_more$alcohol_units, alcohol_more$ac_pred_ev, use="pairwise.complete.obs") # 0.3637376
r_less_ev <- cor(alcohol_less$alcohol_units, alcohol_less$ac_pred_ev, use="pairwise.complete.obs") # 0.3328542

r_usual_ud_filt <- cor(alcohol_usual$alcohol_units, alcohol_usual$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.4900849
r_more_ud_filt <- cor(alcohol_more$alcohol_units, alcohol_more$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.3940938
r_less_ud_filt <- cor(alcohol_less$alcohol_units, alcohol_less$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.3084806

r_usual_ev_filt <- cor(alcohol_usual$alcohol_units, alcohol_usual$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.502323
r_more_ev_filt <- cor(alcohol_more$alcohol_units, alcohol_more$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.3923329
r_less_ev_filt <- cor(alcohol_less$alcohol_units, alcohol_less$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.319828

# Log
r_usual_ud <- cor(alcohol_usual$alcohol_units_log, alcohol_usual$ac_pred_ud, use="pairwise.complete.obs") # 0.4496721
r_more_ud <- cor(alcohol_more$alcohol_units_log, alcohol_more$ac_pred_ud, use="pairwise.complete.obs") # 0.3616474
r_less_ud <- cor(alcohol_less$alcohol_units_log, alcohol_less$ac_pred_ud, use="pairwise.complete.obs") # 0.3037579

r_usual_ev <- cor(alcohol_usual$alcohol_units_log, alcohol_usual$ac_pred_ev, use="pairwise.complete.obs") # 0.472635
r_more_ev <- cor(alcohol_more$alcohol_units_log, alcohol_more$ac_pred_ev, use="pairwise.complete.obs") # 0.3798706
r_less_ev <- cor(alcohol_less$alcohol_units_log, alcohol_less$ac_pred_ev, use="pairwise.complete.obs") # 0.3145444

r_usual_ud_filt <- cor(alcohol_usual$alcohol_units_log, alcohol_usual$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.4799474
r_more_ud_filt <- cor(alcohol_more$alcohol_units_log, alcohol_more$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.3894215
r_less_ud_filt <- cor(alcohol_less$alcohol_units_log, alcohol_less$ac_pred_ud_filt, use="pairwise.complete.obs") # 0.2907736

r_usual_ev_filt <- cor(alcohol_usual$alcohol_units_log, alcohol_usual$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.4969878
r_more_ev_filt <- cor(alcohol_more$alcohol_units_log, alcohol_more$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.4017776
r_less_ev_filt <- cor(alcohol_less$alcohol_units_log, alcohol_less$ac_pred_ev_filt, use="pairwise.complete.obs") # 0.3102138

# Incremental R^2? (Non-Log)
null_w4_usual_ud <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ud <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ud - null_w4_usual_ud), 3) # 19.966
null_w4_more_ud <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ud <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_more))$r.squared
round(100*(full_w4_more_ud - null_w4_more_ud), 3) # 12.068
null_w4_less_ud <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ud <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_less))$r.squared
round(100*(full_w4_less_ud - null_w4_less_ud), 3) # 10.3

null_w4_usual_ev <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ev <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ev - null_w4_usual_ev), 3) # 20.34
null_w4_more_ev <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ev <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_more))$r.squared
round(100*(full_w4_more_ev - null_w4_more_ev), 3) # 12.43
null_w4_less_ev <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ev <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_less))$r.squared
round(100*(full_w4_less_ev - null_w4_less_ev), 3) # 10.713

null_w4_usual_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ud_filt - null_w4_usual_ud_filt), 3) # 22.458
null_w4_more_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_more))$r.squared
round(100*(full_w4_more_ud_filt - null_w4_more_ud_filt), 3) # 15.103
null_w4_less_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ud_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_less))$r.squared
round(100*(full_w4_less_ud_filt - null_w4_less_ud_filt), 3) # 9.719

null_w4_usual_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ev_filt - null_w4_usual_ev_filt), 3) # 23.051
null_w4_more_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_more))$r.squared
round(100*(full_w4_more_ev_filt - null_w4_more_ev_filt), 3) # 14.603
null_w4_less_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ev_filt <- summary(lm(alcohol_units ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_less))$r.squared
round(100*(full_w4_less_ev_filt - null_w4_less_ev_filt), 3) # 10.169


# Incremental R^2? (Log)
null_w4_usual_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ud - null_w4_usual_ud), 3) # 18.972
null_w4_more_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_more))$r.squared
round(100*(full_w4_more_ud - null_w4_more_ud), 3) # 13.364
null_w4_less_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ud <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud, data=alcohol_less))$r.squared
round(100*(full_w4_less_ud - null_w4_less_ud), 3) # 9.274

null_w4_usual_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ev - null_w4_usual_ev), 3) # 20.103
null_w4_more_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_more))$r.squared
round(100*(full_w4_more_ev - null_w4_more_ev), 3) # 13.748
null_w4_less_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ev <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev, data=alcohol_less))$r.squared
round(100*(full_w4_less_ev - null_w4_less_ev), 3) # 9.615

null_w4_usual_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ud_filt - null_w4_usual_ud_filt), 3) # 21.549
null_w4_more_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_more))$r.squared
round(100*(full_w4_more_ud_filt - null_w4_more_ud_filt), 3) # 14.829
null_w4_less_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ud_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ud_filt, data=alcohol_less))$r.squared
round(100*(full_w4_less_ud_filt - null_w4_less_ud_filt), 3) # 8.69

null_w4_usual_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_usual))$r.squared
full_w4_usual_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_usual))$r.squared
round(100*(full_w4_usual_ev_filt - null_w4_usual_ev_filt), 3) # 22.631
null_w4_more_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_more))$r.squared
full_w4_more_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_more))$r.squared
round(100*(full_w4_more_ev_filt - null_w4_more_ev_filt), 3) # 15.414
null_w4_less_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex), data=alcohol_less))$r.squared
full_w4_less_ev_filt <- summary(lm(alcohol_units_log ~ age + as.factor(sex) + ac_pred_ev_filt, data=alcohol_less))$r.squared
round(100*(full_w4_less_ev_filt - null_w4_less_ev_filt), 3) # 9.621