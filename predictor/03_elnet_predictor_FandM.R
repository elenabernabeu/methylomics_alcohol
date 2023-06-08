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

df <- readRDS("/Local_Data/methylation/GS_20k/mvals.rds")
alcohol <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/alcohol_16717.tsv")
rownames(alcohol) <- alcohol$ID
probes <- read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)$V1

# Filter and process meth
df <- df[probes, alcohol$ID]
df <- t(df)
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]
df <- df[,colnames(df) %in% rownames(common_anno)] # 386399 CpGs left
df <- m2beta(df)
gc()

# Filter to EWAS CpGs
ewas_cpgs_carreras <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/cpgs_ewas_carrerasgallo.txt", header = F)$V1 # 2569
ewas_cpgs_dugue <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/cpgs_ewas_dugue.txt", header = F)$V1 # 1414
ewas_cpgs_liu <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/cpgs_ewas_liu.txt", header = F)$V1 # 363
ewas_cpgs <- union(ewas_cpgs_carreras, union(ewas_cpgs_dugue, ewas_cpgs_liu)) # 3999
ewas_cpgs <- ewas_cpgs[ewas_cpgs %in% colnames(df)]
df <- df[,ewas_cpgs]

# Separate into M and F
alcohol_F <- alcohol[alcohol$sex == "F",]
alcohol_M <- alcohol[alcohol$sex == "M",] 

# Match sample sizes to 6958 (max number of males is 6958), with sex agnos being 3479 of each sex
alcohol_F <- alcohol_F[sample(nrow(alcohol_F), 6958),]
alcohol_M <- alcohol_M[sample(nrow(alcohol_M), 6958),]
alcohol <- rbind(alcohol_F[sample(nrow(alcohol_F), 3479),], alcohol_M[sample(nrow(alcohol_M), 3479),])

df_F <- df[rownames(alcohol_F),]
df_M <- df[rownames(alcohol_M),]
df <- df[rownames(alcohol),]
gc()

df_F <- scale(df_F)
df_M <- scale(df_M)
df <- scale(df)
gc()


## Adjust phenos
##########################################################################

alcohol_F["alcohol_units_log_resid"] <- resid(lm(alcohol_F$alcohol_units_log ~ alcohol_F$age, na.action = na.exclude)) 
alcohol_M["alcohol_units_log_resid"] <- resid(lm(alcohol_M$alcohol_units_log ~ alcohol_M$age, na.action = na.exclude)) 
alcohol["alcohol_units_log_resid"] <- resid(lm(alcohol$alcohol_units_log ~ alcohol$age + alcohol$sex, na.action = na.exclude)) 


## Elnet
##########################################################################

seed <- 1234
folds <- 10

x <- df
y <- alcohol$alcohol_units_log_resid
x_F <- df_F
y_F <- alcohol_F$alcohol_units_log_resid
x_M <- df_M
y_M <- alcohol_M$alcohol_units_log_resid

cv <- cv.glmnet(x, y, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit <- glmnet(x, y, family = "gaussian", alpha = 0.5, lambda = cv$lambda.min)

cv_F <- cv.glmnet(x_F, y_F, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit_F <- glmnet(x_F, y_F, family = "gaussian", alpha = 0.5, lambda = cv_F$lambda.min)

cv_M <- cv.glmnet(x_M, y_M, family = "gaussian", alpha = 0.5, nfolds = folds, seed = seed) 
fit_M <- glmnet(x_M, y_M, family = "gaussian", alpha = 0.5, lambda = cv_M$lambda.min)

coefs <- coef(fit) # Extract coefficients 
coefs_F <- coef(fit_F) # Extract coefficients 
coefs_M <- coef(fit_M) # Extract coefficients 

coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coefficients that are 0
coefs_F <- as.data.frame(coefs_F[which(coefs_F!=0),]) # Remove coefficients that are 0 
coefs_M <- as.data.frame(coefs_M[which(coefs_M!=0),]) # Remove coefficients that are 0 

names(coefs)[1] <- "Coefficient" # Tidy naming 
names(coefs_F)[1] <- "Coefficient" # Tidy naming 
names(coefs_M)[1] <- "Coefficient" # Tidy naming 

coefs$CpG_Site <- rownames(coefs)
coefs_F$CpG_Site <- rownames(coefs_F) # Create cpg column
coefs_M$CpG_Site <- rownames(coefs_M) # Create cpg column

coefs <- coefs[c(2,1)] # 423
coefs_F <- coefs_F[c(2,1)] # 445
coefs_M <- coefs_M[c(2,1)] # 370

write.table(coefs, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
write.table(coefs_F, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_F_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)
write.table(coefs_M, "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_M_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", row.names = F, sep = "\t", quote = F)


## Predictions in LBC
####################################################################

coefs_F <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_F_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", header = T)
coefs_M <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_M_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", header = T)
coefs_all <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictor_everyone_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", header = T)

rownames(coefs_F) <- coefs_F$CpG_Site
rownames(coefs_M) <- coefs_M$CpG_Site
rownames(coefs_all) <- coefs_all$CpG_Site

coefs_ni_F <- coefs_F[2:nrow(coefs_F),]
coefs_ni_M <- coefs_M[2:nrow(coefs_M),]
coefs_ni_all <- coefs_all[2:nrow(coefs_all),]
#intercept <- coefs[1,2]
#intercept_all <- coefs_all[1,2]

# Import LBC data
lbc_target <- read.table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", sep = "\t", header = T, row.names = 1)
lbc_mvals <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_mvals_3489.rds")
lbc_mvals <- t(lbc_mvals)
lbc_mvals <- m2beta(lbc_mvals)

# Prep pheno
lbc_target$alcunitsupw <- as.numeric(lbc_target$alcunitsupw)
lbc_target <- lbc_target[!is.na(lbc_target$alcunitsupw),]
lbc_target$alcunitsupw_log <- log(lbc_target$alcunitsupw + 1)

# Replace missing values with 0
lbc_mvals[is.na(lbc_mvals)] <- 0

# Divide LBC data into LBC21 and LBC36
lbc_target_36 <- lbc_target[lbc_target["cohort"] == "LBC36",] # 2796, waves 1, 2, 3, and 4
lbc_target_21 <- lbc_target[lbc_target["cohort"] == "LBC21",] # 692, waves 1, 3, and 4

# Keep just wave 1
lbc_target_21 <- lbc_target_21[lbc_target_21$WAVE == 1,] # 436
lbc_target_36 <- lbc_target_36[lbc_target_36$WAVE == 1,] # 895

# Match CpGs
coefs_ni_lbc_F <- coefs_ni_F[which(rownames(coefs_ni_F) %in% colnames(lbc_mvals)),] # 645/672
coefs_ni_lbc_M <- coefs_ni_M[which(rownames(coefs_ni_M) %in% colnames(lbc_mvals)),] # 392/406
coefs_ni_lbc_all <- coefs_ni_all[which(rownames(coefs_ni_all) %in% colnames(lbc_mvals)),] # 1226/1269

lbc_mvals_F <- lbc_mvals[,rownames(coefs_ni_lbc_F)]
lbc_mvals_M <- lbc_mvals[,rownames(coefs_ni_lbc_M)]
lbc_mvals_all <- lbc_mvals[,rownames(coefs_ni_lbc_all)]

# Subset
lbc_mvals_36_F <- scale(lbc_mvals_F[rownames(lbc_target_36),])
lbc_mvals_21_F <- scale(lbc_mvals_F[rownames(lbc_target_21),])
lbc_mvals_36_M <- scale(lbc_mvals_M[rownames(lbc_target_36),])
lbc_mvals_21_M <- scale(lbc_mvals_M[rownames(lbc_target_21),])
lbc_mvals_36_all <- scale(lbc_mvals_all[rownames(lbc_target_36),])
lbc_mvals_21_all <- scale(lbc_mvals_all[rownames(lbc_target_21),])

# Predictions
pred_36_F <- lbc_mvals_36_F %*% coefs_ni_lbc_F$Coefficient
pred_36_F <- as.data.frame(pred_36_F)
names(pred_36_F) <- c("ac_pred")

pred_36_M <- lbc_mvals_36_M %*% coefs_ni_lbc_M$Coefficient
pred_36_M <- as.data.frame(pred_36_M)
names(pred_36_M) <- c("ac_pred")

pred_36_all <- lbc_mvals_36_all %*% coefs_ni_lbc_all$Coefficient
pred_36_all <- as.data.frame(pred_36_all)
names(pred_36_all) <- c("ac_pred")

pred_21_F <- lbc_mvals_21_F %*% coefs_ni_lbc_F$Coefficient
pred_21_F <- as.data.frame(pred_21_F)
names(pred_21_F) <- c("ac_pred")

pred_21_M <- lbc_mvals_21_M %*% coefs_ni_lbc_M$Coefficient
pred_21_M <- as.data.frame(pred_21_M)
names(pred_21_M) <- c("ac_pred")

pred_21_all <- lbc_mvals_21_all %*% coefs_ni_lbc_all$Coefficient
pred_21_all <- as.data.frame(pred_21_all)
names(pred_21_all) <- c("ac_pred")

# Export
write.table(data.frame(basename = rownames(pred_36_F), pred_36_F), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonefemalespredictor_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21_F), pred_21_F), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonefemalespredictor_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_36_M), pred_36_M), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonemalespredictor_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(pred_21_M), pred_21_M), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonemalespredictor_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)

# Add stuff to data frame
lbc_target_36["ac_pred_sexagnos"] <- pred_36_all["ac_pred"]
lbc_target_36[lbc_target_36["sex"] == "F", "ac_pred_samesex"] <- pred_36_F[lbc_target_36["sex"] == "F", "ac_pred"]
lbc_target_36[lbc_target_36["sex"] == "M", "ac_pred_samesex"] <- pred_36_M[lbc_target_36["sex"] == "M", "ac_pred"]
lbc_target_36[lbc_target_36["sex"] == "F", "ac_pred_opposex"] <- pred_36_M[lbc_target_36["sex"] == "F", "ac_pred"]
lbc_target_36[lbc_target_36["sex"] == "M", "ac_pred_opposex"] <- pred_36_F[lbc_target_36["sex"] == "M", "ac_pred"]

lbc_target_21["ac_pred_sexagnos"] <- pred_21_all["ac_pred"]
lbc_target_21[lbc_target_21["sex"] == "F", "ac_pred_samesex"] <- pred_21_F[lbc_target_21["sex"] == "F", "ac_pred"]
lbc_target_21[lbc_target_21["sex"] == "M", "ac_pred_samesex"] <- pred_21_M[lbc_target_21["sex"] == "M", "ac_pred"]
lbc_target_21[lbc_target_21["sex"] == "F", "ac_pred_opposex"] <- pred_21_M[lbc_target_21["sex"] == "F", "ac_pred"]
lbc_target_21[lbc_target_21["sex"] == "M", "ac_pred_opposex"] <- pred_21_F[lbc_target_21["sex"] == "M", "ac_pred"]

# Export
write.table(data.frame(basename = rownames(lbc_target_21), lbc_target_21), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_sexspecific_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)
write.table(data.frame(basename = rownames(lbc_target_36), lbc_target_36), "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_sexspecific_w1w3w4_noadjustments_subset6958_scaledmeth_filteredcpgs.tsv", sep = "\t", row.names = F, quote = F)

lbc_target_21 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_usualdrinkerspredictor_sexspecific_w1w3w4_noadjustments_subset4410_scaledmeth.tsv", header = T, row.names = 1)
lbc_target_36 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_usualdrinkerspredictor_sexspecific_w1w3w4_noadjustments_subset4410_scaledmeth.tsv", header = T, row.names = 1)


## Test performance
####################################################################

# Correlation
r_21_sexagnos <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_sexagnos, use="pairwise.complete.obs") # 0.4123318
r_21_oppo <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_opposex, use="pairwise.complete.obs") # 0.4161835
r_21_same <- cor(lbc_target_21$alcunitsupw, lbc_target_21$ac_pred_samesex, use="pairwise.complete.obs") # 0.4072733

r_36_sexagnos <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_sexagnos, use="pairwise.complete.obs") # 0.4488094
r_36_oppo <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_opposex, use="pairwise.complete.obs") # 0.4534403
r_36_same <- cor(lbc_target_36$alcunitsupw, lbc_target_36$ac_pred_samesex, use="pairwise.complete.obs") # 0.4780237

# Incremental DNAm R2
null_21 <- summary(lm(alcunitsupw ~ age , data=lbc_target_21))$r.squared
full_21_sexagnos <- summary(lm(alcunitsupw ~ age + sex + ac_pred_sexagnos, data=lbc_target_21))$r.squared
full_21_oppo <- summary(lm(alcunitsupw ~ age + sex + ac_pred_opposex, data=lbc_target_21))$r.squared
full_21_same <- summary(lm(alcunitsupw ~ age + sex + ac_pred_samesex, data=lbc_target_21))$r.squared
round(100*(full_21_sexagnos - null_21), 3) # 24.184
round(100*(full_21_oppo - null_21), 3) # 23.425
round(100*(full_21_same - null_21), 3) # 23.227

null_36 <- summary(lm(alcunitsupw ~ age + sex, data=lbc_target_36))$r.squared
full_36_sexagnos <- summary(lm(alcunitsupw ~ age + sex + ac_pred_sexagnos, data=lbc_target_36))$r.squared
full_36_oppo <- summary(lm(alcunitsupw ~ age + sex + ac_pred_opposex, data=lbc_target_36))$r.squared
full_36_same <- summary(lm(alcunitsupw ~ age + sex + ac_pred_samesex, data=lbc_target_36))$r.squared
round(100*(full_36_sexagnos - null_36), 3) # 19.246
round(100*(full_36_oppo - null_36), 3) # 17.988
round(100*(full_36_same - null_36), 3) # 20.026

# P-vals
summary(lm(alcunitsupw ~ age + sex, data=lbc_target_21)) 
summary(lm(alcunitsupw ~ age + sex + ac_pred_sexagnos, data=lbc_target_21)) 
summary(lm(alcunitsupw ~ age + sex + ac_pred_opposex, data=lbc_target_21)) 
summary(lm(alcunitsupw ~ age + sex + ac_pred_samesex, data=lbc_target_21)) 

summary(lm(alcunitsupw ~ age + sex, data=lbc_target_36)) 
summary(lm(alcunitsupw ~ age + sex + ac_pred_sexagnos, data=lbc_target_36)) 
summary(lm(alcunitsupw ~ age + sex + ac_pred_opposex, data=lbc_target_36)) 
summary(lm(alcunitsupw ~ age + sex + ac_pred_samesex, data=lbc_target_36)) 