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
library(foreign) 
library("survival")
library("ggplot2")

output_dir <- "/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/associations/"


## Import data
#################################################################################

lbc_target <- read.table("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv", header = T)

# Prep pheno
lbc_target$alcunitsupw <- as.numeric(lbc_target$alcunitsupw)
lbc_target <- lbc_target[!is.na(lbc_target$alcunitsupw),]
lbc_target$alcunitsupw_log <- log(lbc_target$alcunitsupw + 1)

# Subset
lbc_target_21 <- lbc_target[(lbc_target$WAVE == 1) & (lbc_target$cohort == "LBC21"),] # 436
lbc_target_36 <- lbc_target[(lbc_target$WAVE == 1) & (lbc_target$cohort == "LBC36"),] # 895
rownames(lbc_target_21) <- lbc_target_21$Basename
rownames(lbc_target_36) <- lbc_target_36$Basename

# LBC1921
dead_21 <- read.spss(file="/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/LBC1921_MultiOmicsOfAlcoholConsumption_RM_EB_12DEC2022.sav", to.data.frame=T)
dead_21$event <- ifelse(dead_21$dead=="ALIVE", 0, 1)
dead_21$age_event <- ifelse(dead_21$event==0, dead_21$agedaysApx_LastCensor/365.25, dead_21$agedays_death/365.25)
rownames(dead_21) <- dead_21$studyno
dead_21 <- dead_21[lbc_target_21$ID_raw,]
dead_21$Basename <- lbc_target_21$Basename
rownames(dead_21) <- dead_21$Basename

# LBC1936
dead_36 <- read.spss(file="/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/LBC1936_MultiOmicsOfAlcoholConsumption_RM_EB_13DEC2022.sav", to.data.frame=T)
dead_36$event <- 0
dead_36$event[dead_36$dead=="DEAD"] <- 1
dead_36$age_event <- ifelse(dead_36$event==0, dead_36$AgedaysApx_LastCensor/365.25, dead_36$agedays_death/365.25)
rownames(dead_36) <- dead_36$lbc36no
dead_36 <- dead_36[lbc_target_36$ID_raw,]
dead_36$Basename <- lbc_target_36$Basename
rownames(dead_36) <- dead_36$Basename

brain_36 <- read.spss(file="/Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR_EpiScore/data/LBC/AleksandraChybowska_BloodAndBrainBasedEWASofSmoking_AC_05APR2023.sav", to.data.frame=T)
rownames(brain_36) <- brain_36$lbc36no
brain_36 <- brain_36[lbc_target_36$ID_raw,]
brain_36$Basename <- lbc_target_36$Basename
rownames(brain_36) <- brain_36$Basename
brain_36 <- brain_36[,c("ICV_mm3_wX", "brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2")] # Intracreaneal volume, total brain volume, white matter volume, grey matter volume, normal appearing white matter volume
dead_36 <- cbind(dead_36, brain_36)

# Predictions
lbc21_en_usual <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1921_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = 1, header = T)
lbc36_en_usual <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/usualdrinkers_test/predictions_lbc_1936_everyonepredictor_w1w3w4_noadjustments_scaledmeth_filteredcpgs.tsv", row.names = 1, header = T)
lbc21_en_usual <- lbc21_en_usual[rownames(lbc_target_21),,drop=FALSE]
lbc36_en_usual <- lbc36_en_usual[rownames(lbc_target_36),,drop=FALSE]

lbc_target_36["acpred_en_ud"] <- lbc36_en_usual["ac_pred"]
lbc_target_21["acpred_en_ud"] <- lbc21_en_usual["ac_pred"]

lbc_target_36 <- cbind(lbc_target_36, dead_36)
lbc_target_21 <- cbind(lbc_target_21, dead_21)

lbc_target_21$tte <- lbc_target_21$age_event - lbc_target_21$age
lbc_target_36$tte <- lbc_target_36$age_event - lbc_target_36$age


## Variables of interest
#################################################################################

# Lifestyle/covariate
#common_21 <- c("nopday", "smoker", "bmi", "yrseduc", "fev", "sixmtime", "gripstr", "soccl") # Number cigarette smoked per day, current/ex/never smoker, BMI, years education, lung function, walking speed, grip strength, occupational social class
#common_36 <- c("smoknumcigs_w1", "smokcat_w1", "bmi_w1", "yrsedu_w1", "fev_w1", "sixmwk_w1", "griprh_w1", , "griplh_w1", "hmsonum_w1") # umber cigarette smoked per day, current/ex/never smoker, BMI, years education, lung function, walking speed, grip strength (right and left hand), occupational social class
common_21 <- c("nopday", "smoker", "bmi", "yrseduc", "sixmtime", "gripstr", "soccl") # Number cigarette smoked per day, current/ex/never smoker, BMI, years education, lung function, walking speed, grip strength, occupational social class
common_36 <- c("smoknumcigs_w1", "smokcat_w1", "bmi_w1", "yrsedu_w1", "sixmwk_w1", "grip", "hmsonum_w1") # umber cigarette smoked per day, current/ex/never smoker, BMI, years education, lung function, walking speed, grip strength (right and left hand), occupational social class

# Diseases (removing parkinsons in 1921 as only 2 cases and dementia from 36)
disease_21 <- c("cvhist", "crvhist", "cahist", "hyphist", "diabhist", "thyhist", "demhist", "parkinhist", "HADSD", "HADSA") # Cardiovascular, cerebrovascular, neoplasia, hypertension, diabetes, thyroid disfunction, dementia, parkinson, depressionscore, anxiety score
disease_36 <- c("cvdhist_w1", "stroke_w1", "neoplas_w1", "hibp_w1", "diab_w1", "thyroid_w1", "dement_w1", "parkin_w1", "hadsd_w1", "hadsa_w1")

# Survival
surv_21 <- c("event", "age_event", "tte") # Dead, age at death
surv_36 <- c("event", "age_event", "tte") # Dead, age at death

# Biomarkers
biom_21 <- c("chol", "triglyc")
biom_36 <- c("bld_choles_w1", "bld_triglyc_w1")

# Inflammatory proteins? (Did not ask for these)
prot_36 <- c()

# MRI variables
mri_36 <- c("ICV_mm3_wX", "brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2") # Intracreaneal volume, total brain volume, white matter volume, grey matter volume, normal appearing white matter volume


## Fix some variables
#################################################################################

# Lifestyle/covars
lbc_target_36$grip <- pmax(lbc_target_36$griplh_w1, lbc_target_36$griprh_w1) # Grip is max of left and right hands
lbc_target_36$depind_w1[lbc_target_36$depind_w1==99999] <- NA

# Social class
lbc_target_21$soccl[grepl("\\bI\\b|\\bI[a-z]\\b", lbc_target_21$soccl)] <- 1
lbc_target_21$soccl[grepl("\\bII\\b|\\bII[a-z]\\b", lbc_target_21$soccl)] <- 2
lbc_target_21$soccl[grepl("\\bIII\\b|\\bIII[a-z]\\b", lbc_target_21$soccl)] <- 3
lbc_target_21$soccl[grepl("\\bIV\\b|\\bIV[a-z]\\b", lbc_target_21$soccl)] <- 4
lbc_target_21$soccl[grepl("\\bV\\b|\\bV[a-z]\\b", lbc_target_21$soccl)] <- 5
lbc_target_21$soccl[!(lbc_target_21$soccl %in% c(1,2, 3, 4, 5))] <- NA
lbc_target_21$soccl <- as.numeric(lbc_target_21$soccl)

# Adjust lung function for age, sex, and height
# dead_36$lung <- resid(lm(fev_w1 ~ age + sex + height_w1, na.action=na.exclude, data=dead_36))

# Make smoker ordinal (0 = non-smoker, 1 = ex smoker, 2 = current smoker)
lbc_target_21$smoker <- as.character(lbc_target_21$smoker)
lbc_target_21$smoker[lbc_target_21$smoker == "no, never smoker"] <- 0
lbc_target_21$smoker[lbc_target_21$smoker == "yes, current smoker"] <- 2
lbc_target_21$smoker[lbc_target_21$smoker == "ex smoker"] <- 1

lbc_target_36$smokcat_w1 <- as.character(lbc_target_36$smokcat_w1)
lbc_target_36$smokcat_w1[lbc_target_36$smokcat_w1 == "never smoked"] <- 0
lbc_target_36$smokcat_w1[lbc_target_36$smokcat_w1 == "current smoker"] <- 2
lbc_target_36$smokcat_w1[lbc_target_36$smokcat_w1 == "ex-smoker"] <- 1

# Binarize disease histories
for (disease in disease_21[1:8]) {
	print(disease)
	lbc_target_21[,disease] <- as.character(lbc_target_21[,disease])
	lbc_target_21[,disease][grepl("yes|ye|angina", lbc_target_21[,disease], ignore.case = T)] <- 1
	lbc_target_21[,disease][grepl("no|maybe|unsure|\\?", lbc_target_21[,disease], ignore.case = T)] <- 0
	lbc_target_21[,disease][grepl(NA, lbc_target_21[,disease], ignore.case = T)] <- NA
	lbc_target_21[,disease] <- as.numeric(lbc_target_21[,disease])
	print(unique(lbc_target_21[,disease]))
}

# Binarize disease histories
for (disease in disease_36[1:8]) {
	print(disease)
	lbc_target_36[,disease] <- as.character(lbc_target_36[,disease])
	lbc_target_36[,disease][grepl("yes|ye|angina", lbc_target_36[,disease], ignore.case = T)] <- 1
	lbc_target_36[,disease][grepl("no|maybe|unsure|\\?", lbc_target_36[,disease], ignore.case = T)] <- 0
	lbc_target_36[,disease][grepl(NA, lbc_target_36[,disease], ignore.case = T)] <- NA
	lbc_target_36[,disease] <- as.numeric(lbc_target_36[,disease])
	print(unique(lbc_target_36[,disease]))
}


## Associations to lifestyle/covariates 
#################################################################################

### Regression of covariates on EpiScore
write.table(summary(lm(scale(acpred_en_ud) ~ scale(age) + as.factor(sex) + scale(nopday) + as.factor(smoker) + scale(bmi) + scale(yrseduc) + scale(as.numeric(fev)) + scale(as.numeric(sixmtime)) + scale(as.numeric(gripstr)) + scale(soccl), data = lbc_target_21))$coefficients[,c(1,2,4)], paste0(output_dir, "lifestyle_assocs_episcore_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(summary(lm(scale(alcunitsupw) ~ scale(age) + as.factor(sex) + scale(nopday) + as.factor(smoker) + scale(bmi) + scale(yrseduc) + scale(as.numeric(fev)) + scale(as.numeric(sixmtime)) + scale(as.numeric(gripstr)) + scale(soccl), data = lbc_target_21))$coefficients[,c(1,2,4)], paste0(output_dir, "lifestyle_assocs_measuredalc_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")

write.table(summary(lm(scale(acpred_en_ud) ~ scale(age) + as.factor(sex) + scale(smoknumcigs_w1) + as.factor(smokcat_w1) + scale(bmi_w1) + scale(yrsedu_w1) + scale(as.numeric(fev_w1)) + scale(as.numeric(sixmwk_w1)) + scale(as.numeric(grip)) + scale(hmsonum_w1), data = lbc_target_36))$coefficients[,c(1,2,4)], paste0(output_dir, "lifestyle_assocs_episcore_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(summary(lm(scale(alcunitsupw) ~ scale(age) + as.factor(sex) + scale(smoknumcigs_w1) + as.factor(smokcat_w1) + scale(bmi_w1) + scale(yrsedu_w1) + scale(as.numeric(fev_w1)) + scale(as.numeric(sixmwk_w1)) + scale(as.numeric(grip)) + scale(hmsonum_w1), data = lbc_target_36))$coefficients[,c(1,2,4)], paste0(output_dir, "lifestyle_assocs_measuredalc_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")

## One at a time
df_21 <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
df_21_upw <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
for (trait in common_21) {
	print(trait)
	trait_df <- data.frame("trait" = lbc_target_21[,trait], lbc_target_21[,c("age", "sex", "acpred_en_ud", "alcunitsupw")])
	df_trait_21_epi <- summary(lm(scale(as.numeric(trait)) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = trait_df))$coefficients[,c(1,2,4)]
	df_trait_21_upw <- summary(lm(scale(as.numeric(trait)) ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = trait_df))$coefficients[,c(1,2,4)]	
	print(nobs(lm(scale(as.numeric(trait)) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = trait_df)))
	df_21[trait,] <- df_trait_21_epi[nrow(df_trait_21_epi),]
	df_21_upw[trait,] <- df_trait_21_upw[nrow(df_trait_21_upw),]
}
write.table(df_21[2:nrow(df_21),], paste0(output_dir, "trait_assocs_episcore_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(df_21_upw[2:nrow(df_21_upw),], paste0(output_dir, "trait_assocs_measuredalc_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")

df_36 <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
df_36_upw <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
for (trait in common_36) {
	print(trait)
	trait_df <- data.frame("trait" = lbc_target_36[,trait], lbc_target_36[,c("age", "sex", "acpred_en_ud", "alcunitsupw")])
	df_trait_36_epi <- summary(lm(scale(as.numeric(trait)) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = trait_df))$coefficients[,c(1,2,4)]
	df_trait_36_upw <- summary(lm(scale(as.numeric(trait)) ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = trait_df))$coefficients[,c(1,2,4)]	
	print(nobs(lm(scale(as.numeric(trait)) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = trait_df)))
	df_36[trait,] <- df_trait_36_epi[nrow(df_trait_36_epi),]
	df_36_upw[trait,] <- df_trait_36_upw[nrow(df_trait_36_upw),]
}
write.table(df_36[2:nrow(df_36),], paste0(output_dir, "trait_assocs_episcore_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(df_36_upw[2:nrow(df_36_upw),], paste0(output_dir, "trait_assocs_measuredalc_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")


## Associations to disease (one per disease)
#################################################################################

df_21 <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
df_21_upw <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
for (disease in disease_21) {
	print(disease)
	disease_df <- data.frame("disease" = lbc_target_21[,disease], lbc_target_21[,c("age", "sex", "acpred_en_ud", "alcunitsupw")])
	if (endsWith(disease, "hist")) {
		df_disease_21_epi <- summary(glm(disease ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df, family = binomial))$coefficients[,c(1,2,4)]
		df_disease_21_upw <- summary(glm(disease ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = disease_df, family = binomial))$coefficients[,c(1,2,4)]
		print(nobs(glm(disease ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df, family = binomial)))
		print(sum(disease_df$disease, na.rm = T))
	} else {
		df_disease_21_epi <- summary(lm(scale(disease) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df))$coefficients[,c(1,2,4)]
		df_disease_21_upw <- summary(lm(scale(disease) ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = disease_df))$coefficients[,c(1,2,4)]		
		print(nobs(lm(scale(disease) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df)))
	}
	df_21[disease,] <- df_disease_21_epi[nrow(df_disease_21_epi),]
	df_21_upw[disease,] <- df_disease_21_upw[nrow(df_disease_21_upw),]
}
write.table(df_21[2:nrow(df_21),], paste0(output_dir, "disease_assocs_episcore_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(df_21_upw[2:nrow(df_21_upw),], paste0(output_dir, "disease_assocs_measuredalc_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")

df_36 <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
df_36_upw <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
for (disease in disease_36) {
	print(disease)
	disease_df <- data.frame("disease" = lbc_target_36[,disease], lbc_target_36[,c("age", "sex", "acpred_en_ud", "alcunitsupw")])
	if (!(disease %in% c("hadsd_w1", "hadsa_w1"))) {
		df_disease_36_epi <- summary(glm(disease ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df, family = binomial))$coefficients[,c(1,2,4)]
		df_disease_36_upw <- summary(glm(disease ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = disease_df, family = binomial))$coefficients[,c(1,2,4)]
		print(nobs(glm(disease ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df, family = binomial)))
		print(sum(disease_df$disease, na.rm = T))
	} else {
		df_disease_36_epi <- summary(lm(scale(disease) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df))$coefficients[,c(1,2,4)]
		df_disease_36_upw <- summary(lm(scale(disease) ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = disease_df))$coefficients[,c(1,2,4)]		
		print(nobs(lm(scale(disease) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = disease_df)))
	}
	df_36[disease,] <- df_disease_36_epi[nrow(df_disease_36_epi),]
	df_36_upw[disease,] <- df_disease_36_upw[nrow(df_disease_36_upw),]
}
write.table(df_36[2:nrow(df_36),], paste0(output_dir, "disease_assocs_episcore_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(df_36_upw[2:nrow(df_36_upw),], paste0(output_dir, "disease_assocs_measuredalc_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")


## Associations to biomarkers (linear model, one per marker)
#################################################################################

df_21 <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
df_21_upw <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
for (biomarker in biom_21) {
	print(biomarker)
	biom_df <- data.frame("biomarker" = lbc_target_21[,biomarker], lbc_target_21[,c("age", "sex", "acpred_en_ud", "alcunitsupw")])
	df_biom_21_epi <- summary(lm(scale(biomarker) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = biom_df))$coefficients[,c(1,2,4)]
	df_biom_21_upw <- summary(lm(scale(biomarker) ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = biom_df))$coefficients[,c(1,2,4)]		
	df_21[biomarker,] <- df_biom_21_epi[nrow(df_biom_21_epi),]
	df_21_upw[biomarker,] <- df_biom_21_upw[nrow(df_biom_21_upw),]
	print(nobs(lm(scale(biomarker) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = biom_df)))
}

write.table(df_21[2:nrow(df_21),], paste0(output_dir, "biomarker_assocs_episcore_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(df_21_upw[2:nrow(df_21_upw),], paste0(output_dir, "biomarker_assocs_measuredalc_lbc1921.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")

df_36 <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
df_36_upw <- data.frame("Beta" = 0, "SE" = 0, "p" = 0)
for (biomarker in biom_36) {
	print(biomarker)
	biom_df <- data.frame("biomarker" = lbc_target_36[,biomarker], lbc_target_36[,c("age", "sex", "acpred_en_ud", "alcunitsupw")])
	df_biom_36_epi <- summary(lm(scale(biomarker) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = biom_df))$coefficients[,c(1,2,4)]
	df_biom_36_upw <- summary(lm(scale(biomarker) ~ scale(age) + as.factor(sex) + scale(alcunitsupw), data = biom_df))$coefficients[,c(1,2,4)]		
	df_36[biomarker,] <- df_biom_36_epi[nrow(df_biom_36_epi),]
	df_36_upw[biomarker,] <- df_biom_36_upw[nrow(df_biom_36_upw),]
	print(nobs(lm(scale(biomarker) ~ scale(age) + as.factor(sex) + scale(acpred_en_ud), data = biom_df)))
}

write.table(df_36[2:nrow(df_36),], paste0(output_dir, "biomarker_assocs_episcore_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(df_36_upw[2:nrow(df_36_upw),], paste0(output_dir, "biomarker_assocs_measuredalc_lbc1936.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")


## Associations to survival (Cox models)
#################################################################################

survival_21 <- coxph(Surv(tte, event) ~ age + factor(sex) + scale(acpred_en_ud), data=lbc_target_21)
survival_21_upw <- coxph(Surv(tte, event) ~ age + factor(sex) + scale(alcunitsupw), data=lbc_target_21)

survival_36 <- coxph(Surv(tte, event) ~ age + factor(sex) + scale(acpred_en_ud), data=lbc_target_36)
survival_36_upw <- coxph(Surv(tte, event) ~ age + factor(sex) + scale(alcunitsupw), data=lbc_target_36)


## Associations to MRI variables (linear modes, one per variable) - based on Ola's --> Maybe project on wave 2?
#################################################################################

## Ola's
neuroimaging <- lbc_target_36[c("ID",  "age", "sex", "acpred_en_ud", "alcunitsupw", "ICV_mm3_wX", "brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2")]
neuroimaging <- na.omit(neuroimaging) # 533  10

columns <- c("brain_mm3_w2", "wmh_mm3_w2", "gm_mm3_w2", "nawm_mm3_w2", "ICV_mm3_wX", "acpred_en_ud", "alcunitsupw")

pdf(paste0(output_dir, "outcome_distribution_neuroimaging.pdf"))
par(mfrow = c(3, 2)) vgtfc
for(col in columns) {
  hist(neuroimaging[,col], main = col, breaks=20, xlab="")
}
dev.off()

print("Scaling variables, not tranforming them - watch out for wmh.")
neuroimaging$brain_mm3_w2_sc_t = scale(neuroimaging$brain_mm3_w2)
neuroimaging$wmh_mm3_w2_sc_t = scale(transform(neuroimaging$wmh_mm3_w2))
neuroimaging$gm_mm3_w2_sc_t = scale(neuroimaging$gm_mm3_w2)
neuroimaging$nawm_mm3_w2_sc_t = scale(neuroimaging$nawm_mm3_w2)
neuroimaging$ICV_mm3_wX_sc_t = scale(neuroimaging$ICV_mm3_wX)
neuroimaging$acpred_en_ud_sc_t = scale(neuroimaging$acpred_en_ud)
neuroimaging$alcunitsupw_sc_t = scale(neuroimaging$alcunitsupw)

pdf(paste0(output_dir, "outcome_distribution_neuroimaging_scaled.pdf"))
par(mfrow = c(3, 2))
for(col in columns) {
  hist(neuroimaging[,paste0(col, "_sc_t")], main = paste0(col, "_sc_t"), breaks=20, xlab="")
}
dev.off()

out_cont <- data.frame(Outcome=character(), 
	Predictor=character(),
	n=double(),
	Beta=double(),
	SE=double(),
	P=double(), 
	LCI=double(),
	UCI=double(),
	stringsAsFactors=FALSE)

my_data = neuroimaging
y <- c("brain_mm3_w2_sc_t", "wmh_mm3_w2_sc_t", "gm_mm3_w2_sc_t", "nawm_mm3_w2_sc_t")
x <- c("acpred_en_ud_sc_t", "alcunitsupw_sc_t") 

count <- 0

for(i in 1:length(y)) {
	for(j in 1:length(x)) {
	count = count + 1
	
	outcome <- y[i]
	predictor <- x[j]
	
	model <- lm(my_data[,y[i]] ~ my_data[,x[j]] + scale(age) + as.factor(sex) + ICV_mm3_wX_sc_t, data=my_data) # to be edited manually for covariates 
	coefs <- summary(model)$coefficients[2,c(1,2,4)] # need to make sure that the coefficient of interest is first in the model
	int <- confint(model)[2,1:2]
	n = nobs(model)

	out_cont[count, ] <- c(outcome, predictor, n, signif(coefs, 2), signif(int, 2)) 
	}
}

out_cont <- as.data.frame(out_cont)
out_cont$Beta <- as.numeric(out_cont$Beta)
out_cont$SE <- as.numeric(out_cont$SE)
out_cont$LCI <- as.numeric(out_cont$LCI)
out_cont$UCI <- as.numeric(out_cont$UCI)
write.table(out_cont, paste0(output_dir, "brainmri_assocs_lbc1936.tsv"), sep = "\t", row.names = F, quote = F)


## Prep plot
####################################################################
out_cont$Outcome = gsub("brain_mm3_w2_sc_t", "Total brain\nvolume", out_cont$Outcome)
out_cont$Outcome = gsub("wmh_mm3_w2_sc_t", "White matter\nhyperintensity\nvolume", out_cont$Outcome)
out_cont$Outcome = gsub("gm_mm3_w2_sc_t", "Gray matter\nvolume", out_cont$Outcome)
out_cont$Outcome = gsub("nawm_mm3_w2_sc_t", "Normal appearing\nwhite matter\nvolume", out_cont$Outcome)
out_cont$Predictor = gsub("acpred_en_ud_sc_t", "EpiScore", out_cont$Predictor)
out_cont$Predictor = gsub("alcunitsupw_sc_t", "Alcohol UPW", out_cont$Predictor)

My_Theme = theme(
  panel.border = element_rect(colour="black",size=1, fill = NA),
  axis.title.x = element_text(size = 16), # controls HR label size 
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 16, face = "bold"),
  legend.text=element_text(size=16),
	    legend.position=c(0.8,0.1),
	    legend.background = element_rect(size=0.5, linetype="solid", colour="black"),
	    legend.title=element_blank(),
  axis.title=element_text(size=16))

stacked = ggplot(out_cont,aes(y=Beta, x=Outcome, group=Predictor, colour=Predictor)) + 
  geom_point(size = 2, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1)+
  ylab("Effect Size (Beta)")+ 
  xlab ("") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5), axis.text.y = element_text(size = 8), legend.position = "right",
        plot.title = element_text(size = 8))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + My_Theme

pdf(paste0(output_dir, "forestplot_neuroimaging.pdf"))
stacked
dev.off()