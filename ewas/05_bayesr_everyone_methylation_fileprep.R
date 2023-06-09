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

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)


## (1) Pheno file prep
############################################################################

# Import prepped alcohol data
alcohol <- read.delim(paste0(datadir, "alcohol_16717.tsv"), row.names = 2)

# Import covariates to residualize pheno by: age, sex, smoking episcore
epismoker <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/epismoker_20k.csv", row.names = 1)
age_sex <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
covars <- merge(as.data.frame(epismoker[,"smokingScore",drop=FALSE]), as.data.frame(age_sex[,c("age", "sex", "Batch")]), by='row.names', all=TRUE)
covars[is.na(covars["smokingScore"]), "smokingScore"] <- mean(covars$smokingScore, na.rm = T) # Mean impute missing smoking episcores
rownames(covars) <- covars[,1]

# Residualize
alcohol["alcohol_units_log_resid"] <- resid(lm(alcohol$alcohol_units_log ~ covars[rownames(alcohol),]$age + factor(covars[rownames(alcohol),]$sex) + covars[rownames(alcohol),]$smokingScore, na.action = na.exclude)) 

# Export residualized
write.table(x = t(as.matrix(as.numeric(scale(alcohol$alcohol_units_log_resid)))),file = paste0(datadir, "allindividuals_16717_logalcohol_residualized_everyone.csvphen") ,quote = F, sep = ",", row.names = F, col.names = F)


## (2) Methylation files prep
############################################################################

meth_covars <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/covars_18413.tsv", row.names = 1)
all_covars <- data.frame(covars[rownames(alcohol),], meth_covars[rownames(alcohol), c(8:ncol(meth_covars))])
probes <- read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)

# Set up stuff for paralellizing
cores <- detectCores() # 128
cl <- makeCluster(6, outfile="parallel_test.txt") # Let's start with 6 since they use up a lot of RAM
registerDoParallel(cl)

# Iterate per chromosome
foreach(i=1:22, .packages = "data.table") %dopar% { 
  print(paste0("Working on chromosome ",i))

  # Import
  meth <- readRDS(paste0(localdir, "GS20k_chr", i, "_mvals.rds"))
  
  # Subset
  meth <- meth[,which(colnames(meth) %in% rownames(alcohol))] # Subset to those with phenotype in question 
  meth <- meth[which(rownames(meth) %in% probes$V1),] # Subset to probes passing QC 
  
  # Match order of IDs in phenotype and methylation file 
  meth <- meth[,match(rownames(alcohol), colnames(meth))]
  
  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth, 1, meanimpute)
  meth <- t(meth)
  
  # Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~as.factor(sex) + age + as.factor(Batch) + smokingScore, data=all_covars)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)),]
  rm(fit.resid)
  gc()
  
  # Scale 
  meth <- t(apply(meth,1,scale))

  # Save out residualised file scaled
  fwrite(meth, paste0(datadir, "methylation/GS20k_chr", i, "_resid_mvals_everyone_nogeno_somecovars.txt"),row.names=F)
  
  # Write out CpGs 
  cpgs <- as.data.frame(row.names(meth))
  names(cpgs)[1] <- "CpG"
  fwrite(cpgs, paste0(datadir, "methylation/GS20k_chr", i, "_cpgs.txt"),row.names=F)
  
  # Remove methylation object and clean up environment 
  rm(meth)
  gc()
} 

# End parallel
stopCluster(cl)

# Fuse methylation files - scaled
files <- list.files(path = paste0(datadir, "methylation/"), pattern = "_resid_mvals_everyone_nogeno_somecovars.txt") # Get files
files <- files[order(files)]
data <- rbindlist(lapply(paste0(datadir, "methylation/", files),fread))
gc()

# Export fused methylation file - scaled
fwrite(x = as.matrix(data), paste0(datadir, "methylation/GS20k_allchrom_everyone_nogeno_somecovars.csv"), sep = ",", row.names = F, col.names = F, quote = F) 
gc()

