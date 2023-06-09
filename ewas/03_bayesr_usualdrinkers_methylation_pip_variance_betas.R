#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena (based on Rob's and Daniel's scripts)

library("data.table")
library("foreach")
library("doParallel")


# (1) Thin out files
###################################################################

setwd("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/")

## Thin out betas
beta_everyone <- fread("beta/allindividuals_10506_logalcohol_residualized_usualdrinkers.csv", sep = ",")
beta_everyone <- as.data.frame(beta_everyone)
beta_everyone <- beta_everyone[4901:nrow(beta_everyone),]
beta_everyone <- beta_everyone[seq(1, nrow(beta_everyone), 5),]
fwrite(beta_everyone, "beta/allindividuals_10506_logalcohol_residualized_usualdrinkers_processed.csv", sep = ",", row.names = F, col.names = F)

## Thin out comp
comp_everyone <- fread("comp/allindividuals_10506_logalcohol_residualized_usualdrinkers.csv", sep = ",")
comp_everyone <- as.data.frame(comp_everyone)
comp_everyone <- comp_everyone[4901:nrow(comp_everyone),]
comp_everyone <- comp_everyone[seq(1, nrow(comp_everyone), 5),]
fwrite(comp_everyone, "comp/allindividuals_10506_logalcohol_residualized_usualdrinkers_processed.csv", sep = ",", row.names = F, col.names = F)

## Thin out sigma
sigma_everyone <- fread("sigma/allindividuals_10506_logalcohol_residualized_usualdrinkers.csv", sep = ",")
sigma_everyone <- as.data.frame(sigma_everyone)
sigma_everyone <- sigma_everyone[4901:nrow(sigma_everyone),]
sigma_everyone <- sigma_everyone[seq(1, nrow(sigma_everyone), 5),]
fwrite(sigma_everyone, "sigma/allindividuals_10506_logalcohol_residualized_usualdrinkers_processed.csv", sep = ",", row.names = F, col.names = F)


# (2) Calculate AC variance explained by epigenetic probes
###################################################################

setwd("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/")
loop <- list.files("comp/", pattern = "_usualdrinkers_processed.csv")
names <- read.csv(paste0("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/methylation/GS20k_allchrom_cpg_list.txt"),header=T) # Same for the F and M files so just need this one for all 3 runs

for(file in loop){  

  output <- matrix(nrow = 1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Model"
  output$Model <- gsub("_processed.csv*", "", file)
  
  ### (1) Mean variance explained by epigenetic probes and credible intervals
  sigma <- read.csv(paste("sigma/", file, sep = ""))  
  output$Epigenetic_Mean_Variance_Explained <- mean(sigma[,2]/rowSums(sigma[,1:2]), na.rm = T)
  output$Epigenetic_Low_CI <- quantile(sigma[,2]/rowSums(sigma[,1:2]), na.rm = T, prob = 0.025)
  output$Epigenetic_High_CI <- quantile(sigma[,2]/rowSums(sigma[,1:2]), na.rm = T, prob = 0.975)
  
  ### (2) Calculate the proportion of variance that is attributable to small, medium and large effects - 1,2,3
  betas <- fread(paste("beta/", file, sep = ""))
  betas <- as.data.frame(betas)
  comp <- fread(paste("comp/", file, sep = ""))  
  comp <- as.data.frame(comp) 
  names(comp) <- names$Marker
  names(betas) <- names$Marker
  list <- apply(comp, 1, function(x) which(!x %in% c(1,2,3)))  
  
  x <- as.matrix(0, ncol = 1, nrow = 1000) 
  x <- as.data.frame(x) 
  for (i in 1:1000) { 
    x[[i]] <- length(list[[i]]) == ncol(comp) ##### NUMBER OF PROBES/MARKERS STUDIED HERE ###### 
  } 
  
  if (length(which(x[1,] %in% "TRUE")) > 0) { 
    comp <- comp[-which(x %in% "TRUE"),] 
  } else { 
    comp <- comp
  } 
  
  # Probes of small effect
  t <- vector() 
  list <- apply(comp, 1, function(x) which(x %in% 1)) 
  for (i in 1:1000) { 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true <- which(t %in% "TRUE")
  ind_false <- which(t %in% "FALSE")
  list_true <- list[ind_true]
  list_false <- list[ind_false] 
  n <- length(list_true) 
  m1_1_true <- matrix(0, ncol = 1, nrow = n)
  m1_1_true <- as.data.frame(m1_1_true) 
  m1_1_true$ind <- ind_true
  x <- vector()
  for (j in m1_1_true$ind) { 
    x[j] <- sum((betas[j, list[[j]]])^2) 
  } 
  m1 <- as.data.frame(x) 
  m1$x[is.na(m1$x)] <- 0 
  names(m1) <- "Variance_Small_Effects" 
  
  # Probes of medium effect
  t <- vector() 
  list <- apply(comp, 1, function(x) which(x %in% 2)) 
  for (i in 1:1000) { 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true <- which(t %in% "TRUE")
  ind_false <- which(t %in% "FALSE")
  list_true <- list[ind_true]
  list_false <- list[ind_false] 
  n <- length(list_true) 
  m2_true <- matrix(0, ncol = 1, nrow = n)
  m2_true <- as.data.frame(m2_true) 
  m2_true$ind <- ind_true
  x <- vector()
  for (j in m2_true$ind) { 
    x[j] <- sum((betas[j, list[[j]]])^2) 
  } 
  m2 <- as.data.frame(x) 
  m2$x[is.na(m2$x)] <- 0 
  names(m2) <- "Variance_Medium_Effects"
  
  # Probes of large effect
  t <- vector() 
  list <- apply(comp, 1, function(x) which(x %in% 3)) 
  for (i in 1:1000){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true <- which(t %in% "TRUE")
  ind_false <- which(t %in% "FALSE")
  list_true <- list[ind_true]
  list_false <- list[ind_false] 
  n <- length(list_true) 
  m3_true <- matrix(0, ncol = 1, nrow = n)
  m3_true <- as.data.frame(m3_true) 
  m3_true$ind <- ind_true
  x <- vector()
  for (j in m3_true$ind) { 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m3 <- as.data.frame(x) 
  m3$x[is.na(m3$x)] <- 0 
  names(m3) <- "Variance_Large_Effects"
  
  # Fuse
  m1$num <- row.names(m1) 
  m2$num <- row.names(m2) 
  m3$num <- row.names(m3) 
  all <- merge(m1, m2, by = "num", all = T) 
  var <- merge(all, m3, by = "num", all = T) 
  var[is.na(var)] <- 0 
  var$num <- NULL
  var$Total_Variance <- var[,1] + var[,2] + var[,3]
  var$Proportion_Small_Effects <- var[,1]/var[,4]
  var$Proportion_Medium_Effects <- var[,2]/var[,4]
  var$Proportion_Large_Effects <- var[,3]/var[,4]
  output$Proportion_Small_Effects <- mean(var$Proportion_Small_Effects) 
  output$Proportion_Medium_Effects <- mean(var$Proportion_Medium_Effects) 
  output$Proportion_Large_Effects <- mean(var$Proportion_Large_Effects) 
  
  # Export
  write.table(var, file = paste("summary/", gsub("_processed.csv*", "", file), "_varianceexplained_periteration.tsv", sep = ""), row.names = F, quote = F, sep = "\t") 
  write.table(output, file = paste("summary/", gsub("_processed.csv*", "", file), "_varianceexplained.tsv", sep = ""), row.names = F, quote = F, sep = "\t") 

}


# (3) Calculate Posterior Inclusion Probability (PIP) of CpGs as well as median/mean betas
##########################################################################################

cpgs <- names$Marker
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
anno <- anno[cpgs,c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]

for (i in loop) { 

  tmp <- fread(paste0("beta/", i))
  tmp <- as.data.frame(tmp)
  print(i)
  A <- gsub("_processed.csv*", "", i)
  median <- apply(tmp, 2, median)
  median_LLCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.025))
  median_LCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.05))
  median_HCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.95))
  median_HHCI <- apply(tmp, 2, function(x) quantile(x, probs =  0.975))
  names(median) <- cpgs
  names(median_LLCI) <- cpgs
  names(median_LCI) <- cpgs
  names(median_HCI) <- cpgs
  names(median_HHCI) <- cpgs
  median <- as.matrix(median)
  median_LCI <- as.matrix(median_LCI)
  median_LLCI <- as.matrix(median_LLCI)
  median_HCI <- as.matrix(median_HCI)
  median_HHCI <- as.matrix(median_HHCI)
  median_LCI <- cbind(median_LLCI, median_LCI)
  median_HCI <- cbind(median_HCI, median_HHCI)
  betas <- cbind(median,cbind(median_LCI,median_HCI))
  betas <- as.data.frame(betas)
  betas$CpG <- row.names(betas)
  betas <- betas[,c(6,1,2,3,4,5)]
  names(betas)[2:6] <- c("Median_Beta", "Beta_2.5", "Beta_5", "Beta_95", "Beta_97.5")
  comp <- fread(paste0("./comp/", i))
  comp <- as.data.frame(comp)
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- cpgs
  betas1 <- cbind(anno, betas[,c(2,3)], pip[,2])
  print(A)
  write.table(betas1, paste0("summary/", A, "_medianbeta_pip.tsv"), row.names = F, sep = "\t", quote = F)

} 

# Calculate Mean Betas
for (i in loop) {

  tmp <- fread(paste0("beta/", i))
  tmp <- as.data.frame(tmp)
  print(i)
  A <- gsub("_processed.csv*","",i)
  means <- apply(tmp, 2, mean)
  ses <- apply(tmp, 2, function(x){sd(x)/sqrt(length(x))})
  names(means) <- cpgs
  names(ses) <- cpgs
  means <- as.matrix(means)
  ses <- as.matrix(ses)  
  betas <- cbind(means,ses)
  betas <- as.data.frame(betas)
  betas$CpG <- row.names(betas)
  betas <- betas[,c(3,1,2)]
  names(betas)[2:3] <- c("Mean_Beta", "SE")
  comp <- fread(paste0("comp/", i))
  comp <- as.data.frame(comp)
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- cpgs
  betas1 <- cbind(anno, betas[,c(2,3)], pip[,2])
  write.table(betas1, paste0("summary/", A, "_meanbeta_pip.tsv"), row.names = F, quote = F, sep = "\t")

  # Sig CpGs
  for (j in c(0.95, 0.90, 0.85, 0.80)) {
    sig_cpgs <- betas1[betas1$PIP > j,]
    sig_cpgs <- sig_cpgs[order(-sig_cpgs$PIP),]
    write.table(sig_cpgs, paste0("summary/", A, "_meanbeta_pip", j,".tsv"), row.names = F, quote = F, sep = "\t")
  }

} 

