#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("missMethyl")

### PIP > 0.95
###################################################################################

# Import data
every <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_16717_logalcohol_residualized_everyone_meanbeta_pip.tsv", sep = "\t", header = T) 
usual <- read.table("/Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/summary/allindividuals_10506_logalcohol_residualized_usualdrinkers_onelesssample_meanbeta_pip.tsv", sep = "\t", header = T) 
every_s <- every[every$PIP > 0.95,]
usual_s <- usual[usual$PIP > 0.95,]

# CpGs
cpgs <- union(every_s$Name, usual_s$Name) # 8 total CpGs
all <- every$Name # 752722

# Annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Enrichment
go <- gometh(sig.cpg = cpgs, all.cpg = all, collection = "GO", prior.prob = TRUE, anno = ann)
kegg <- gometh(sig.cpg = cpgs, all.cpg = all, collection = "KEGG", prior.prob = TRUE, anno = ann)

# Total number of categories significant at 5% FDR
table(go$FDR<0.05) # None
table(kegg$FDR<0.05) # None

# Table of top results
topGSA(go)
topGSA(kegg)


### PIP > 0.8
###################################################################################

every_s <- every[every$PIP > 0.8,]
usual_s <- usual[usual$PIP > 0.8,]

# CpGs
cpgs <- union(every_s$Name, usual_s$Name) #21 total CpGs

# Enrichment
go <- gometh(sig.cpg = cpgs, all.cpg = all, collection = "GO", prior.prob = TRUE, anno = ann)
kegg <- gometh(sig.cpg = cpgs, all.cpg = all, collection = "KEGG", prior.prob = TRUE, anno = ann)

# Total number of categories significant at 5% FDR
table(go$FDR<0.05) # None
table(kegg$FDR<0.05) # None

# Table of top results
topGSA(go)
topGSA(kegg)

