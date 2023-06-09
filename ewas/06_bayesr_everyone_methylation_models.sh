#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Elena

export LD_LIBRARY_PATH=/opt/gcc/lib64


############################################################################
# Script to run BayesR+ on everyone
############################################################################

# Seed 1
cd /Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/methylation/chain1/
/Cluster_Filespace/Marioni_Group/BayesRRcmd/src/brr --data-file /Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/methylation/GS20k_allchrom_everyone_nogeno_somecovars.csv --pheno /Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/allindividuals_16717_logalcohol_residualized_everyone.csvphen --analysis-type preprocess --thread 24 --thread-spawned 24 --marker-cache --seed 1 
/Cluster_Filespace/Marioni_Group/BayesRRcmd/src/brr --data-file /Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/methylation/GS20k_allchrom_everyone_nogeno_somecovars.csv --pheno /Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/data/allindividuals_16717_logalcohol_residualized_everyone.csvphen --analysis-type ppbayes --chain-length 10000 --burn-in 100 --thin 1  --S "0.01,0.001,0.0001" --mcmc-samples allindividuals_16717_logalcohol_residualized_everyone.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1


############# Parse output

# Parse output
cd /Cluster_Filespace/Marioni_Group/Elena/alcohol_consumption/results/bayesr/usual_drinkers/output/methylation/

for i in ./chain1/*.csv
do
    
sigma1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'NR==1') 
sigma2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'END{print $NF}') 
A=$( echo $i | cut -d"/" -f3 )
B=$( echo $A | cut -d"." -f1 )
cat $i | cut -d ',' -f $sigma1-$sigma2 > ../sigma/${B}.csv

beta1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
beta2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3 )
B=$( echo $A | cut -d"." -f1 )
cat $i | cut -d ',' -f $beta1-$beta2 > ../beta/${B}.csv

comp1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
comp2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3 )
B=$( echo $A | cut -d"." -f1 )
cat $i | cut -d ',' -f $comp1-$comp2 > ../comp/${B}.csv

done 
