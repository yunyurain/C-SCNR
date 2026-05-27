#!/bin/bash

Dir=/home/r9user7/Documents/PRS/Simulation_Split_Conformal_CalPred

for h2 in {0.2,0.5,0.8}; do
for causal in {0.01,0.1,0.5}; do
for phi in {0.05,0.15}; do
for rep in {1..10}; do

/home/r9user7/plink/plink2 --bfile $Dir/UKBB_individuals/individuals.chr1 \
                           --pheno $Dir/simulated_data/h2_${h2}_poly_${causal}_phi_${phi}_rep_${rep}_subset_5_removed_cont_pheno.txt \
                           --pheno-name PHENO \
                           --pheno-quantile-normalize \
                           --covar $Dir/UKBB_individuals/covar/covar_remove_subset_5.txt \
                           --covar-name SEX,AGE,PC1 \
                           --covar-variance-standardize \
                           --glm hide-covar \
                           --threads 40 \
                           --out $Dir/sumstat/h2_${h2}_poly_${causal}_phi_${phi}_rep_${rep}_subset_5_removed

done
done
done
done