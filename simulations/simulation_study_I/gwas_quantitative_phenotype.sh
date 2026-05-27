#!/bin/bash

Dir=/home/r9user7/Documents/PRS/Simulation_CV+_Split_Conformal

for causal in {0.1,0.5}; do
for rep in {1..10}; do
for k in {1..5}; do

/home/r9user7/plink/plink2 --bfile $Dir/UKBB_individuals/individuals.chr1 \
                           --pheno $Dir/simulated_data/h2_0.5_poly_${causal}_rep_${rep}_subset_${k}_removed_cont_pheno.txt \
                           --pheno-name PHENO \
                           --pheno-quantile-normalize \
                           --glm allow-no-covars \
                           --threads 20 \
                           --out $Dir/sumstat/h2_0.5_poly_${causal}_rep_${rep}_subset_${k}_removed

done
done
done