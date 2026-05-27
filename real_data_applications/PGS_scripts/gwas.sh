#!/bin/bash

Dir=/home/r9user7/Documents/PRS/lassosum
workDir=/home/r9user7/Documents/PRS/lassosum/HDL/example_HDL_covar

for i in {1..5}; do
for k in {1..5}; do
for chr in {1..22}; do

/home/r9user7/plink/plink2 --bfile $Dir/UKBB_bfile/chr${chr} \
                           --pheno $workDir/training_test_HDL_covar/training_${i}/training_pheno_remove_subset_${k}.txt \
                           --pheno-name HDL \
                           --pheno-quantile-normalize \
                           --covar $workDir/training_test_HDL_covar/training/training_covar_${i}.txt \
                           --covar-name \
                             SEX,AGE,PC1-PC20 \
                           --covar-variance-standardize \
                           --glm hide-covar \
                           --threads 60 \
                           --out $workDir/GWAS_training_${i}/remove_subset_${k}/training_${i}_rm_${k}_chr${chr}

done
done
done