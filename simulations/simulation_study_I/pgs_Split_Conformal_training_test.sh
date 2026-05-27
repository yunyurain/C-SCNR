#!/bin/bash

Dir=/home/r9user7/Documents/PRS/Simulation_CV+_Split_Conformal

for causal in {0.01,0.1,0.5}; do
for rep in {1..10}; do

for k in {1..5}; do

/home/r9user7/plink/plink --bfile $Dir/UKBB_individuals/subset_bfile/subset_${k} \
                          --score $Dir/weights/h2_0.5_poly_${causal}_rep_${rep}_subset_5_removed.txt 2 4 6 sum \
                          --threads 20 \
                          --out $Dir/pgs/Split_Conformal/training/h2_0.5_poly_${causal}_rep_${rep}_subset_${k}

done

/home/r9user7/plink/plink --bfile $Dir/UKBB_individuals/test \
                          --score $Dir/weights/h2_0.5_poly_${causal}_rep_${rep}_subset_5_removed.txt 2 4 6 sum \
                          --threads 20 \
                          --out $Dir/pgs/Split_Conformal/test/h2_0.5_poly_${causal}_rep_${rep}_subset_5_removed

done
done