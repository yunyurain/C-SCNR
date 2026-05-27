#!/bin/bash

Dir=/home/r9user7/Documents/PRS/lassosum/HDL/example_HDL_covar

for i in {1..5}; do
for k in {1..5}; do

/home/r9user7/plink/plink --bfile $Dir/subset_bfile/training_${i}/subset_${k} \
                          --score $Dir/weights/weights_training_${i}_remove_subset_${k}.txt 2 4 6 sum \
                          --out $Dir/PGS_CV/training/PGS_training_${i}_subset_${k}

/home/r9user7/plink/plink --bfile $Dir/test_bfile/test_${i} \
                          --score $Dir/weights/weights_training_${i}_remove_subset_${k}.txt 2 4 6 sum \
                          --out $Dir/PGS_CV/test/PGS_test_${i}_remove_subset_${k}

done
done
