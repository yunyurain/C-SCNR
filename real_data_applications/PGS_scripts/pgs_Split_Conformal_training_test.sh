#!/bin/bash

Dir=/home/r9user7/Documents/PRS/lassosum/HDL/example_HDL_covar

for i in {1..5}; do

/home/r9user7/plink/plink --bfile $Dir/subset_bfile/training_${i}/subset_5 \
                          --score $Dir/weights/weights_training_${i}_remove_subset_5.txt 2 4 6 sum \
                          --out $Dir/PGS/training/PGS_training_${i}_subset_5
                           
/home/r9user7/plink/plink --bfile $Dir/test_bfile/test_${i} \
                          --score $Dir/weights/weights_training_${i}_remove_subset_5.txt 2 4 6 sum \
                          --out $Dir/PGS/test/PGS_test_${i}_remove_subset_5

done
