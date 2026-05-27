library(data.table)
path0 = "~/Documents/PRS/lassosum/HDL/example_HDL_covar"

for (i in 1:5) {
  for (k in 1:5) {
    #### Predicted Y in training subsets 
    subset.pheno.name = paste0(path0,"/training_test_HDL_covar/training_",i,"/training_pheno_subset_",k,".txt")
    subset.pheno = fread(subset.pheno.name, header = T)
    subset.covar = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",k,".txt"), header = T)
    subset.pgs.name = paste0(path0,"/PGS_CV/training/PGS_training_",i,"_subset_",k,".profile")
    subset.pgs = fread(subset.pgs.name, header = T)
    subset.pgs = subset.pgs[match(subset.pheno$IID,subset.pgs$IID),]
    
    subset.data = subset.covar[,3:14]
    subset.data$PGS = subset.pgs$SCORESUM
    subset.data$HDL = subset.pheno$HDL
    subset.data$AGE_2 = (subset.data$AGE)^2
    subset.data$AGE_SEX = subset.data$AGE * subset.data$SEX
    subset.data$SEX = as.factor(subset.data$SEX)
    
    subset.model = lm(HDL ~ ., data = subset.data)
    subset.pheno$predicted_HDL = as.numeric(subset.model$fitted.values)
    
    ## calculate Residual_CV
    subset.pheno$residual_CV = abs(subset.pheno$HDL - subset.pheno$predicted_HDL)
    
    write.table(subset.pheno, paste0(path0,"/predicted/PredInterval/training_",i,"_subset_",k,".txt"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    
    #### Predicted Y in test 
    test.pheno.name = paste0(path0,"/training_test_HDL_covar/test/test_pheno_",i,".txt")
    test.pheno = fread(test.pheno.name, header = T)
    test.covar = fread(paste0(path0,"/training_test_HDL_covar/test/test_covar_",i,".txt"), header = T)
    test.pgs.name = paste0(path0,"/PGS_CV/test/PGS_test_",i,"_remove_subset_",k,".profile")
    test.pgs = fread(test.pgs.name, header = T)
    test.pgs = test.pgs[match(test.pheno$IID,test.pgs$IID),]
    
    test.data = test.covar[,3:14]
    test.data$PGS = test.pgs$SCORESUM
    test.data$HDL = test.pheno$HDL
    test.data$AGE_2 = (test.data$AGE)^2
    test.data$AGE_SEX = test.data$AGE * test.data$SEX
    test.data$SEX = as.factor(test.data$SEX)
    
    test.model = lm(HDL ~ ., data = test.data)
    test.pheno$predicted_HDL = as.numeric(test.model$fitted.values)
    
    write.table(test.pheno, paste0(path0,"/predicted/PredInterval/test_",i,"_subset_",k,"_removed.txt"),
                sep = "\t", quote = F, row.names = F, col.names = T)
  }
  cat("i =",i,"done","\n")
}

