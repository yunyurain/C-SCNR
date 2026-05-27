library(data.table)

for (h2 in c(0.5)) {
  for (causal in c(0.01,0.1,0.5)) {
      for (rep in 1:10) {
        for (k in 1:5) {
          #### Predicted Y in training subsets 
          subset.pheno.name = paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, "_cont_pheno.txt")
          subset.pheno = fread(subset.pheno.name, header = T)
          subset.covar = fread(paste0("UKBB_individuals/covar/covar_subset_",k,".txt"), header = T)
          subset.pgs.name = paste0("pgs/CV+/training/","h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, ".profile")
          subset.pgs = fread(subset.pgs.name, header = T)
          subset.pgs = subset.pgs[match(subset.pheno$IID,subset.pgs$IID),]
          
          subset.data = subset.covar[,c("SEX","AGE","PC1")]
          subset.data$SEX = as.factor(subset.data$SEX)
          subset.data$PGS = subset.pgs$SCORESUM
          subset.data$pheno = subset.pheno$PHENO
          
          model = lm(pheno ~ PGS, data=subset.data)
          subset.pheno$predicted_PHENO = as.numeric(model$fitted.values)
          
          ## calculate Residual_CV
          subset.pheno$residual_CV = abs(subset.pheno$PHENO-subset.pheno$predicted_PHENO)
          
          write.table(subset.pheno, paste0("predicted/PredInterval/","h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, ".txt"),
                      sep = "\t", quote = F, row.names = F, col.names = T)
          
          #### Predicted Y in test 
          test.pheno.name = paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_cont_pheno_test.txt")
          test.pheno = fread(test.pheno.name, header = T)
          test.covar = fread(paste0("UKBB_individuals/covar/covar_test.txt"), header = T)
          test.pgs.name = paste0("pgs/CV+/test/","h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, "_removed.profile")
          test.pgs = fread(test.pgs.name, header = T)
          test.pgs = test.pgs[match(test.pheno$IID,test.pgs$IID),]
          
          test.data = test.covar[,c("SEX","AGE","PC1")]
          test.data$SEX = as.factor(test.data$SEX)
          test.data$PGS = test.pgs$SCORESUM
          test.data$pheno = test.pheno$PHENO
          
          test.model = lm(pheno ~ PGS, data=test.data)
          test.pheno$predicted_PHENO = as.numeric(test.model$fitted.values)
          
          write.table(test.pheno, paste0("predicted/PredInterval/","h2_", h2, "_poly_", causal, "_test_rep_", rep, "_subset_", k, "_removed.txt"),
                      sep = "\t", quote = F, row.names = F, col.names = T)
        }
        cat("replicate",rep,"done","\n")
      }
  }
}
