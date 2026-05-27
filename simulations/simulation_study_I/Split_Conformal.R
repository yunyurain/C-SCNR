library(data.table)

conf_level = 0.95
alpha = 1 - conf_level

for (h2 in c(0.5)) {
  for (causal in c(0.01,0.1,0.5)) {
      cat("heritability =",h2,", polygenicity =",causal,"...\n")
      for (rep in 1:10) {
        ## calibration data
        k = 5
        
        #### Residuals on Calibration data
        calib.pheno = fread(paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, "_cont_pheno.txt"),header=T)
        calib.covar = fread(paste0("UKBB_individuals/covar/covar_subset_",k,".txt"),header = T)
        calib.pgs = fread(paste0("pgs/Split_Conformal/training/","h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, ".profile"),header=T)
        calib.pgs = calib.pgs[match(calib.pheno$IID,calib.pgs$IID),]
        
        calib.data = calib.covar[,c("SEX","AGE","PC1")]
        calib.data$PGS = calib.pgs$SCORESUM
        calib.data$pheno = calib.pheno$PHENO
        
        model = lm(pheno ~ PGS, data = calib.data)
        calib.pheno$predicted_PHENO = model$fitted.values
        
        calib.pheno$residual = abs(calib.pheno$PHENO-calib.pheno$predicted_PHENO)
        ordered_R_i = sort(calib.pheno$residual)
        
        #### Test data
        test.pheno = fread(paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_cont_pheno_test.txt"),header=T)
        test.covar = fread("UKBB_individuals/covar/covar_test.txt",header = T)
        test.PGS = fread(paste0("pgs/Split_Conformal/test/","h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_", k, "_removed.profile"),header=T)
        test.PGS = test.PGS[match(test.pheno$IID,test.PGS$IID),]
        
        test.data = test.covar[,c("SEX","AGE","PC1")]
        test.data$PGS = test.PGS$SCORESUM
        test.data$pheno = test.pheno$PHENO
        
        test.model = lm(pheno ~ PGS, data = test.data)
        test.pheno$predicted_PHENO = test.model$fitted.values
        
        ## construct prediction intervals
        test.pheno$lower = test.pheno$predicted_PHENO - ordered_R_i[floor((1-alpha)*(nrow(calib.data)+1))]
        test.pheno$upper = test.pheno$predicted_PHENO + ordered_R_i[floor((1-alpha)*(nrow(calib.data)+1))]

        write.table(test.pheno, paste0("prediction_interval/Split_Conformal/","h2_", h2, "_poly_", causal, "_PI_test_rep_", rep, ".txt"), 
                    sep = "\t", row.names = F, quote = F, col.names = T)
        
        cat("replicate",rep,"done","\n")
      }
  }
}

