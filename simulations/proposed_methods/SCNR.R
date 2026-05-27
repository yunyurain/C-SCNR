library(data.table)

predict <- function(mean_mat, sd_mat, mean_coef, sd_coef) {
  y_mean <- mean_mat %*% mean_coef
  y_var <- exp(sd_mat %*% sd_coef)
  return(data.frame(
    mean = y_mean,
    sd = sqrt(y_var)
  ))
}

## Specify alpha
conf_level = 0.95
alpha = 1 - conf_level

for (h2 in c(0.2,0.5,0.8)) {
  for (causal in c(0.01,0.1,0.5)) {
    for (phi in c(0.05,0.15)) {
      cat("heritability =",h2,",polygenicity =",causal,",phi =",phi,"...\n")
      for (rep in 1:10) {
        ## Calibration data
        k = 5
        
        #### Training data
        training.pheno = fread(paste0("simulated_data/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, "_removed_cont_pheno.txt"),header=T)
        training.covar = fread(paste0("UKBB_individuals/covar/covar_remove_subset_",k,".txt"),header = T)
        training.pgs = fread(paste0("pgs/training/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, ".profile"),header=T)
        training.pgs = training.pgs[FALSE,]
        for (j in setdiff(1:5,k)) {
          temp = fread(paste0("pgs/training/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", j, ".profile"),header=T)
          training.pgs = rbind(training.pgs,temp)
        };rm(temp)
        training.pgs = training.pgs[match(training.pheno$IID,training.pgs$IID),]
        
        training.data = training.covar[,c("SEX","AGE","PC1")]
        training.data$PGS = training.pgs$SCORESUM
        training.data$pheno = training.pheno$PHENO
        
        cat("start fitting...")
        
        fit = lm(pheno ~ PGS+SEX+AGE+PC1, data = training.data)
        training.data$pheno_hat = fit$fitted.values
        training.data$epsilon = fit$residuals
        fit2 = glm(epsilon^2 ~ SEX+AGE+PC1, data = training.data, family = gaussian(link = "log"))
        
        cat("done.","\n")
        
        #### Normalized Residuals on Calibration data
        calib.pheno = fread(paste0("simulated_data/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, "_cont_pheno.txt"),header=T)
        calib.covar = fread(paste0("UKBB_individuals/covar/covar_subset_",k,".txt"),header = T)
        calib.pgs = fread(paste0("pgs/training/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, ".profile"),header=T)
        calib.pgs = calib.pgs[match(calib.pheno$IID,calib.pgs$IID),]
        
        calib.data = calib.covar[,c("SEX","AGE","PC1")]
        calib.data$PGS = calib.pgs$SCORESUM
        calib.data$pheno = calib.pheno$PHENO
        calib.data$Intercept = 1
        
        mean_mat_cal = as.matrix(calib.data[,c("Intercept","PGS","SEX","AGE","PC1")])
        sd_mat_cal = as.matrix(calib.data[,c("Intercept","SEX","AGE","PC1")])
        
        pred.cal <- predict(mean_mat=mean_mat_cal, sd_mat=sd_mat_cal, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))

        R_i = abs(calib.data$pheno - pred.cal$mean) / pred.cal$sd
        ordered_R_i = sort(R_i)
        
        #### Test data
        test.pheno = fread(paste0("simulated_data/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_cont_pheno_test.txt"),header=T)
        test.covar = fread("UKBB_individuals/covar/covar_test.txt",header = T)
        test.PGS = fread(paste0("pgs/test/h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, "_removed.profile"),header=T)
        test.PGS = test.PGS[match(test.pheno$IID,test.PGS$IID),]
        
        test.data = test.covar[,c("SEX","AGE","PC1")]
        test.data$PGS = test.PGS$SCORESUM
        test.data$pheno = test.pheno$PHENO
        test.data$Intercept = 1
        
        mean_mat_test = as.matrix(test.data[,c("Intercept","PGS","SEX","AGE","PC1")])
        sd_mat_test = as.matrix(test.data[,c("Intercept","SEX","AGE","PC1")])
        
        pred.test <- predict(mean_mat=mean_mat_test, sd_mat=sd_mat_test, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
        
        ## construct prediction intervals
        test.pheno$point = pred.test$mean
        test.pheno$lower = pred.test$mean - pred.test$sd * ordered_R_i[floor((1-alpha)*(nrow(calib.data)+1))]
        test.pheno$upper = pred.test$mean + pred.test$sd * ordered_R_i[floor((1-alpha)*(nrow(calib.data)+1))]
        
        #### Output prediction intervals
        test.pheno$SEX = test.data$SEX
        test.pheno$AGE = test.data$AGE
        test.pheno$PC1 = test.data$PC1
        
        out.path = paste0("~/Documents/PRS/Simulation_Split_Conformal_CalPred/prediction_interval/Split_Conformal_Normalized_Residual/h2_", h2, "_poly_", causal, "_phi_", phi)
        write.table(test.pheno, paste0(out.path,"_PI_test_rep_",rep,".txt"), 
                    sep = "\t", row.names = F, quote = F, col.names = T)
        
        cat("replicate",rep,"done","\n")
      }
    }
  }
}
