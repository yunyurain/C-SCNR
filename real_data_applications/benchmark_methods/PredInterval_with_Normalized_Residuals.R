library(data.table)
path0 = "~/Documents/PRS/lassosum/HDL/example_HDL_covar"

predict <- function(mean_mat, sd_mat, mean_coef, sd_coef) {
  y_mean <- mean_mat %*% mean_coef
  y_var <- exp(sd_mat %*% sd_coef)
  return(data.frame(
    mean = y_mean,
    sd = sqrt(y_var)
  ))
}

conf_level = 0.95
alpha = 1 - conf_level

for (i in 1:5) {
  ## CV to obtain normalized residuals across subsets
  test_point = test_lower = test_upper = rep(list(NULL), 60000)
  for (k in 1:5) {
    #### Training data
    training.pheno = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_pheno_remove_subset_",k,".txt"), header = T)
    
    training.covar = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",k,".txt"), header = T)
    training.covar = training.covar[FALSE,]
    for (j in setdiff(1:5,k)) {
      temp = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",j,".txt"), header = T)
      training.covar = rbind(training.covar,temp)
    };rm(temp)
    
    training.pgs = fread(paste0(path0,"/PGS_CV/training/PGS_training_",i,"_subset_",k,".profile"), header = T)
    training.pgs = training.pgs[FALSE,]
    for (j in setdiff(1:5,k)) {
      temp = fread(paste0(path0,"/PGS_CV/training/PGS_training_",i,"_subset_",j,".profile"), header = T)
      training.pgs = rbind(training.pgs,temp)
    };rm(temp)
    training.pgs = training.pgs[match(training.pheno$IID,training.pgs$IID),]
    
    training.data = training.covar[,3:14]
    training.data$PGS = training.pgs$SCORESUM
    training.data$HDL = training.pheno$HDL
    training.data$AGE_2 = (training.data$AGE)^2
    training.data$AGE_SEX = training.data$AGE * training.data$SEX
    training.data$SEX = as.factor(training.data$SEX)
    
    cat("start fitting...")
    
    fit = lm(HDL ~ PGS+AGE+AGE_2+AGE_SEX+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = training.data)
    training.data$HDL_hat = fit$fitted.values
    training.data$epsilon = fit$residuals
    fit2 = glm(epsilon^2 ~ AGE+SEX+PC1+PC2, data = training.data, family = gaussian(link = "log"))
    
    cat("done.","\n")
    
    #### Normalized Residuals on Calibration data
    calib.pheno = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_pheno_subset_",k,".txt"), header = T)
    calib.covar = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",k,".txt"), header = T)
    calib.pgs = fread(paste0(path0,"/PGS_CV/training/PGS_training_",i,"_subset_",k,".profile"), header = T)
    calib.pgs = calib.pgs[match(calib.pheno$IID,calib.pgs$IID),]
    
    calib.data = calib.covar[,3:14]
    calib.data$PGS = calib.pgs$SCORESUM
    calib.data$HDL = calib.pheno$HDL
    calib.data$AGE_2 = (calib.data$AGE)^2
    calib.data$AGE_SEX = calib.data$AGE * calib.data$SEX
    calib.data$Intercept = 1
    
    mean_mat_cal = as.matrix(calib.data[,c("Intercept","PGS","AGE","AGE_2","AGE_SEX","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
    sd_mat_cal = as.matrix(calib.data[,c("Intercept","AGE","SEX","PC1","PC2")])
    
    pred.cal <- predict(mean_mat=mean_mat_cal, sd_mat=sd_mat_cal, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
    
    R_i = abs(calib.data$HDL - pred.cal$mean) / pred.cal$sd
    
    #### Test data
    test.pheno = fread(paste0(path0,"/training_test_HDL_covar/test/test_pheno_",i,".txt"), header = T)
    test.covar = fread(paste0(path0,"/training_test_HDL_covar/test/test_covar_",i,".txt"), header = T)
    test.PGS = fread(paste0(path0,"/PGS_CV/test/PGS_test_",i,"_remove_subset_",k,".profile"), header = T)
    test.PGS = test.PGS[match(test.pheno$IID,test.PGS$IID),]
    
    test.data = test.covar[,3:14]
    test.data$PGS = test.PGS$SCORESUM
    test.data$HDL = test.pheno$HDL
    test.data$AGE_2 = (test.data$AGE)^2
    test.data$AGE_SEX = test.data$AGE * test.data$SEX
    test.data$Intercept = 1
    
    mean_mat_test = as.matrix(test.data[,c("Intercept","PGS","AGE","AGE_2","AGE_SEX","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
    sd_mat_test = as.matrix(test.data[,c("Intercept","AGE","SEX","PC1","PC2")])
    
    pred.test <- predict(mean_mat=mean_mat_test, sd_mat=sd_mat_test, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
    
    for (j in 1:nrow(test.pheno)) {
      test_point[[j]] = c(test_point[[j]], pred.test$mean[j])
      test_lower[[j]] = c(test_lower[[j]], pred.test$mean[j] - pred.test$sd[j] * R_i)
      test_upper[[j]] = c(test_upper[[j]], pred.test$mean[j] + pred.test$sd[j] * R_i) 
    }
  }
  
  point = lower_bound = upper_bound = c()
  ## Calculate lower and upper bounds for each individual in test
  for (j in 1:nrow(test.pheno)) {
    ordered_lower = sort(test_lower[[j]])
    lb = ordered_lower[ceiling(alpha*(length(ordered_lower)+1))]
    
    ordered_upper = sort(test_upper[[j]])
    ub = ordered_upper[ceiling((1-alpha)*(length(ordered_upper)+1))]
    
    point <- c(point, mean(test_point[[j]]))
    lower_bound <- c(lower_bound, lb)
    upper_bound <- c(upper_bound, ub)
  }
  
  ## output prediction intervals
  test.pheno$point = point
  test.pheno$lower = lower_bound
  test.pheno$upper = upper_bound
  
  test.pheno$SEX = test.data$SEX
  test.pheno$AGE = test.data$AGE
  test.pheno$PC1 = test.data$PC1
  
  out.path = paste0(path0,"/prediction_interval/PredInterval_NR")
  write.table(test.pheno, paste0(out.path,"/PI_test_",i,".txt"), 
              sep = "\t", row.names = F, quote = F, col.names = T)
  
  cat("i =",i,"done","\n")
}

      
