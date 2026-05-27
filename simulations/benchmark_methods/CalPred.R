library(data.table)

train <- function(mean_mat, sd_mat, y, tol = 1e-6, maxit = 100) {
  fit <- statmod::remlscore(y = y, X = mean_mat, Z = sd_mat, tol = tol, maxit = maxit)
  mean_coef <- as.vector(fit$beta)
  sd_coef <- as.vector(fit$gamma)
  mean_se <- fit$se.beta
  sd_se <- fit$se.gam
  
  names(mean_coef) <- colnames(mean_mat)
  names(sd_coef) <- colnames(sd_mat)
  names(mean_se) <- colnames(mean_mat)
  names(sd_se) <- colnames(sd_mat)
  
  return(list(
    mean_coef = mean_coef,
    mean_se = mean_se,
    sd_coef = sd_coef,
    sd_se = sd_se
  ))
}

predict <- function(mean_mat, sd_mat, mean_coef, sd_coef) {
  if (!all.equal(colnames(mean_mat), names(mean_coef))) {
    stop("colnames(mean_mat) != names(mean_coef)")
  }
  if (!all.equal(colnames(sd_mat), names(sd_coef))) {
    stop("colnames(sd_mat) != names(sd_coef)")
  }
  y_mean <- mean_mat %*% mean_coef
  y_var <- exp(sd_mat %*% sd_coef)
  return(data.frame(
    mean = y_mean,
    sd = sqrt(y_var)
  ))
}

normalize_reference <- function(x) {
  q <- qnorm((rank(x) - 0.5) / length(x))
  # remove duplicated values in x
  idx <- !duplicated(x)
  x <- x[idx]
  q <- q[idx]
  q2x <- approxfun(q, x, rule = 2)
  x2q <- approxfun(x, q, rule = 2)
  return(list(q2x = q2x, x2q = x2q))
}

conf_level = 0.95
alpha = 1 - conf_level

for (h2 in c(0.2,0.5,0.8)) {
  for (causal in c(0.01,0.1,0.5)) {
    for (phi in c(0.05,0.15)) {
      cat("heritability =",h2,",polygenicity =",causal,",phi =",phi,"...\n")
      for (rep in 1:10) {
        ## calibration data
        k = 5
        
        #### Calibration data
        calib.pheno = fread(paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, "_cont_pheno.txt"),header=T)
        calib.covar = fread(paste0("UKBB_individuals/covar/covar_subset_",k,".txt"),header = T)
        calib.pgs = fread(paste0("pgs/training/","h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, ".profile"),header=T)
        calib.pgs = calib.pgs[match(calib.pheno$IID,calib.pgs$IID),]
        
        calib.data = calib.covar[,c("SEX","AGE","PC1")]
        calib.data$PGS = calib.pgs$SCORESUM
        calib.data$pheno = calib.pheno$PHENO
        calib.data$Intercept = 1
        
        #### Test data
        test.pheno = fread(paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_cont_pheno_test.txt"),header=T)
        test.covar = fread("UKBB_individuals/covar/covar_test.txt",header = T)
        test.PGS = fread(paste0("pgs/test/","h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_subset_", k, "_removed.profile"),header=T)
        test.PGS = test.PGS[match(test.pheno$IID,test.PGS$IID),]
        
        test.data = test.covar[,c("SEX","AGE","PC1")]
        test.data$PGS = test.PGS$SCORESUM
        test.data$pheno = test.pheno$PHENO
        test.data$Intercept = 1
        
        ## build reference to convert quantile <-> value using the calibration dataset
        qref <- normalize_reference(calib.data$pheno)
        calib.data$pheno_q = qref$x2q(calib.data$pheno)

        #### Use calibration data to train the model
        mean_mat_cal = as.matrix(calib.data[,c("Intercept","PGS","SEX","AGE","PC1")])
        sd_mat_cal = as.matrix(calib.data[,c("Intercept","SEX","AGE","PC1")])
        pheno_cal = calib.data$pheno_q
        
        cat("start fitting...")
        
        fit <- train(mean_mat=mean_mat_cal, sd_mat=sd_mat_cal, y=pheno_cal)
        
        cat("done.","\n")
     
        ## construct prediction interval on test data
        mean_mat_test = as.matrix(test.data[,c("Intercept","PGS","SEX","AGE","PC1")])
        sd_mat_test = as.matrix(test.data[,c("Intercept","SEX","AGE","PC1")])
        
        pred.test <- predict(mean_mat=mean_mat_test, sd_mat=sd_mat_test, mean_coef=fit$mean_coef, sd_coef=fit$sd_coef)
        
        ## CalPred
        test.pheno$calpred_point = qref$q2x(pred.test$mean)
        test.pheno$calpred_lower = qref$q2x(pred.test$mean - pred.test$sd * qnorm(1-alpha/2))
        test.pheno$calpred_upper = qref$q2x(pred.test$mean + pred.test$sd * qnorm(1-alpha/2))
        
        #### Output prediction intervals
        test.pheno$SEX = test.data$SEX
        test.pheno$AGE = test.data$AGE
        test.pheno$PC1 = test.data$PC1
        
        out.path = paste0("~/Documents/PRS/Simulation_Split_Conformal_CalPred/prediction_interval/CalPred/h2_", h2, "_poly_", causal, "_phi_", phi)
        write.table(test.pheno, paste0(out.path,"_PI_test_rep_",rep,".txt"), 
                    sep = "\t", row.names = F, quote = F, col.names = T)
        
        cat("replicate",rep,"done","\n")
      }
    }
  }
}
