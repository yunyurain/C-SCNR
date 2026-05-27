library(data.table)
library(tidyverse)
library(NbClust)

path0 = "~/Documents/PRS/lassosum/HDL/example_HDL_covar"

predict <- function(mean_mat, sd_mat, mean_coef, sd_coef) {
  y_mean <- mean_mat %*% mean_coef
  y_var <- exp(sd_mat %*% sd_coef)
  return(data.frame(
    mean = y_mean,
    sd = sqrt(y_var)
  ))
}

## 2*6 Interaction implemented in the study
## Base groups can be defined by the user given the available covariates
get_base_groups <- function(df) {
  df %>%
    mutate(
      AGE_Group = cut(AGE, 
                      breaks = c(-Inf, 45, 50, 55, 60, 65, Inf),
                      labels = c(1:6),
                      right = FALSE),
      Base_Group = paste0("S", SEX, "A", AGE_Group)
    )
}

conf_level = 0.95
alpha = 1 - conf_level

for (i in 1:5) {
  ## Calibration data
  k = 5
  
  #### Training data
  training.pheno = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_pheno_remove_subset_",k,".txt"), header = T)
  
  training.covar = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",k,".txt"), header = T)
  training.covar = training.covar[FALSE,]
  for (j in setdiff(1:5,k)) {
    temp = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",j,".txt"), header = T)
    training.covar = rbind(training.covar,temp)
  };rm(temp)
  
  training.pgs = fread(paste0(path0,"/PGS/training/PGS_training_",i,"_subset_",k,".profile"), header = T)
  training.pgs = training.pgs[FALSE,]
  for (j in setdiff(1:5,k)) {
    temp = fread(paste0(path0,"/PGS/training/PGS_training_",i,"_subset_",j,".profile"), header = T)
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
  
  #### Studentized Residuals on Calibration data
  calib.pheno = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_pheno_subset_",k,".txt"), header = T)
  calib.covar = fread(paste0(path0,"/training_test_HDL_covar/training_",i,"/training_covar_subset_",k,".txt"), header = T)
  calib.pgs = fread(paste0(path0,"/PGS/training/PGS_training_",i,"_subset_",k,".profile"), header = T)
  calib.pgs = calib.pgs[match(calib.pheno$IID,calib.pgs$IID),]
  
  calib.data = calib.covar[,2:14]
  calib.data$PGS = calib.pgs$SCORESUM
  calib.data$HDL = calib.pheno$HDL
  calib.data$AGE_2 = (calib.data$AGE)^2
  calib.data$AGE_SEX = calib.data$AGE * calib.data$SEX
  calib.data$Intercept = 1
  
  mean_mat_cal = as.matrix(calib.data[,c("Intercept","PGS","AGE","AGE_2","AGE_SEX","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
  sd_mat_cal = as.matrix(calib.data[,c("Intercept","AGE","SEX","PC1","PC2")])
  pred.cal <- predict(mean_mat=mean_mat_cal, sd_mat=sd_mat_cal, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
  calib.data$R_i = abs(calib.data$HDL - pred.cal$mean) / pred.cal$sd
  
  for (rep in 1:20) {
  #### Clustering on Calibration data
  ## Base Groups  
  calib.data = get_base_groups(calib.data)
  
  ## Split into Set A and B
  idx_a = sample(1:nrow(calib.data), size = floor(nrow(calib.data)/2))
  cal_a <- calib.data[idx_a,]
  cal_b <- calib.data[-idx_a,]
  
  ## Feature Vectors for Base Groups
  feature_matrix_df <- cal_a %>%
    group_by(Base_Group) %>%
    summarise(
      q50 = quantile(R_i, 0.50),
      tail_ratio = quantile(R_i, conf_level) / q50
    ) %>%
    column_to_rownames("Base_Group") %>% 
    as.matrix()
  
  feature_matrix_scaled <- scale(feature_matrix_df)
  
  ## Clustering on Set A
  nb_result <- NbClust(feature_matrix_scaled, 
                       distance = "euclidean", 
                       min.nc = 1, max.nc = 4, 
                       method = "ward.D2", 
                       index = "silhouette")
  
  best_k <- nb_result$Best.nc[1]
  cat("Best Number of Clusters:", best_k, "\n")
  
  hc <- hclust(dist(feature_matrix_scaled), method = "ward.D2")
  cluster_assignments <- cutree(hc, k = best_k)
  
  group_mapping <- data.frame(
    Base_Group = names(cluster_assignments),
    Cluster_ID = as.factor(cluster_assignments)
  )
  
  ## Cluster-wise Quantile on Set B
  cal_b <- cal_b %>%
    left_join(group_mapping, by = "Base_Group")
  
  cluster_cutoffs <- cal_b %>%
    group_by(Cluster_ID) %>%
    summarise(
      q_k = quantile(R_i, probs = 1 - alpha, names = FALSE)
    )
  
  calibration_results <- group_mapping %>%
    left_join(cluster_cutoffs, by = "Cluster_ID")
  
  #### Test data
  test.pheno = fread(paste0(path0,"/training_test_HDL_covar/test/test_pheno_",i,".txt"), header = T)
  test.covar = fread(paste0(path0,"/training_test_HDL_covar/test/test_covar_",i,".txt"), header = T)
  test.PGS = fread(paste0(path0,"/PGS/test/PGS_test_",i,"_remove_subset_",k,".profile"), header = T)
  test.PGS = test.PGS[match(test.pheno$IID,test.PGS$IID),]
  
  test.data = test.covar[,2:14]
  test.data$PGS = test.PGS$SCORESUM
  test.data$HDL = test.pheno$HDL
  test.data$AGE_2 = (test.data$AGE)^2
  test.data$AGE_SEX = test.data$AGE * test.data$SEX
  test.data$Intercept = 1
  
  mean_mat_test = as.matrix(test.data[,c("Intercept","PGS","AGE","AGE_2","AGE_SEX","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
  sd_mat_test = as.matrix(test.data[,c("Intercept","AGE","SEX","PC1","PC2")])
  pred.test <- predict(mean_mat=mean_mat_test, sd_mat=sd_mat_test, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
  
  ## Cluster-wise Mapping
  test.data = get_base_groups(test.data)
  test.data <- test.data %>%
    left_join(calibration_results, by = "Base_Group")
  
  ## Construct Prediction Intervals
  test.pheno$point = pred.test$mean
  test.pheno$lower = pred.test$mean - pred.test$sd * test.data$q_k
  test.pheno$upper = pred.test$mean + pred.test$sd * test.data$q_k
  
  #### Output Test Results
  test.pheno$SEX = test.data$SEX
  test.pheno$AGE = test.data$AGE
  test.pheno$PC1 = test.data$PC1
  test.pheno$Base_Group = test.data$Base_Group
  test.pheno$Cluster_ID = test.data$Cluster_ID
  
  out.path = paste0(path0,"/prediction_interval/C-SCNR")
  write.table(calibration_results, paste0(out.path,"/C-SCNR_",i,"_rep_",rep,".txt"),
              sep = "\t", row.names = F, quote = F, col.names = T)
  write.table(test.pheno, paste0(out.path,"/PI_test_",i,"_rep_",rep,".txt"), 
              sep = "\t", row.names = F, quote = F, col.names = T)
  
  cat("rep =",rep,"done","\n")
  }
  
  cat("i =",i,"done","\n")
}


