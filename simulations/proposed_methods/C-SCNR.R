library(data.table)
library(tidyverse)
library(NbClust)

predict <- function(mean_mat, sd_mat, mean_coef, sd_coef) {
  y_mean <- mean_mat %*% mean_coef
  y_var <- exp(sd_mat %*% sd_coef)
  return(data.frame(
    mean = y_mean,
    sd = sqrt(y_var)
  ))
}

## Simulation Study IV (Scenario 1)
## Base groups can be flexibly defined by the user given the available covariates
get_base_groups <- function(df) {
  df %>%
    mutate(
      AGE_Group = if_else(AGE < 55, 1, 2),
      Base_Group = paste0("S", SEX, "A", AGE_Group) 
    )
}

conf_level = 0.95
alpha = 1 - conf_level

L = 20

for (rep in 1:20) {
  ## Calibration data
  k = 5
  
  #### Training data
  training.pheno = fread(paste0("simulated_data/rep_", rep, "_subset_", k, "_removed_cont_pheno.txt"),header=T)
  training.covar = fread(paste0("UKBB_individuals/covar/covar_remove_subset_",k,".txt"),header = T)
  
  training.data = training.covar[,c("SEX","AGE","PC1")]
  training.data$Y = training.pheno$Y_hat
  training.data$pheno = training.pheno$PHENO
  
  cat("start fitting...")
  
  fit = lm(pheno ~ Y, data = training.data)
  training.data$pheno_hat = fit$fitted.values
  training.data$epsilon = fit$residuals
  fit2 = glm(epsilon^2 ~ SEX+AGE+PC1, data = training.data, family = gaussian(link = "log"))
  
  cat("done.","\n")
  
  #### Studentized Residuals on Calibration data
  calib.pheno = fread(paste0("simulated_data/rep_", rep, "_subset_", k, "_cont_pheno.txt"),header=T)
  calib.covar = fread(paste0("UKBB_individuals/covar/covar_subset_",k,".txt"),header = T)

  calib.data = calib.covar[,c("SEX","AGE","PC1")]
  calib.data$Y = calib.pheno$Y_hat
  calib.data$pheno = calib.pheno$PHENO
  calib.data$Intercept = 1
  
  mean_mat_cal = as.matrix(calib.data[,c("Intercept","Y")])
  sd_mat_cal = as.matrix(calib.data[,c("Intercept","SEX","AGE","PC1")])
  
  pred.cal <- predict(mean_mat=mean_mat_cal, sd_mat=sd_mat_cal, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
  
  calib.data$R_i = abs(calib.data$pheno - pred.cal$mean) / pred.cal$sd
  
  for (iter in 1:L) {
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
                         min.nc = 1, max.nc = 3, 
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
    test.pheno = fread(paste0("simulated_data/rep_", rep, "_cont_pheno_test.txt"),header=T)
    test.covar = fread("UKBB_individuals/covar/covar_test.txt",header = T)
    
    test.data = test.covar[,c("SEX","AGE","PC1")]
    test.data$Y = test.pheno$Y_hat
    test.data$pheno = test.pheno$PHENO
    test.data$Intercept = 1
    
    mean_mat_test = as.matrix(test.data[,c("Intercept","Y")])
    sd_mat_test = as.matrix(test.data[,c("Intercept","SEX","AGE","PC1")])
    
    pred.test <- predict(mean_mat=mean_mat_test, sd_mat=sd_mat_test, mean_coef=as.numeric(fit$coefficients), sd_coef=as.numeric(fit2$coefficients))
    
    ## Cluster-wise Mapping
    test.data = get_base_groups(test.data)
    test.data <- test.data %>%
      left_join(calibration_results, by = "Base_Group")
    
    ## Construct prediction intervals
    test.pheno$point = pred.test$mean
    test.pheno$lower = pred.test$mean - pred.test$sd * test.data$q_k
    test.pheno$upper = pred.test$mean + pred.test$sd * test.data$q_k
    
    #### Output results
    test.pheno$SEX = test.data$SEX
    test.pheno$AGE = test.data$AGE
    test.pheno$PC1 = test.data$PC1
    test.pheno$Base_Group = test.data$Base_Group
    test.pheno$Cluster_ID = test.data$Cluster_ID
    
    out.path = paste0("~/Documents/PRS/Simulation_clustering_effective/prediction_interval/C-SCNR")
    write.table(test.pheno, paste0(out.path,"/PI_test_rep_",rep,"_iter_",iter,".txt"), 
                sep = "\t", row.names = F, quote = F, col.names = T)
    
    cat("iter =",iter,"done","\n")
  }
  
  cat("replicate",rep,"done","\n")
}
      
