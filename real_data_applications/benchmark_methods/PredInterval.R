library(data.table)
path0 = "~/Documents/PRS/lassosum/HDL/example_HDL_covar"

## Specifying alpha
conf_level = 0.95
alpha = 1 - conf_level

for (i in 1:5) {
  subset = test = list()
  ## Load Data
  for (k in 1:5) {
    subset.pheno = fread(paste0(path0,"/predicted/PredInterval/training_",i,"_subset_",k,".txt"), header = T)
    subset[[k]] = subset.pheno
    
    test.pheno = fread(paste0(path0,"/predicted/PredInterval/test_",i,"_subset_",k,"_removed.txt"), header = T)
    test[[k]] = test.pheno
  }
  
  lower_bound = upper_bound = c()
  point = c()
  ## Calculate lower and upper bounds for each individual in test
  for (j in 1:nrow(test[[1]])) {
    mu_minus_R_i = mu_plus_R_i = c()
    mu = c()
    for (k in 1:5) {
      mu_minus_R_i = c(mu_minus_R_i, test[[k]]$predicted_HDL[j] - subset[[k]]$residual_CV)
      mu_plus_R_i = c(mu_plus_R_i, test[[k]]$predicted_HDL[j] + subset[[k]]$residual_CV)
      mu = c(mu,test[[k]]$predicted_HDL[j])
    }
    # lower bound
    ordered_mu_minus_R_i <- mu_minus_R_i[order(mu_minus_R_i)]
    lb <- ordered_mu_minus_R_i[ceiling(alpha*(length(ordered_mu_minus_R_i)+1))]
    
    # upper bound
    ordered_mu_plus_R_i <- mu_plus_R_i[order(mu_plus_R_i)]
    ub <- ordered_mu_plus_R_i[ceiling((1-alpha)*(length(ordered_mu_plus_R_i)+1))]
    
    # output
    point <- append(point, mean(mu))
    lower_bound <- append(lower_bound, lb)
    upper_bound <- append(upper_bound, ub)
  }
  
  ## Output confidence intervals for all individuals in test
  output_data <- data.frame(IID = test[[1]]$IID, HDL = test[[1]]$HDL, point, lower_bound, upper_bound)
  write.table(output_data, paste0(path0,"/prediction_interval/PredInterval/PI_test_",i,".txt"),
              row.names=F, col.names=T, quote=F, sep="\t")
  
  cat("i =",i,"done","\n")
}
      

