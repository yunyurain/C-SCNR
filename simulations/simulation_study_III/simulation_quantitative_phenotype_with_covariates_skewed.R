library(data.table)
path = "/home/r9user7/Documents/PRS/Simulation_Split_Conformal_CalPred_skewed"

for (rep in 1:20) {
  ## get dimensions
  geno_fam = fread("UKBB_individuals/individuals.chr1.fam", header = F)
  n_sample = dim(geno_fam)[1]
  
  #### simulate phenotype Y ####
  ## simulate Y hat
  y_hat = rnorm(n_sample, mean = 0, sd = 1)
  
  ## get covariates and quantile normalize
  sample_covar = fread("UKBB_individuals/covar/covar.txt", header = T)
  covar = sample_covar[,c("SEX","AGE","PC1")]
  scaled_covar = scale(covar)
  scaled_covar = as.data.frame(scaled_covar)
  colnames(scaled_covar) = c("SEX","AGE","PC1")
  
  ## simulate epsilon
  eps = rexp(n_sample, rate = 1)
  epsilon = rep(0, n_sample)
  for (k in 1:n_sample) {
    var = exp(log(7/3) + 0.2*scaled_covar$SEX[k] + 0.25*scaled_covar$AGE[k] + 0.15*scaled_covar$PC1[k])
    epsilon[k] = -sqrt(var) + sqrt(var)*eps[k]
  }
  
  ## calculate the final simulated Y
  y <- y_hat + epsilon
  
  ## output simulated phenotype
  pheno_output <- data.frame(geno_fam$V1, geno_fam$V2, y_hat, y)
  colnames(pheno_output) <- c("FID", "IID", "Y_hat", "PHENO")
  pheno_output_name <- paste0(path, "/simulated_data/rep_", rep, "_cont_pheno_all_samples.txt")
  write.table(pheno_output, pheno_output_name, row.names=F, col.names=T, quote=F, sep=" ")
}

      


