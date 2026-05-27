library(data.table)
library(dplyr)
path = "/home/r9user7/Documents/PRS/Simulation_clustering_effective"

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

  ## Scenario 1
  covar <- covar %>%
    mutate(
      AGE_Group = if_else(AGE < 55, 1, 2),
      Base_Group = paste0("S", SEX, "A", AGE_Group) 
    )
  
  ## Scenario 2
  ## mutate(AGE_Group = case_when(AGE < 50 ~ 1, AGE >= 50 & AGE <= 60 ~ 2,AGE > 60 ~ 3),
  ##        Base_Group = paste0("S", SEX, "A", AGE_Group))
  ## Scenario 3
  ## mutate(AGE_Group = if_else(AGE < 55, 1, 2), PC1_Group = if_else(PC1 < median(PC1), 1, 2),
  ##        Base_Group = paste0("S", SEX, "A", AGE_Group, "P", PC1_Group))
  
  scaled_covar = scale(covar[,c("SEX","AGE","PC1")])
  scaled_covar = as.data.frame(scaled_covar)
  colnames(scaled_covar) = c("SEX","AGE","PC1")
  scaled_covar$Base_Group = covar$Base_Group

  ## Scenario 1
  group_effect_map <- data.frame(
    Base_Group = c("S0A1","S0A2","S1A1","S1A2"),
    cluster_multiplier = c(0.5, 2.0, 0.5, 2.0)
  )
  
  ## Scenario 2
  ## Base_Group = c("S0A1","S0A2","S0A3","S1A1","S1A2","S1A3"), 
  ## cluster_multiplier = c(0.5, 1.0, 2.0, 0.5, 1.0, 2.0)
  ## Scenario 3
  ## Base_Group = c("S0A1P1", "S0A1P2", "S0A2P1", "S0A2P2","S1A1P1", "S1A1P2", "S1A2P1", "S1A2P2"),
  ## cluster_multiplier = c(0.5, 1.0, 1.0, 2.0,0.5, 1.0, 1.0, 2.0)
  
  scaled_covar <- scaled_covar %>%
    left_join(group_effect_map, by = "Base_Group")
  
  ## simulate epsilon
  epsilon = rep(0, n_sample)
  log_sigma_sq <- log(7/3) + 0.2*scaled_covar$SEX + 0.25*scaled_covar$AGE + 0.15*scaled_covar$PC1
  for (k in 1:n_sample) {
    epsilon[k] = rnorm(1, mean=0, sd=exp(0.5 * log_sigma_sq[k]) * sqrt(scaled_covar$cluster_multiplier[k]))
  }
  
  ## calculate the final simulated Y
  y <- y_hat + epsilon
  
  ## output simulated phenotype
  pheno_output <- data.frame(geno_fam$V1, geno_fam$V2, y_hat, y)
  colnames(pheno_output) <- c("FID", "IID", "Y_hat", "PHENO")
  pheno_output_name <- paste0(path, "/simulated_data/rep_", rep, "_cont_pheno_all_samples.txt")
  write.table(pheno_output, pheno_output_name, row.names=F, col.names=T, quote=F, sep=" ")
}

      


