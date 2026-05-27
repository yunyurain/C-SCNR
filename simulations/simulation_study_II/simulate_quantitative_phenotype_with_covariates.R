library(data.table)
path = "/home/r9user7/Documents/PRS/Simulation_Split_Conformal_CalPred"

## set parameters
for (h2 in c(0.2,0.5,0.8)) {
  for (causal in c(0.01,0.1,0.5)) {
    for (phi in c(0.05,0.15)) {
      for (rep in 1:10) {
        ## get dimensions
        geno_bim = fread("UKBB_individuals/individuals.chr1.bim", header = F)
        geno_fam = fread("UKBB_individuals/individuals.chr1.fam", header = F)
        n_snp = dim(geno_bim)[1]
        n_sample = dim(geno_fam)[1]
        
        #### simulate SNP effect sizes ####
        total_index = 1:n_snp
        causal_index = sample(total_index, size = round(causal*n_snp), replace = F)
        beta = rep(0,n_snp)
        for (i in 1:length(causal_index)) {
          beta[causal_index[i]] <- rnorm(1, mean=0, sd=sqrt(h2/round(causal*n_snp)))
        }
        
        maf_data = fread("UKBB_individuals/individuals.chr1.frq", header = T)
        
        ## calculate non-scaled beta
        beta_noscl <- beta / sqrt(2*maf_data$ALT_FREQS*(1-maf_data$ALT_FREQS))
        
        ## output the simulated beta using PLINK score format
        beta_output <- data.frame(maf_data$ID, maf_data$ALT, beta_noscl)
        beta_name <- paste0(path,"/simulated_data/", "h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_beta.txt")
        write.table(beta_output, beta_name, row.names=F, col.names=F, quote=F, sep=" ")
        
        #### simulate phenotype Y ####
        ## calculate X*beta by PLINK score
        score_name <- paste0(path,"/simulated_data/", "h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_X_times_beta")
        
        system(paste0("~/plink/plink --bfile ", "~/Documents/PRS/Simulation_Split_Conformal_CalPred/UKBB_individuals/individuals.chr1",  
                      " --score ",  beta_name, " 1 2 3 sum",
                      " --out ",    score_name))
        
        x_beta_prod <- fread(paste0(score_name, ".profile"))
        
        ## get covariates and quantile normalize
        sample_covar = fread("UKBB_individuals/covar/covar.txt", header = T)
        covar = sample_covar[,c("SEX","AGE","PC1")]
        scaled_covar = scale(covar)
        scaled_covar = as.data.frame(scaled_covar)
        colnames(scaled_covar) = c("SEX","AGE","PC1")
        
        ## simulate covariate effects
        gamma = rnorm(3, mean = 0, sd = sqrt(phi/3))
        
        ## simulate epsilon
        beta_sigma_0 = log(1-h2-phi)-0.2*mean(scaled_covar$SEX)-0.25*mean(scaled_covar$AGE)-0.15*mean(scaled_covar$PC1)
        epsilon <- rep(0, n_sample)
        for (k in 1:n_sample){
          epsilon[k] <- rnorm(1, mean=0, sd=sqrt(exp(beta_sigma_0+0.2*scaled_covar$SEX[k]+0.25*scaled_covar$AGE[k]+0.15*scaled_covar$PC1[k])))
        }
        
        ## calculate the final simulated Y
        unscaled_pheno <- x_beta_prod$SCORESUM + as.matrix(scaled_covar) %*% gamma + epsilon
        #scaled_pheno <- scale(unscaled_pheno) #scale 
        
        ## output simulated phenotype
        pheno_output <- data.frame(x_beta_prod$FID, x_beta_prod$IID, unscaled_pheno)
        colnames(pheno_output) <- c("FID", "IID", "PHENO")
        pheno_output_name <- paste0(path, "/simulated_data/", "h2_", h2, "_poly_", causal, "_phi_", phi, "_rep_", rep, "_cont_pheno_all_samples.txt")
        write.table(pheno_output, pheno_output_name, row.names=F, col.names=T, quote=F, sep=" ")
      }
    }
  }
}


