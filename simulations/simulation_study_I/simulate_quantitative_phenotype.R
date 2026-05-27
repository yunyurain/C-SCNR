library(data.table)
path = "/home/r9user7/Documents/PRS/Simulation_CV+_Split_Conformal"

## set parameters
for (h2 in c(0.5)) {
  for (causal in c(0.01,0.1,0.5)) {
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
        beta_name <- paste0(path,"/simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_beta.txt")
        write.table(beta_output, beta_name, row.names=F, col.names=F, quote=F, sep=" ")
        
        #### simulate phenotype Y ####
        ## calculate X*beta by PLINK score
        score_name <- paste0(path,"/simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_X_times_beta")
        
        system(paste0("~/plink/plink --bfile ", "~/Documents/PRS/Simulation_CV+_Split_Conformal/UKBB_individuals/individuals.chr1",  
                      " --score ",  beta_name, " 1 2 3 sum",
                      " --out ",    score_name))
        
        x_beta_prod <- fread(paste0(score_name, ".profile"))
        
        ## simulate epsilon
        epsilon <- rep(0, n_sample)
        for (k in 1:n_sample){
          epsilon[k] <- rnorm(1, mean=0, sd=sqrt(1-h2))
        }

        ## calculate the final simulated Y
        unscaled_pheno <- x_beta_prod$SCORESUM + epsilon
        scaled_pheno <- scale(unscaled_pheno) #scale 
        
        ## output simulated phenotype
        pheno_output <- data.frame(x_beta_prod$FID, x_beta_prod$IID, scaled_pheno)
        colnames(pheno_output) <- c("FID", "IID", "PHENO")
        pheno_output_name <- paste0(path, "/simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_cont_pheno_all_samples.txt")
        write.table(pheno_output, pheno_output_name, row.names=F, col.names=T, quote=F, sep=" ")
      }
  }
}


