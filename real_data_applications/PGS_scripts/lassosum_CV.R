library(lassosum)
library(data.table)
library(parallel)

for (i in 1:5) {
  for (k in 1:5) {
    hdl.ss <- fread(paste0("GWAS_format/HDL_training_",i,"_rm_",k,".sumstats.txt"))
    ref.bfile <- "ref_panel"
    test.bfile <- paste0("subset_bfile/training_",i,"/subset_",k)
    cor <- p2cor(p = hdl.ss$P, n = hdl.ss$N[1], sign=hdl.ss$BETA)
    LDblocks <- "EUR.hg19"
    
    cl <- makeCluster(30) 
    out <- lassosum.pipeline(cor=cor, chr=hdl.ss$CHR, pos=hdl.ss$POS, 
                             A1=hdl.ss$A1, A2=hdl.ss$A2, 
                             ref.bfile=ref.bfile, test.bfile=test.bfile, 
                             LDblocks = LDblocks, cluster = cl)
    stopCluster(cl)
    
    test_pheno = fread(paste0("training_test_HDL_covar/training_",i,"/training_pheno_subset_",k,".txt"))
    v <- validate(out, pheno = test_pheno)
    out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)
    
    test.bim <- fread(paste0("subset_bfile/training_",i,"/subset_",k,".bim"),header = F)
    snp <- test.bim[out2$also.in.refpanel,]
    
    weights <- out2$sumstats
    weights$chr <- as.numeric(weights$chr)
    weights$beta <- out2$beta[[1]][,1]
    weights$snp <- snp$V2
    weights <- weights[,c(1,8,2:4,7)]
    colnames(weights) = c("CHR","SNP","POS","A1","A2","BETA")
    write.table(weights, paste0("weights/weights_training_",i,"_remove_subset_",k,".txt"), 
                sep = "\t", quote = F, row.names = F, col.names = T)
  }
}
