library(lassosum)
library(data.table)
library(parallel)

for (h2 in c(0.5)) {
for (causal in c(0.01,0.1,0.5)) {
for (rep in 1:10) {
  for (k in 1:5) {
    ss <- fread(paste0("sumstat/h2_",h2,"_poly_",causal,"_rep_",rep,"_subset_",k,"_removed.PHENO.glm.linear"))
    ss$P <- as.numeric(ss$P)
    ss <- ss[ss$P != 0 & !is.na(ss$P),]
    ss$P <- pmax(ss$P,.Machine$double.xmin)
    
    ref.bfile <- "ref_panel/ref_panel"
    test.bfile <- paste0("UKBB_individuals/subset_bfile/subset_",k)
    cor <- p2cor(p = ss$P, n = 20000, sign=ss$BETA)
    LDblocks <- "EUR.hg19"
    
    cl <- makeCluster(20) 
    out <- lassosum.pipeline(cor=cor, chr=ss$`#CHROM`, pos=ss$POS, 
                             A1=ss$ALT, A2=ss$REF, 
                             ref.bfile=ref.bfile, test.bfile=test.bfile, 
                             LDblocks = LDblocks, cluster = cl)
    stopCluster(cl)
    
    test_pheno = fread(paste0("simulated_data/", "h2_", h2, "_poly_", causal, "_rep_", rep, "_subset_",k,"_cont_pheno.txt"))
    v <- validate(out, pheno = test_pheno)
    out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)
    
    test.bim <- fread(paste0("UKBB_individuals/subset_bfile/subset_",k,".bim"),header = F)
    snp <- test.bim[out2$also.in.refpanel,]
    
    weights <- out2$sumstats
    weights$chr <- as.numeric(weights$chr)
    weights$beta <- out2$beta[[1]][,1]
    weights$snp <- snp$V2
    weights <- weights[,c(1,8,2:4,7)]
    colnames(weights) = c("CHR","SNP","POS","A1","A2","BETA")
    
    write.table(weights, paste0("weights/","h2_",h2,"_poly_",causal,"_rep_",rep,"_subset_",k,"_removed.txt"), 
                sep = "\t", quote = F, row.names = F, col.names = F)
  }
}
}
}

