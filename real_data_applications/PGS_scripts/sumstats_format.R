library(data.table)
path0 <- getwd()

for (i in 1:5) {
for (k in 1:5) {
gwas <- fread(paste0(path0,"/GWAS_training_",i,"/remove_subset_",k,"/training_",i,"_rm_",k,"_chr1.HDL.glm.linear"),header = T)

for (chr in 2:22) {
  chr.sumstat <- fread(paste0(path0,"/GWAS_training_",i,"/remove_subset_",k,"/training_",i,"_rm_",k,"_chr",chr,".HDL.glm.linear"),header = T)
  gwas <- rbind(gwas,chr.sumstat)
  cat("chr",chr,"done \n")
}

gwas$P <- as.numeric(gwas$P)
gwas <- gwas[gwas$P != 0 & !is.na(gwas$P),]
gwas$P <- pmax(gwas$P,.Machine$double.xmin)
colnames(gwas)[1] <- "CHR"
colnames(gwas)[3] <- "SNP"
colnames(gwas)[4] <- "A2"
colnames(gwas)[5] <- "A1"

gwas <- gwas[,c(1:3,5,4,9,11:15)]

training_pheno = fread(paste0("~/Documents/PRS/lassosum/HDL/example_HDL_covar/training_test_HDL_covar/training_",i,"/training_pheno_remove_subset_",k,".txt"),header = T)
gwas$N <- nrow(training_pheno)

write.table(gwas, file = paste0(path0,"/GWAS_format/HDL_training_",i,"_rm_",k,".sumstats.txt"), sep = "\t",
            quote = F, row.names = F, col.names = T)
}
}
