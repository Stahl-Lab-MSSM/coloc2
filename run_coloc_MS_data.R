args <- commandArgs(trailingOnly=TRUE)
chromosome = args[1]
start = args[2]
end = args[3]
gwas_file = args[4]
eqtl_file = args[5]
outfolder = args[6]

# bed_input_file <- "/u/project/pasaniuc/pasaniucdata/DATA/PICKRELL_LOCI/EUR/fourier_ls-all.bed"
# bed = read.table(bed_input_file, header=T, sep="\t", stringsAsFactors=F)
# for (i in 1:nrow(bed)) {
# gwas_file <- "/u/project/pasaniuc/pasaniucdata/DATA/All_Summary_Statistics/2_Final/ALZHEIMERS_2013.sumstats"
# eqtl_file <- "/u/project/eeskin/pasaniuc/clagiamb/trans/chr1_10583_1892607/temp//chr1_10583_1892607_all_pval.tab"
# outfolder= "/u/project/eeskin/pasaniuc/clagiamb/trans_coloc/"

source("/u/project/eeskin/pasaniuc/clagiamb/scripts/coloc/coloc_functions.R")
library(data.table)

chromosome = as.numeric(gsub("chr", "", chromosome))
start = as.numeric(start)
end = as.numeric(end)

biom.df <- fread(gwas_file, header=T)
  biom.df <- data.frame(biom.df)
eqtl.df <- fread(eqtl_file, header=T)
  eqtl.df <- data.frame(eqtl.df)

#names(biom.df) = c("CHR", "POS", "SNP", "N", "A1", "A2", "Z", "BETA", "SE", "N_CONTROLS")
names(biom.df)[which(names(biom.df)=="SNPID")] = "SNP"
#biom.df$Ncases = biom.df$N-biom.df$N_CONTROLS
#biom.df$PVAL = 2*pnorm(-abs(biom.df$Z))

names(eqtl.df)[which(names(eqtl.df)=="SNPID")] = "SNP"
names(eqtl.df)[which(names(eqtl.df)=="REF")] = "A1"
names(eqtl.df)[which(names(eqtl.df)=="ALT")] = "A2"
names(eqtl.df)[which(names(eqtl.df)=="beta")] = "BETA"
names(eqtl.df)[which(names(eqtl.df)=="se.beta")] = "SE"


prefix = paste("merged", chromosome, start, end, sep="_")
merged.data = prepare.eqtl.biom(eqtl.df, biom.df, chromosome, start, end,  useBETA=TRUE, outfolder,prefix=prefix, bychrpos=TRUE, allele_merge=T, anyOverlap=TRUE)
# If this takes too long, submit one script per gene by specifying character vector for list.probes
prefix = paste("coloc", chromosome, start, end, sep="_")
coloc = coloc.eqtl.biom(merged.data, list.probes=NULL, useBETA=TRUE, outfolder, prefix= prefix, min_snps=50, p1=1e-4, p2=1e-4, p12=1e-6, moloc=FALSE)
