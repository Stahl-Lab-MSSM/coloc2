# coloc2 


file with all functions here: functions_coloc_likelihood_summary_integrated.R


commands to run coloc2 for one tissue/trait combination

format the GWAS summary statistics file:

biom.df = formatColoc(fname = [GWAS_sumstats file], type=[case/control or quantitative trait; "cc" or "quant"], N=[GWAS sample size], Ncases=[number GWAS cases if case/control], info_filter=[info filter], maf_filter=[maf filter], fread=T, eqtl=FALSE)

format the eQTL summary statistics file:

eqtl.df = formatColoc(fname = [eQTL_sumstats file], type="quant", N=[eQTL sample size], Ncases=NA, info_filter=[info filter], maf_filter=[maf filter], fread=T, eqtl=TRUE)

run co-localization analysis:

res = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=TRUE, outfolder=[outfolder], prefix=[prefix], plot=FALSE, save.coloc.output=FALSE, match_snpid=T)


example commands:

biom.df = formatColoc(fname = "CAD_sumstats.txt", type="cc", N=100000, Ncases=30000, info_filter=0.6, maf_filter=0.05, fread=T, eqtl=FALSE)

eqtl.df = formatColoc(fname = "Adipose.txt", type="quant", N=400, Ncases=NA, info_filter=0.6, maf_filter=0.05, fread=T, eqtl=TRUE)

res = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=TRUE, outfolder="results/", prefix="CAD_Adipose", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
