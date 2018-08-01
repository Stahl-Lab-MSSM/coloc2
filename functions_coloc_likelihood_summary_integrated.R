#### original code by claudia giambartolomei
#### modifications by james boocock, eli stahl, and amanda dobbyn

library(foreach)
library(doParallel)
library(GenomicRanges)
library(rtracklayer)
library(data.table)

Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}


get_n_effective  <-  function(f,N,s){
  # estimate var_y
  data_set$neff = (var_y - data_set$var *data_set$b^2)/(data_set$var *data_set$se^2) + 1

}

logsum <- function(x) {
  my.max <- max(x)                              ## get the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

logdiff <- function(x,y) {
  my.max <- max(x,y)                              ## get the maximum value in log form
  if (x>y) {
    my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  } 
  if (x<y) {
    my.res <- my.max + log(exp(y-my.max) - exp(x - my.max ))
  }
  return(my.res)
}


process.dataset <- function(d, suffix, ave=TRUE, estimate_sdy=TRUE,correlation=0) {
  message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop('The variable type must be set, otherwise the Bayes factors cannot be computed')

  if("beta" %in% nd && "varbeta" %in% nd && ("MAF" %in% nd || "sdY" %in% nd)) {
    if(length(d$beta) != length(d$varbeta))
      stop("Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("Length of snp names and beta vectors must match")
    if(estimate_sdy){
        if(d$type == 'quant' & !('sdY' %in% nd)){
           if("N" %in% nd){ 
                d$sdY <- sdY.est(d$varbeta, d$MAF, d$N,d$beta)
                print(d$sdY)
           }else{
                d$sdY  <- 1
           }
        }
    }
    print("SdY")
    print(d$sdY)
    d$sdY = 1
    print(d$sdY)
    if(correlation != 0){
        df = data.frame(Z,V,sdy) 
        df$Z=d$beta/sqrt(d$varbeta) 
        df$V =d$varbeta
        df$sdY=d$sdY
        if(!is.null(suffix))
            colnames(ret) <- paste(colnames(ret), suffix, sep=".")
        df$snp <- as.character(d$snp)
        return(df)
    }
    # if there are no negative values, assume OR
	if(d$type == "cc"){
		if (length(d$beta[d$beta<0])==0) d$beta = log(d$beta)
	} 
    
    if (d$type == "quant") {
        if(!ave) {
            df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                                      V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
        }else if(ave){
            df <- approx.bf.estimates.ave.pvalue(z=d$beta/sqrt(d$varbeta),
                                                 n=d$N,maf=d$MAF,type=d$type, suffix=suffix)
        }
    }else if(d$type =="cc"){
        if (!ave) {
            df <- approx.bf.estimates(z=(d$beta)/sqrt(d$varbeta),
                                      V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
        }
        if (ave) {
            df <- approx.bf.estimates.ave.pvalue(z=(d$beta)/sqrt(d$varbeta),
                                                 n=d$N, maf=d$MAF,s=d$s, type=d$type, suffix=suffix) 
        }
    }
    df$snp <- as.character(d$snp)
    return(df)
  }

  if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (length(d$pvalues) != length(d$MAF))
      stop('Length of the P-value vectors and MAF vector must match')
    if(d$type=="cc" & !("s" %in% nd))
      stop("Must specify s if type=='cc' and you want to use approximate Bayes Factors")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues,
                     MAF = d$MAF,
                     snp=as.character(d$snp))    
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
    df <- subset(df, df$MAF>0 & df$pvalues>0) # all p values and MAF > 0
    if(ave){
        abf <- approx.bf.estimates.ave.pvalue(pnorm(df$pvalues/2,lower.tail=F), maf=df$MAF,type=d$type,s=d$s,N=d$N,suffix=suffix)
    }else{
        abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    }
    df <- cbind(df, abf)
    return(df)  
  }

  stop("Must give, as a minimum, either (beta, varbeta, type) or (pvalues, MAF, N, type)")
}


coloc.abf <- function(dataset1, dataset2, MAF=NULL,
                      p1=1e-4, p2=1e-4, p12=1e-5, correlation=0) {
  if(!is.list(dataset1) || !is.list(dataset2))
    stop("dataset1 and dataset2 must be lists.")
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
    
  
  df1 <- process.dataset(d=dataset1, suffix="df1", ave=FALSE,correlation=correlation)
  df2 <- process.dataset(d=dataset2, suffix="df2", ave=TRUE,correlation=correlation)
  merged.df <- merge(df1,df2)
   

   if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  # if there are no columns with lABF computed from different sds, internal sum is just simple sum of lABF:
  if(correlation != 0){
    merged.df$approx.bf.estimates.ave(d$Z.df1,d$V.df1,d$Z.df2,d$V.df2,d$sdY.df1,d$sdY.df2)   
  }else{
   if (length(grep("lABF_sd1.df1", names(merged.df), value=T))==0) { 
   	merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  } else {
     merged.df$lABF_ave.df1 = apply(cbind(merged.df$lABF_sd1.df1, merged.df$lABF_sd2.df1, merged.df$lABF_sd3.df1),1, function(x) logsum(sum(x)) - log(3))
     merged.df$lABF_ave.df2 = apply(cbind(merged.df$lABF_sd1.df2, merged.df$lABF_sd2.df2, merged.df$lABF_sd3.df2),1, function(x) logsum(sum(x)) - log(3))

     merged.df$internal.sum.lABF <- apply(cbind(merged.df$lABF_sd1.df1 + merged.df$lABF_sd1.df2, merged.df$lABF_sd2.df1 + merged.df$lABF_sd2.df2, merged.df$lABF_sd3.df1 + merged.df$lABF_sd3.df2), 1, function(x) logsum(x) - log(3))
  }
  }

  
  my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
  
############################## 

  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
  
  output<-list(summary=results, results=merged.df)
  return(output)
}

approx.bf.p <- function(p,f,type, N, s, suffix=NULL) {
  if(type=="quant") {
    sd.prior <- 0.44 #0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.44 #0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}


lABF.fn <- function (z, V, sd.prior=0.15) {
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  return(lABF)
}

approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  lABF <- lABF.fn(z, V, sd.prior)
  ret <- data.frame(V, z, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

var_mle_from_z = function(z,n,maf){
    var_mle= 1/(2*maf*(1-maf) * ( n + z^2))
    return(var_mle)
}

var_mle_from_z_cc = function(z,n,maf,s){
    var_mle= (s*(1-s))/(2*maf*(1-maf) * ( n + z^2))
    return(var_mle)
}

var_mle_from_z_cc2 = function(z,n,maf,s){
    var_mle= 1/((2*maf*(1-maf) * ( n + z^2))*(s*(1-s)))
    return(var_mle)
}


approx.bf.estimates.ave.pvalue  <- function(z, n,type, maf, suffix=NULL, sdY=1, s=NA){
    z[is.nan(z)] = 0
    if(type == "cc"){
        V = var_mle_from_z_cc2(z,n,maf,s)
        return(approx.bf.estimates.ave(z,V,type,suffix,s))
    }else{
        V = var_mle_from_z(z,n,maf)
        return(approx.bf.estimates.ave(z,V,type,suffix))
    }
} 


approx.bf.estimates.ave <- function (z, V, type, suffix=NULL, sdY=1) {
  listVec <- list(lABF_sd1 = lABF.fn(z, V, sd.prior=sqrt(0.01)*sdY), lABF_sd2 = lABF.fn(z, V, sd.prior=sqrt(0.1)*sdY), lABF_sd3 = lABF.fn(z, V, sd.prior=sqrt(0.5)*sdY))
  m <- do.call(cbind, listVec)
  lABF <- apply(m, 1, function(x) logsum(x) -log(3))
  ret <- data.frame(V, z, m, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

lABF.fn.corr  <-  function(z1,V1,z2,sd.prior){
    r = 1 
}

lABF.fn.corr.coloc <-  function(z1,V1,z2,V2,sd.prior1,sd.prior2){
    r = 1
}



approx.bf.estimates.ave.corr <- function (z1, V1,z2,V2, type, suffix=NULL, sdY1=1,sdY2=1,correlation=0) {
  listVec <- list(lABF_sd1 = lABF.fn.corr(z1, V1,z2, sd.prior1=sqrt(0.01*sdY1),correlation=correlation), lABF_sd2 = lABF.fn.corr(z1,V1,z2,sd.prior1=sqrt(0.1*sdY1),correlation=correlation),lABF_sd3 = lABF.fn.corr(z1, V1,z2, sd.prior1=sqrt(0.5*sdY1),correlation=correlation))
  m <- do.call(cbind, listVec)
  lABF1<- apply(m, 1, function(x) logsum(x) -log(3))
  
  listVec <- list(lABF_sd1 = lABF.fn.corr(z2, V2,z1, sd.prior1=sqrt(0.01*sdY2),correlation=correlation), lABF_sd2 = lABF.fn.corr(z2,V2,z1,sd.prior1=sqrt(0.1*sdY2),correlation=correlation),lABF_sd3 = lABF.fn.corr(z2, V2,z1, sd.prior1=sqrt(0.5*sdY2),correlation=correlation))
  m <- do.call(cbind, listVec)
  lABF2 = apply(m, 1, function(x) logsum(x) -log(3))
  
  
  listVec <- list(lABF_sd1 = lABF.fn.corr.coloc(z1, V1,z2,V2, sd.prior1=sqrt(0.01*sdY1),sd.prior2=sqrt(0.01*sdY2),correlation=correlation), lABF_sd2 = lABF.fn.corr(z1,V1,z2,V2,sd.prior1=sqrt(0.1*sdY1),sd.prior2=sqrt(0.1*sdY2),correlation=correlation),lABF_sd3 = lABF.fn.corr(z1, V1,z2,V2,sd.prior1=sqrt(0.5*sdY1),sd.prior2=sqrt(0.5*sdY2),correlation=correlation))

  m <- do.call(cbind, listVec)
  lABF3 = apply(m, 1, function(x) logsum(x) -log(3))

  ret <- data.frame(V, z, lABF1,lABF2,lABF3)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}



complement_snp = function(x){
    as = x =="A"
    ts = x == "T"
    gs = x == "G"
    cs = x == "C"
    x[as] = "T"
    x[ts] = "A"
    x[gs] = "C"
    x[cs] = "G"
    return(x)
}

merge_results  <- function(a, b){
    if(is.null(a) & is.null(b)){
        return(NULL)
    }else if(is.null(a)){
        return(b)
    }else if(is.null(b)){
        return(a)
    }else{
        return(rbind(a,b))
    }
}

sdY.est <- function(vbeta, maf, n, beta) {
  vars = 2 *maf * ( 1- maf) * n * vbeta * (n -1 ) + 2 *maf * ( 1- maf) * n * beta^2 
  return(sqrt(median(vars/(n-1))))
}


# Remove duplicated IDs (duplicated snps and indels, keep only snps if have allele info)
remove_dupl = function(data) {
    n_occur <- data.frame(table(data$SNP))
    dupl = data[data$SNP %in% n_occur$Var1[n_occur$Freq > 1],]
         if (nrow(dupl)>0) {
          if (all(c("A1","A2") %in% names(data))) {
             dupl <- transform(dupl, n=nchar(as.character(dupl$A1)) + nchar(as.character(dupl$A2)))
             dupl=dupl[order(dupl$n, decreasing=T),]
          } else {
             dupl=dupl[order(dupl$MAF, decreasing=T),]
          }
          toremove = rownames(dupl[ !duplicated(dupl$SNP), ])
          removed_list <- data.frame(Marker_removed = dupl$SNP[!duplicated(dupl$SNP)], reason = "Duplicated SNPs")
          data = data[!(rownames(data) %in% toremove),]
          message("Removed ", length(toremove), " duplicated SNP names")
          }  else {
          removed_list <- data.frame(Marker_removed = NA, reason = "Duplicated SNPs")
          }
    return(list(data, removed_list))
}

combine.abf.locus <- function(l0, l1, l2, l3, l4, a0, a1, a2, a3, a4) {

  lH0.abf  <- log(a0) + l0
  lH1.abf  <- log(a1) + l1
  lH2.abf  <- log(a2) + l2
  lH3.abf  <- log(a3) + l3
  lH4.abf  <- log(a4) + l4

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  #print(signif(pp.abf,3))
  #print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}

combine.abf <- function(l1, l2, p1, p2, p12) {
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  #print(signif(pp.abf,3))
  #print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}


fn.pw.gwas = function(p, data) {
  a0 = p[1]
  a1 = p[2]
  a2 = p[3]
  a3 = p[4]
  a4 = p[5]
  #print(nrow(data))
  suma = sum(exp(c(a0,a1,a2,a3,a4)))
  #print(paste("Alphas:" , exp(a0)/suma,exp(a1)/suma,exp(a2)/suma,exp(a3)/suma,exp(a4)/suma), sep=" ")
  lkl.frame.temp <- as.matrix(data)
  lkl.frame.temp[,1] <- log(exp(a0)/suma)
  lkl.frame.temp[,2] <- lkl.frame.temp[,2] + log(exp(a1)/suma)
  lkl.frame.temp[,3] <- lkl.frame.temp[,3] + log(exp(a2)/suma)
  lkl.frame.temp[,4] <- lkl.frame.temp[,4] + log(exp(a3)/suma)
  lkl.frame.temp[,5] <- lkl.frame.temp[,5] + log(exp(a4)/suma)
  #print(lkl.frame.temp[1,])
  #print(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
  #print(log(sum(exp(lkl.frame.temp[1,]))))
  sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
  #print(sumlkl)
  return(sumlkl)
}


coloc.eqtl.biom <- function(eqtl.df, biom.df, p12=1e-6, useBETA=TRUE, plot=FALSE, outfolder, prefix= "pref", save.coloc.output=FALSE, match_snpid=TRUE,cores=20,bootstrap=F,no_bootstraps=1000, min_snps=50, bed_input_file=NULL){
  if (class(eqtl.df$ProbeID)!="character") stop("When reading the data frame, make sure class of ProbeID in eQTL data is a character")

  source("/sc/orga/projects/psychgen/resources/COLOC2/COLOC_scripts/scripts/optim_function.R")
  print(head(eqtl.df))
  print(head(biom.df))
  eqtl.df.rs = eqtl.df[grep("rs",eqtl.df$SNPID),]
  print(head(eqtl.df.rs)) 

# Estimate trait variance. 

if (!file.exists(outfolder)) dir.create(outfolder)
if (plot) {
   plot.fld = paste(outfolder, "plot/", sep="")
   pval.fld= paste(outfolder, "pval/", sep="")
   dir.create(file.path(plot.fld), showWarnings = FALSE)
   dir.create(file.path(pval.fld), showWarnings = FALSE)

   # For plotting with locuszoom
   refFlat_path = "/hpc/users/giambc02/scripts/locuszoom/refFlat.RData"
   source("/hpc/users/giambc02/scripts/locuszoom/call_locuszoom3_temp2.R")
   load(refFlat_path)
   refFlatRaw <- refFlatRaw.VP
}

  # Set Variables 
  maf_filter = 0.001 #MAF filter applied to datasets
  rsq_filter = 0.6 #Imputation quality filter applied to datasets

###
if (unique(biom.df$type) == "cc") cc=TRUE else cc=FALSE

maf.eqtl = ifelse("MAF" %in% names(eqtl.df), TRUE, FALSE)
maf.biom = ifelse("MAF" %in% names(biom.df), TRUE, FALSE)
if (!maf.eqtl & !maf.biom) message("There is no MAF information in either datasets, looking for frequency column") 


if (useBETA) {
   cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE", "ProbeID","N","A1","A2") # We need N only if we estimate sdYest
   cols.biom = c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE","N","type","A1","A2") # We need N only if we estimate sdYest
   }
if (!useBETA) {
   cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N","A1","A2")
   cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N","A1","A2")
   }

if (cc) cols.biom = c(cols.biom, "Ncases")
no_allele_merge =F
if (!all(cols.eqtl %in% names(eqtl.df))){
    if(!all(cols.eqtl[-which(cols.eqtl %in% c("A1","A2"))] %in% names(eqtl.df))){
        print(cols.eqtl[-which(cols.eqtl %in% c("A1","A2"))])
        stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
    }else{
        no_allele_merge = T
        cols.eqtl = cols.eqtl[-which(cols.eqtl %in% c("A1","A2"))]

        message("Warning: allele columns missing, will sometimes merge SNPs with indels")
    }
}


if (!all(  cols.biom %in% names(biom.df))){ 
    if(!all( cols.biom[-which(cols.biom %in% c("A1","A2"))] %in% names(biom.df))){
        print(all(c("A1","A2") %in% names(biom.df)))
        stop("These columns are missing from the biomarker data: ", cols.biom[!cols.biom %in% names(biom.df)])
    }else{
        no_allele_merge = T
        cols.biom = cols.biom[-which(cols.biom %in% c("A1","A2"))]
        message("Warning: allele columns missing, will sometimes merge SNPs with indels")
    }
}


#####################
# Filter by imputation quality if column exists
info.columns <- grep( names(biom.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    biom.df = subset(biom.df, biom.df[,info.columns] > rsq_filter)
    }
info.columns <- grep( names(eqtl.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    eqtl.df = subset(eqtl.df, eqtl.df[,info.columns] > rsq_filter)
    }

# use only one of the MAFs from the two datasets
# First check if there is a MAF in eQTL data and use this, if not take the one in biom data
# Filter by MAF


if (maf.eqtl) {
   cols.eqtl = c(cols.eqtl, "MAF")
    maf.eqtl = TRUE
   eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
}else if("F" %in% names(eqtl.df)){
    eqtl.df$MAF = ifelse(eqtl.df$F<0.5, eqtl.df$F, 1-eqtl.df$F)
    eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
    maf.eqtl = TRUE
    cols.eqtl = c(cols.eqtl, "F")
}else if (maf.biom) {
   cols.biom = c(cols.biom, "MAF")
   biom.df = subset(biom.df, biom.df$MAF > maf_filter)
   maf.biom= TRUE
}else if("F" %in% names(biom.df)){
   biom.df$MAF = ifelse(biom.df$F<0.5, biom.df$F, 1-biom.df$F)
   biom.df = subset(biom.df, biom.df$MAF > maf_filter)
   maf.biom= TRUE
   cols.biom = c(cols.biom, "F")
}

if (!maf.eqtl & !maf.biom) stop("There is no MAF information in either dataset")


################## SNPID MATCHING
# The reasoning here is that the user can give either only "SNPID", or "SNPID" and "input_names" (or we can find it as rsid or chr:pos and add it in the data as input_names)
# if there is no "input_names" column and the SNPID is not already chr:pos format, then the software will find the chr:pos and see if it matches better with the eQTL data than the SNPID provided


  hasChr=ifelse(any(grep("chr",eqtl.df$CHR[1:2]))>0, TRUE,FALSE)
    if (hasChr) (eqtl.df$CHR=gsub("chr", "", eqtl.df$CHR))
  hasChr=ifelse(any(grep("chr",biom.df$CHR[1:2]))>0, TRUE,FALSE)
    if (hasChr) (biom.df$CHR=gsub("chr", "", biom.df$CHR))


  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biom.df$SNPID))!=nrow(biom.df)) addChrposBiom = TRUE
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", eqtl.df$SNPID))!=nrow(eqtl.df)) addChrposEQTL = TRUE

  if (addChrposBiom) {
    biom.df$chrpos = paste(biom.df$CHR, biom.df$POS, sep=":")
  }

  if (addChrposEQTL) {
    eqtl.df$chrpos = paste(eqtl.df$CHR, eqtl.df$POS, sep=":")
  }

if (!match_snpid) {
# Find the combinations of SNPID that matches the most SNPs between the two datasets
biomSNPID = unique(biom.df$SNPID)
eqtlSNPID = unique(eqtl.df$SNPID)
match_snpid = max(length(biomSNPID[biomSNPID %in% eqtlSNPID]), length(eqtlSNPID[eqtlSNPID %in% biomSNPID]))
match_chrpos_snpid = 0
match_chrpos = 0

if (addChrposBiom) {
  biomchrpos = unique(biom.df$chrpos)
  match_chrpos_snpid = max(length(biomchrpos[biomchrpos %in% eqtlSNPID]), length(biomchrpos[biomchrpos %in% eqtlSNPID]))
}

if (addChrposEQTL) {
   eqtlchrpos = unique(eqtl.df$chrpos)
   if (!addChrposBiom) {
   match_chrpos_snpid = max(length(eqtlchrpos[eqtlchrpos %in% biomSNPID]), length(eqtlchrpos[eqtlchrpos %in% biomSNPID]))
   }
}

if (addChrposBiom & addChrposEQTL) match_chrpos = max(length(biomchrpos[biomchrpos %in% eqtlchrpos]), length(biomchrpos[biomchrpos %in% eqtlchrpos]))


find_best_column = which.max(c(match_snpid, match_chrpos_snpid, match_chrpos))

if (find_best_column==1) message("Best combination is SNPID: do not change column names")
if (find_best_column==2) {
    message("Best combination is SNPID in one dataset and input_name in the other: change one of the column names from input_name to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="chrpos"] <- "SNPID"
}
if (find_best_column==3) {
    message("Best combination is chrpos in both the datasets: change both of the column names from chrpos to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="chrpos"] <- "SNPID"
    names(eqtl.df)[names(eqtl.df)=="SNPID"] <- "SNPID2"
    names(eqtl.df)[names(eqtl.df)=="chrpos"] <- "SNPID"
}
} 
# if !match_snpid
#####################

  biom.df = biom.df[,cols.biom]
  if("F" %in% names(biom.df)){
    colnames(biom.df)[which("F" == names(biom.df))] = "MAF.biom"
  }   
  if("MAF" %in% names(biom.df)){
    colnames(biom.df)[which("MAF" == names(biom.df))] = "MAF.biom"
  }
  eqtl.df = eqtl.df[,cols.eqtl]
  if("F" %in% names(eqtl.df)){
    colnames(eqtl.df)[which("F" ==  names(eqtl.df))] = "MAF.eqtl"
  }   
  if("MAF" %in% names(eqtl.df)){
    colnames(eqtl.df)[which("MAF" == names(eqtl.df))] = "MAF.eqtl"
  }
  # Remove missing data
  eqtl.df = eqtl.df[complete.cases(eqtl.df),]
  biom.df = biom.df[complete.cases(biom.df),]

  res.all <- data.frame()

###################
# Find ProbeIDs overlapping with biom.df
if(!is.null(bed_input_file)){
    message("Reading LD independent bed file")
    bed = import.bed(bed_input_file)
    bed = bed[bed$ProbeID %in% unique(eqtl.df$ProbeID),]
}else{
   library(data.table)
   DT <- as.data.table(eqtl.df)
   expr_table <- DT[, list(CHR= unique(CHR), START = min(POS), STOP = max(POS), minP = min(PVAL)), by = ProbeID]
   expr_table <- data.frame(expr_table)
   message("There are ", nrow(expr_table), " ProbeIDs in the eQTL data")

   decreaseGap = FALSE
   if (decreaseGap) {
      targetGap=100000
      currentGap = (expr_table$STOP - expr_table$START)/2
      if (targetGap < currentGap[1]) {
        expr_table$START = expr_table$START + currentGap
        expr_table$STOP = expr_table$STOP - currentGap
      }
   }
   bed = expr_table
}

if (!all(c("ProbeID", "CHR", "START", "STOP") %in% names(bed))) stop("Bed file is missing info")
bed$ProbeID = as.character(bed$ProbeID)
bed$CHR=as.numeric(as.character(bed$CHR))
bed$START=as.numeric(bed$START)
bed$STOP=as.numeric(bed$STOP)
message("Looping through ", nrow(bed), " genes from the eQTL data")

##################################################### now start the loop
# Go over all regions that overlap between eQTL table and input.data

overlap.df<-merge(biom.df,eqtl.df,by="SNPID",suffixes=c(".biom",".eqtl"))

biom.df$sdY.biom = sdY.est((overlap.df$SE.biom)^2,overlap.df$MAF.eqtl,overlap.df$N.biom,overlap.df$BETA.biom) 
message("Running in parallel")
registerDoParallel(cores=cores)
list.probes = bed$ProbeID
print(class(eqtl.df))
eqtl.dfByProbe = split(seq(nrow(eqtl.df)), eqtl.df$ProbeID)
print(head(eqtl.df))
print(length(eqtl.dfByProbe[[1]]))
print("eqtl.df rows")
#print(head(eqtl.df$ProbeID))
if(!is.null(bed_input_file)){
    message("Reading LD independent bed file")
    bed = import.bed(bed_input_file)
}else{
    bed = NULL
}

duplicated_snp_list = data.frame()
res.all = data.frame()
for(i in 1:length(list.probes)){
       ProbeID = as.character(list.probes[i]) # important for probe names that are numbers
       region.eqtl <- eqtl.df[eqtl.dfByProbe[[ProbeID]],]
       pos.start <- min(region.eqtl$POS)
       pos.end   <- max(region.eqtl$POS)
       chrom = region.eqtl$CHR[1]
       matches <- which(chrom == biom.df$CHR & biom.df$POS >= pos.start & biom.df$POS <= pos.end )
       region.biom <- biom.df[matches, ]
       duplicated_snp_list = rbind(duplicated_snp_list, data.frame(ProbeID = ProbeID, data="biom",remove_dupl(region.biom)[[2]]))
       region.biom = remove_dupl(region.biom)[[1]]
       duplicated_snp_list = rbind(duplicated_snp_list, data.frame(ProbeID = ProbeID, data="eqtl",remove_dupl(region.eqtl)[[2]]))
       region.eqtl = remove_dupl(region.eqtl)[[1]]

       # Loop over each biomarker 
        
      if (cc) {
          type= "cc"
          region.biom$s1 = region.biom$Ncases/region.biom$N
         }
      if (!cc & !useBETA) {
          type = "quant"
          region.biom$s1=rep(0.5, length(region.biom$N)) ## This will be ignored since the type is "quant"
         }
         merged.data <- merge(region.biom, region.eqtl, by = "SNPID",  suffixes=c(".biom", ".eqtl"))
         if(!maf.eqtl){
            merged.data$MAF.eqtl = merged.data$MAF.biom
         }
         # Remove the pvalues at zero, otherwise it gives an error!
      # Check that the alleles match
      if(!no_allele_merge){
          match_correct = toupper(merged.data$A1.biom) == toupper(merged.data$A1.eqtl) & toupper(merged.data$A2.biom)== toupper(merged.data$A2.eqtl)
	        match_flip = toupper(merged.data$A1.biom) == toupper(merged.data$A2.eqtl) & toupper(merged.data$A2.biom) == toupper(merged.data$A1.eqtl)
          match_comp_one = toupper(merged.data$A1.biom) == complement_snp(toupper(merged.data$A1.eqtl)) & toupper(merged.data$A2.biom)== complement_snp(toupper(merged.data$A2.eqtl))
     
          match_comp_two = toupper(merged.data$A1.biom) == complement_snp(toupper(merged.data$A2.eqtl)) & toupper(merged.data$A2.biom) == complement_snp(toupper(merged.data$A1.eqtl))

          snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
          print(merged.data[!snp_allele_match,])
          message(sum(snp_allele_match), " SNPs out of ", length(snp_allele_match), " had the correct alleles, discarding SNPs without the correct alleles")
          merged.data = merged.data[snp_allele_match,]
      }
      if (!useBETA) merged.data = merged.data[merged.data$PVAL.biom>0 & merged.data$PVAL.eqtl>0,]
         n_occur <- data.frame(table(merged.data$SNPID))
         dupl = merged.data[merged.data$SNPID %in% n_occur$Var1[n_occur$Freq > 1],]
         message("There are ", nrow(dupl)/2, " duplicated SNP names in the data")
         if (nrow(dupl)>0) {
          dupl=dupl[order(dupl$MAF.eqtl, decreasing=T),]
          toremove = rownames(dupl[ !duplicated(dupl$SNPID), ])
          merged.data = merged.data[!(rownames(merged.data) %in% toremove),]
         }
         nsnps = nrow(merged.data)
         message(ProbeID, ": ", nsnps, " snps in both biomarker and eQTL data. From: ", pos.start, " To: ", pos.end)
         if (nsnps <= min_snps ) {
             message("There are not enough shared snps in the region")
             next
             }else{
         if(!is.null(bed)){
            merged.ranges = GRanges(seqnames=merged.data$CHR.biom,IRanges(start=merged.data$POS.biom,end=merged.data$POS.biom))
            merged.overlaps = findOverlaps(merged.ranges,bed)
            merged.data$bed_region = bed[merged.overlaps@to]$name 
            split_merged_data = split(merged.data, merged.overlaps@to)
         }
         else{
            split_merged_data = list(merged.data) 
         }
        res.out = data.frame()
        for (i in 1:length(split_merged_data)){
            merged.data = split_merged_data[[i]] 
            if(is.null(merged.data$bed_region)){
                merged.data$bed_region=NA
            }
            nsnps = nrow(merged.data)

           if (!useBETA) {
                dataset.biom = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.biom,
                           N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF.biom)
                dataset.eqtl = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.eqtl,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF.eqtl)
           } else {
                dataset.biom = list(snp = merged.data$SNPID, beta = merged.data$BETA.biom, varbeta= (merged.data$SE.biom)^2,
                           s=merged.data$s1, type = type, MAF=merged.data$MAF.biom,N=merged.data$N.biom, sdY=unique(merged.data$sdY.biom))
                dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$BETA.eqtl, varbeta= (merged.data$SE.eqtl)^2,
                           N = as.numeric(merged.data$N.eqtl), type = "quant", MAF=merged.data$MAF.eqtl)
         }
         (coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12))
         pp0       <- as.numeric(coloc.res$summary[2])
         pp1       <- as.numeric(coloc.res$summary[3])
         pp2       <- as.numeric(coloc.res$summary[4])
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])
         snp.biom <- merged.data[which.min(merged.data$PVAL.biom), "SNPID"]
         snp.eqtl <- merged.data[which.min(merged.data$PVAL.eqtl), "SNPID"]
         min.pval.biom <- min(merged.data$PVAL.biom)
         min.pval.eqtl <- min(merged.data$PVAL.eqtl)
         best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])
         ## Per locus likelihood
         # Take the logsum of the 4 models
         l1 = coloc.res$results$lABF.df1
         l2 = coloc.res$results$lABF.df2
         lsum <- coloc.res$results$internal.sum.lABF # lsum = l1 + l2
         lH0.abf <- 0
         lH1.abf <-  logsum(l1) - log(nsnps)
         lH2.abf <-  logsum(l2) - log(nsnps)
         lH3.abf <- logdiff(logsum(l1) + logsum(l2), logsum(lsum))  - log(nsnps^2)
         lH4.abf <- logsum(lsum) -log(nsnps)
         all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
         message(unique(merged.data$bed_region))
         res.temp = data.frame(ProbeID = ProbeID, Chr = chrom, pos.start=pos.start, pos.end=pos.end, nsnps = nsnps, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, PP0.coloc.priors=pp0, PP1.coloc.priors=pp1, PP2.coloc.priors=pp2, PP3.coloc.priors = pp3, PP4.coloc.priors=pp4, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, plotFiles=NA, files=NA, bed_region=unique(merged.data$bed_region))
         if (save.coloc.output) {
           coloc.out = paste(outfolder,"/",prefix,".coloc.output.perSNP/", sep="")
           if (!file.exists(coloc.out)) dir.create(coloc.out)
           if(all(is.na(merged.data$bed_region))){
           write.table(x=coloc.res$results, file=paste(coloc.out, ProbeID,'_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')
           }else{
           write.table(x=coloc.res$results, file=paste(coloc.out, ProbeID,"_",unique(merged.data$bed_region), '_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')
           res.temp$files= as.character(coloc.out)
           }
         }


         ############# PLOT
         # For now disable.
         if (plot & (pp4 > 0.2 | pp3 >=0.2) & nsnps > 2 & F) {
         # For plotting, input called chr has to be numeric!
         # Make sure this is the case otherwise no plot is produced (because the locusZoom script coerces this value to numeric and NA is produced if not numeric
         chrom = gsub("chr","", chrom)

                pvalue_BF_df = as.data.frame(coloc.res[2])
                #region_name <- paste(ProbeID,'.', biom.names[j], ".chr", chr.name, "_", pos.start, "_", pos.end, sep= '')
                region_name <- paste(prefix, ".", ProbeID, ".chr", chrom, "_", pos.start, "_", pos.end, sep= '')
                pvalue_BF_file <- paste(pval.fld, 'pval_', unique(merged.data$bed_region), '.txt', sep="")

                ### LocusZoom arguments:
                pvalue_BF_df$chr = chrom 
                pvalue_BF_df$pos = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "POS.biom"]
                pvalue_BF_df$locus_snp <- paste(pvalue_BF_df$chr, pvalue_BF_df$pos, sep=":")
                pvalue_BF_df$locus_snp <- paste("chr", pvalue_BF_df$locus_snp, sep="")
                # INDELS FORMAT FOR LOCUSZOOM: chr1:117930794:AAG_A (not rsid)
                pvalue_BF_df$locus_snp <- ifelse(grepl("*[:][:]*", pvalue_BF_df$results.snp), paste("chr", as.character(pvalue_BF_df$results.snp), sep=""), as.character(pvalue_BF_df$locus_snp))
                pvalue_BF_df$results.pvalues.df1 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.biom"]
                pvalue_BF_df$results.pvalues.df2 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.eqtl"]
                image.biom = paste(plot.fld, "/", region_name, '_df1.pdf', sep='')
                image.eqtl = paste(plot.fld, "/", region_name, '_df2.pdf', sep='')

                write.table(x =  pvalue_BF_df , file = pvalue_BF_file, row.names = FALSE, quote = FALSE, sep = '\t')

                message('Output pdf for biomarker: ', image.biom)
                pdf(image.biom, width = 9, height = 9)
                #output of region_ld.ld is in /SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs/ ?
                # If INDEL ALLELES do not match exactly (for ex. are reversed from the reference EUR files in here /cluster/project8/vyp/vincent/toolsVarious/locuszoom/EUR/), skip for now:
                plotted = tryCatch(locuszoom.ms(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] , #rs10877835
                  title = 'A',
                  pvalCol='results.pvalues.df1',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( biomarker P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = chrom,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.biom)

                 message('Output pdf for eQTL: ', image.eqtl)
                pdf(image.eqtl, width = 9, height = 9)
                plotted = tryCatch(locuszoom.ms(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] , #rs10877835
                  title = 'B',
                  pvalCol='results.pvalues.df2',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( expression P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = chrom,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.eqtl)

                 res.temp$plotFiles= paste(as.character(paste0(plot.fld, region_name, "_df1.pdf", sep="")), as.character(paste(plot.fld, region_name, "_df2.pdf", sep="")),sep=",")
         }
         res.out = rbind(res.out,res.temp)
         #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
     }
    }
      #if(nrow(res.out)==0){
      #  return(NULL)
      #}
      #return(res.out)
      if(nrow(res.out)!=0){
        res.all = rbind(res.all, res.out)
      }
}

   outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
   res.all <- data.frame(res.all)
   res.all$ProbeID <- as.character(res.all$ProbeID)
   res.all$snp.eqtl <- as.character(res.all$snp.eqtl)
   res.all$best.causal <- as.character(res.all$best.causal)
   res.all$plotFiles <- as.character(res.all$plotFiles)
   res.all$files <- as.character(res.all$files)
   optim.res =  paste(outfolder, 'maximization_results.txt', sep='') 
   # Optimize to find the best parameters
   lkl.frame = res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")]
   alphas = optim(c(2,-2,-2,-2), fn, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
   optim.alphas = exp(alphas$par)/ sum(exp(c(alphas$par,alphas$par[2] + alphas$par[3])))
   write(paste("Model with 4 parameters: ", prefix, ": ", paste(optim.alphas, collapse =" , "), sep=""), file = optim.res, append=TRUE)
       
   lkl.frame = res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")]

 alphas = optim(c(2, -2, -2, -2, -2), fn.pw.gwas, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
  optim.alphas.mle= exp(alphas$par)/ sum(exp(alphas$par))
 if(bootstrap){
bootstrap.all <-  foreach(i=1:no_bootstraps, .combine=rbind) %dopar% {
     llk.frame.temp = lkl.frame[sample(nrow(lkl.frame), size=nrow(lkl.frame), replace=T),]
     alphas = optim(c(2, -2, -2, -2, -2), fn.pw.gwas, data=llk.frame.temp, method = "Nelder-Mead", control=list(fnscale=-1),hessian=T)
     optim.alphas = exp(alphas$par)/ sum(exp(alphas$par))
     return(optim.alphas)
	}
    boot_strap.out  <-  paste(outfolder, "boostrap_estimates.txt", sep="")
    write.table(bootstrap.all, file=boot_strap.out, quote=F, row.names=F)
    cis = t(apply(bootstrap.all,2,function(x){ quantile(x, probs=c(0.025,0.975))}))
    ml_estimates = data.frame(low=cis[,1],mle=optim.alphas.mle,hi=cis[,2])
    bootstrap.summary = paste(outfolder, "bootstrap_mle.txt", sep="")
    write.table(ml_estimates,file=bootstrap.summary, quote=F, row.names=F)
 }
    write(paste("Model with 5 parameters: ", prefix, ": ", paste(optim.alphas.mle, collapse =" , "), sep=""), file = optim.res, append=TRUE)

  # compute posteriors using the already computed likelihoods per locus (lH1.abf etc) and the optimized parameters
   
   new.coloc = apply(res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")], 1, function(x) combine.abf.locus(x[1],x[2], x[3], x[4], x[5], a0 = optim.alphas.mle[1], a1 = optim.alphas.mle[2], a2 = optim.alphas.mle[3], a3 = optim.alphas.mle[4], a4 = optim.alphas.mle[5]))
   new.coloc=t(new.coloc)

   res.all = cbind.data.frame(res.all, new.coloc)

   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')

   # If Gene.name is missing, use ensemblID instead, then try to retrieve name from biomaRt. 
   if (length(res.all$ProbeID[grep("ENSG", res.all$ProbeID)]) >0  & !("Gene.name" %in% names(res.all))) addGeneName = TRUE
   addGeneName= FALSE
   if (addGeneName) {
   res.all$Gene.name = res.all$ProbeID

   biomart=FALSE
      if (biomart) {
      library(biomaRt)
        mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        res.gn <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = as.character(res.all$Gene.name[grep("ENSG", res.all$Gene.name)]), mart = mart)
        res.gn = res.gn[res.gn$hgnc_symbol!="",]
        res.all$Gene.name = res.gn[match(res.all$ProbeID, res.gn$ensembl_gene_id),"hgnc_symbol"]
     } else {
        geneFileNames = "/sc/orga/projects/roussp01a/resources/Ensembl2HGNC/ENSEMBL_v70_TO_HGNC.tsv"
        genes = read.table(geneFileNames, header=F, stringsAsFactors=FALSE, col.names=c("ensembl_gene_id", "hgnc_symbol"))
        res.all$Gene.name = genes[match(res.all$Gene.name, genes$ensembl_gene_id), "hgnc_symbol"]
    }
   }
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
   return(res.all)
}


coloc.eqtl.biom.nopwgwas <- function(eqtl.df, biom.df, p12=1e-6, useBETA=TRUE, plot=FALSE, outfolder, prefix= "pref", save.coloc.output=FALSE, match_snpid=TRUE,cores=20,bootstrap=F,no_bootstraps=1000, min_snps=50, bed_input_file=NULL){
  if (class(eqtl.df$ProbeID)!="character") stop("When reading the data frame, make sure class of ProbeID in eQTL data is a character")

  source("optim_function.R")
  print(head(eqtl.df))
  print(head(biom.df))
  eqtl.df.rs = eqtl.df[grep("rs",eqtl.df$SNPID),]
  print(head(eqtl.df.rs)) 

# Estimate trait variance. 

if (!file.exists(outfolder)) dir.create(outfolder)
if (plot) {
   plot.fld = paste(outfolder, "plot/", sep="")
   pval.fld= paste(outfolder, "pval/", sep="")
   dir.create(file.path(plot.fld), showWarnings = FALSE)
   dir.create(file.path(pval.fld), showWarnings = FALSE)

   # For plotting with locuszoom
   refFlat_path = "/hpc/users/giambc02/scripts/locuszoom/refFlat.RData"
   source("/hpc/users/giambc02/scripts/locuszoom/call_locuszoom3_temp2.R")
   load(refFlat_path)
   refFlatRaw <- refFlatRaw.VP
}
###
  # Set Variables 
  maf_filter = 0.001 # 0.05  #MAF filter applied to datasets
  rsq_filter = 0.6 #Imputation quality filter applied to datasets

###
if (unique(biom.df$type) == "cc") cc=TRUE else cc=FALSE
#if (all(c("CHR", "POS") %in% names(biom.df))) haveCHRPOS.biom=TRUE else haveCHRPOS.biom=FALSE
#if (all(c("CHR", "POS") %in% names(eqtl.df))) haveCHRPOS.eqtl=TRUE else haveCHRPOS.eqtl=FALSE
maf.eqtl = ifelse("MAF" %in% names(eqtl.df), TRUE, FALSE)
maf.biom = ifelse("MAF" %in% names(biom.df), TRUE, FALSE)
if (!maf.eqtl & !maf.biom) message("There is no MAF information in neither datasets, looking for frequency column") 

## check all columns exist
#if (useBETA) cols.eqtl = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "ProbeID", "N") else cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N")
#if (!all(  cols.eqtl %in% names(eqtl.df))) stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
#if (useBETA) cols.biom = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "N") else cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N")
#if (cc) cols.biom = c(cols.biom, "Ncases")
#if (!all(  cols.biom %in% names(biom.df))) stop("These columns are missing from the biomarker data: ", cols.biom[!cols.biom %in% names(biom.df)])

#if ("PVAL" %in% names(biom.df))

if (useBETA) {
   cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE", "ProbeID","N","A1","A2") # We need the N only if we do the sdYest step...
   cols.biom = c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE","N","type","A1","A2")# We need the N only if we do the sdYest step...
   }
if (!useBETA) {
   cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N","A1","A2")
   cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N","A1","A2")
   }

if (cc) cols.biom = c(cols.biom, "Ncases")
no_allele_merge =F
if (!all(cols.eqtl %in% names(eqtl.df))){
    if(!all(cols.eqtl[-which(cols.eqtl %in% c("A1","A2"))] %in% names(eqtl.df))){
        print(cols.eqtl[-which(cols.eqtl %in% c("A1","A2"))])
        stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
    }else{
        no_allele_merge = T
        cols.eqtl = cols.eqtl[-which(cols.eqtl %in% c("A1","A2"))]

        message("Warning: allele columns missing, will sometimes merge SNPs with indels")
    }
}


if (!all(  cols.biom %in% names(biom.df))){ 
    if(!all( cols.biom[-which(cols.biom %in% c("A1","A2"))] %in% names(biom.df))){
        print(all(c("A1","A2") %in% names(biom.df)))
        stop("These columns are missing from the biomarker data: ", cols.biom[!cols.biom %in% names(biom.df)])
    }else{
        no_allele_merge = T
        cols.biom = cols.biom[-which(cols.biom %in% c("A1","A2"))]
        message("Warning: allele columns missing, will sometimes merge SNPs with indels")
    }
}


#####################
# Filter by imputation quality if column exists
info.columns <- grep( names(biom.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    biom.df = subset(biom.df, biom.df[,info.columns] > rsq_filter)
    }
info.columns <- grep( names(eqtl.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    eqtl.df = subset(eqtl.df, eqtl.df[,info.columns] > rsq_filter)
    }

# use only one of the MAFs from the two datasets
# First check if there is a MAF in eQTL data and use this, if not take the one in biom data
# Filter by MAF



if (maf.eqtl) {
   cols.eqtl = c(cols.eqtl, "MAF")
    maf.eqtl = TRUE
   eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
}else if("F" %in% names(eqtl.df)){
    eqtl.df$MAF = ifelse(eqtl.df$F<0.5, eqtl.df$F, 1-eqtl.df$F)
    eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
    maf.eqtl = TRUE
    cols.eqtl = c(cols.eqtl, "F")
}else if (maf.biom) {
   cols.biom = c(cols.biom, "MAF")
   biom.df = subset(biom.df, biom.df$MAF > maf_filter)
   maf.biom= TRUE
}else if("F" %in% names(biom.df)){
   biom.df$MAF = ifelse(biom.df$F<0.5, biom.df$F, 1-biom.df$F)
   biom.df = subset(biom.df, biom.df$MAF > maf_filter)
   maf.biom= TRUE
   cols.biom = c(cols.biom, "F")
}

if (!maf.eqtl & !maf.biom) stop("There is no MAF information in either dataset")


################## SNPID MATCHING
# The reasoning here is that the user can give either only "SNPID", or "SNPID" and "input_names" (or we can find it as rsid or chr:pos and add it in the data as input_names)
# if there is no "input_names" column and the SNPID is not already chr:pos format, then the software will find the chr:pos and see if it matches better with the eQTL data than the SNPID provided

# if there is a "chr" in front of SNPID take out:
  #hasChrSNPID=ifelse(any(grep("chr", biom.df$SNPID[1:100]))>0, TRUE,FALSE)
  #  if (hasChrSNPID) biom.df$SNPID = gsub("chr", "", biom.df$SNPID)
  #hasCharSNPID=ifelse(any(grep("_|-|.", biom.df$SNPID[1:100]))>0, TRUE,FALSE)
  #  if (hasCharSNPID) biom.df$SNPID = gsub("_|-|.", ":", biom.df$SNPID)

  # if there is a "chr" in front of CHR column
  hasChr=ifelse(any(grep("chr",eqtl.df$CHR[1:2]))>0, TRUE,FALSE)
    if (hasChr) (eqtl.df$CHR=gsub("chr", "", eqtl.df$CHR))
  hasChr=ifelse(any(grep("chr",biom.df$CHR[1:2]))>0, TRUE,FALSE)
    if (hasChr) (biom.df$CHR=gsub("chr", "", biom.df$CHR))

  #if (!"input_name" %in% colnames(biom.df) && (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biom.df$SNPID))!=nrow(biom.df))) {
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biom.df$SNPID))!=nrow(biom.df)) addChrposBiom = TRUE
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", eqtl.df$SNPID))!=nrow(eqtl.df)) addChrposEQTL = TRUE

  if (addChrposBiom) {
    biom.df$chrpos = paste(biom.df$CHR, biom.df$POS, sep=":")
  }

  #if (!"input_name" %in% colnames(eqtl.df)) {
  if (addChrposEQTL) {
    eqtl.df$chrpos = paste(eqtl.df$CHR, eqtl.df$POS, sep=":")
  }

if (!match_snpid) {
# Find the combinations of SNPID that matches the most SNPs between the two datasets
biomSNPID = unique(biom.df$SNPID)
eqtlSNPID = unique(eqtl.df$SNPID)
match_snpid = max(length(biomSNPID[biomSNPID %in% eqtlSNPID]), length(eqtlSNPID[eqtlSNPID %in% biomSNPID]))
match_chrpos_snpid = 0
match_chrpos = 0

if (addChrposBiom) {
  biomchrpos = unique(biom.df$chrpos)
  match_chrpos_snpid = max(length(biomchrpos[biomchrpos %in% eqtlSNPID]), length(biomchrpos[biomchrpos %in% eqtlSNPID]))
}

if (addChrposEQTL) {
   eqtlchrpos = unique(eqtl.df$chrpos)
   if (!addChrposBiom) {
   match_chrpos_snpid = max(length(eqtlchrpos[eqtlchrpos %in% biomSNPID]), length(eqtlchrpos[eqtlchrpos %in% biomSNPID]))
   }
}

if (addChrposBiom & addChrposEQTL) match_chrpos = max(length(biomchrpos[biomchrpos %in% eqtlchrpos]), length(biomchrpos[biomchrpos %in% eqtlchrpos]))

# match_snpid = max(length(intersect(biomSNPID, eqtlSNPID)), length(intersect(eqtlSNPID, biomSNPID)))
# Match is faster than intersect
find_best_column = which.max(c(match_snpid, match_chrpos_snpid, match_chrpos))

if (find_best_column==1) message("Best combination is SNPID: do not change column names")
if (find_best_column==2) {
    message("Best combination is SNPID in one dataset and input_name in the other: change one of the column names from input_name to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="chrpos"] <- "SNPID"
}
if (find_best_column==3) {
    message("Best combination is chrpos in both the datasets: change both of the column names from chrpos to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="chrpos"] <- "SNPID"
    names(eqtl.df)[names(eqtl.df)=="SNPID"] <- "SNPID2"
    names(eqtl.df)[names(eqtl.df)=="chrpos"] <- "SNPID"
}
} # if !match_snpid
#####################

  biom.df = biom.df[,cols.biom]
  if("F" %in% names(biom.df)){
    colnames(biom.df)[which("F" == names(biom.df))] = "MAF.biom"
  }   
  if("MAF" %in% names(biom.df)){
    colnames(biom.df)[which("MAF" == names(biom.df))] = "MAF.biom"
  }
  eqtl.df = eqtl.df[,cols.eqtl]
  if("F" %in% names(eqtl.df)){
    colnames(eqtl.df)[which("F" ==  names(eqtl.df))] = "MAF.eqtl"
  }   
  if("MAF" %in% names(eqtl.df)){
    colnames(eqtl.df)[which("MAF" == names(eqtl.df))] = "MAF.eqtl"
  }
  # Remove missing data
  eqtl.df = eqtl.df[complete.cases(eqtl.df),]
  biom.df = biom.df[complete.cases(biom.df),]

  res.all <- data.frame()

###################
# Find ProbeIDs overlapping with biom.df
if(!is.null(bed_input_file)){
    message("Reading LD independent bed file")
    bed = import.bed(bed_input_file)
    bed = bed[bed$ProbeID %in% unique(eqtl.df$ProbeID),]
}else{
   library(data.table)
   DT <- as.data.table(eqtl.df)
   expr_table <- DT[, list(CHR= unique(CHR), START = min(POS), STOP = max(POS), minP = min(PVAL)), by = ProbeID]
   expr_table <- data.frame(expr_table)
   message("There are ", nrow(expr_table), " ProbeIDs in the eQTL data")

   decreaseGap = FALSE
   if (decreaseGap) {
      targetGap=100000
      currentGap = (expr_table$STOP - expr_table$START)/2
      if (targetGap < currentGap[1]) {
        expr_table$START = expr_table$START + currentGap
        expr_table$STOP = expr_table$STOP - currentGap
      }
   }
   bed = expr_table
}

if (!all(c("ProbeID", "CHR", "START", "STOP") %in% names(bed))) stop("Bed file is missing info")
bed$ProbeID = as.character(bed$ProbeID)
bed$CHR=as.numeric(as.character(bed$CHR))
bed$START=as.numeric(bed$START)
bed$STOP=as.numeric(bed$STOP)
message("Looping through ", nrow(bed), " genes from the eQTL data")

##################################################### now start the loop
# Now go over all regions that overlap between eQTL table and input.data

overlap.df<-merge(biom.df,eqtl.df,by="SNPID",suffixes=c(".biom",".eqtl"))

biom.df$sdY.biom = sdY.est((overlap.df$SE.biom)^2,overlap.df$MAF.eqtl,overlap.df$N.biom,overlap.df$BETA.biom) 
#biom.df$sdY.biom = sdY.est((biom.df$SE)^2,biom.df$MAF,biom.df$N,biom.df$BETA) 
message("Running in parallel")
registerDoParallel(cores=cores)
list.probes = bed$ProbeID
print(class(eqtl.df))
eqtl.dfByProbe = split(seq(nrow(eqtl.df)), eqtl.df$ProbeID)
print(head(eqtl.df))
print(length(eqtl.dfByProbe[[1]]))
print("eqtl.df rows")
#print(head(eqtl.df$ProbeID))
if(!is.null(bed_input_file)){
    message("Reading LD independent bed file")
    bed = import.bed(bed_input_file)
}else{
    bed = NULL
}

duplicated_snp_list = data.frame()
res.all = data.frame()
for(i in 1:length(list.probes)){
	print(list.probes[[i]])
    #res.all  <-  foreach(i=1:length(list.probes), .combine=merge_results) %dopar% {
       ProbeID = as.character(list.probes[i]) ##the character bit is important for probe names that are numbers
       #region.eqtl <- subset(eqtl.df.chr, ProbeID == as.character(list.probes[i]))
       region.eqtl <- eqtl.df[eqtl.dfByProbe[[ProbeID]],]
       pos.start <- min(region.eqtl$POS)
       pos.end   <- max(region.eqtl$POS)
       #my.chr = unique(region.eqtl$CHR)
       chrom = region.eqtl$CHR[1]
       #matches <- which(biom.df.chr$CHR==my.chr & biom.df.chr$POS > pos.start & biom.df.chr$POS < pos.end )
       matches <- which(chrom == biom.df$CHR & biom.df$POS >= pos.start & biom.df$POS <= pos.end )
       region.biom <- biom.df[matches, ]
       duplicated_snp_list = rbind(duplicated_snp_list, data.frame(ProbeID = ProbeID, data="biom",remove_dupl(region.biom)[[2]]))
       region.biom = remove_dupl(region.biom)[[1]]
       duplicated_snp_list = rbind(duplicated_snp_list, data.frame(ProbeID = ProbeID, data="eqtl",remove_dupl(region.eqtl)[[2]]))
       region.eqtl = remove_dupl(region.eqtl)[[1]]

       # Loop over each biomarker 
       # message(ProbeID, ": ", length(matches), " snps in biomarkers. From: ", pos.start, " To: ", pos.end)
        
      if (cc) {
          type= "cc"
          #  s = proportion of individuals that are cases (cases / N)
          region.biom$s1 = region.biom$Ncases/region.biom$N
         }
      if (!cc & !useBETA) {
          type = "quant"
          region.biom$s1=rep(0.5, length(region.biom$N)) ## This will be ignored since the type is "quant"
         }
         merged.data <- merge(region.biom, region.eqtl, by = "SNPID",  suffixes=c(".biom", ".eqtl"))
         if(!maf.eqtl){
            merged.data$MAF.eqtl = merged.data$MAF.biom
         }
         # Remove the pvalues at zero, otherwise it gives an error!
      # Check that the alleles match
      if(!no_allele_merge){
          match_correct = toupper(merged.data$A1.biom) == toupper(merged.data$A1.eqtl) & toupper(merged.data$A2.biom)== toupper(merged.data$A2.eqtl)
	  match_flip = toupper(merged.data$A1.biom) == toupper(merged.data$A2.eqtl) & toupper(merged.data$A2.biom) == toupper(merged.data$A1.eqtl)
          match_comp_one = toupper(merged.data$A1.biom) == complement_snp(toupper(merged.data$A1.eqtl)) & toupper(merged.data$A2.biom)== complement_snp(toupper(merged.data$A2.eqtl))
     
          match_comp_two = toupper(merged.data$A1.biom) == complement_snp(toupper(merged.data$A2.eqtl)) & toupper(merged.data$A2.biom) == complement_snp(toupper(merged.data$A1.eqtl))

          snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
          print(merged.data[!snp_allele_match,])
          message(sum(snp_allele_match), " SNPs out of ", length(snp_allele_match), " had the correct alleles, discarding SNPs without the correct alleles")
          merged.data = merged.data[snp_allele_match,]
      }
      if (!useBETA) merged.data = merged.data[merged.data$PVAL.biom>0 & merged.data$PVAL.eqtl>0,]
         n_occur <- data.frame(table(merged.data$SNPID))
         dupl = merged.data[merged.data$SNPID %in% n_occur$Var1[n_occur$Freq > 1],]
         message("There are ", nrow(dupl)/2, " duplicated SNP names in the data")
         if (nrow(dupl)>0) {
          dupl=dupl[order(dupl$MAF.eqtl, decreasing=T),]
          #dupl=dupl[order(dupl$MAF.biom, decreasing=T),]
          toremove = rownames(dupl[ !duplicated(dupl$SNPID), ])
          merged.data = merged.data[!(rownames(merged.data) %in% toremove),]
         }
         nsnps = nrow(merged.data)
         message(ProbeID, ": ", nsnps, " snps in both biomarker and eQTL data. From: ", pos.start, " To: ", pos.end)
         if (nsnps <= min_snps ) {
             message("There are not enough shared snps in the region")
             next
             }else{
         if(!is.null(bed)){
            merged.ranges = GRanges(seqnames=merged.data$CHR.biom,IRanges(start=merged.data$POS.biom,end=merged.data$POS.biom))
            merged.overlaps = findOverlaps(merged.ranges,bed)
            merged.data$bed_region = bed[merged.overlaps@to]$name 
            split_merged_data = split(merged.data, merged.overlaps@to)
         }
         else{
            split_merged_data = list(merged.data) 
         }
        res.out = data.frame()
        for (i in 1:length(split_merged_data)){
            merged.data = split_merged_data[[i]] 
            if(is.null(merged.data$bed_region)){
                merged.data$bed_region=NA
            }
            nsnps = nrow(merged.data)

           if (!useBETA) {
                dataset.biom = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.biom,
                           N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF.biom)
                dataset.eqtl = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.eqtl,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF.eqtl)
           } else {
                dataset.biom = list(snp = merged.data$SNPID, beta = merged.data$BETA.biom, varbeta= (merged.data$SE.biom)^2,
                           s=merged.data$s1, type = type, MAF=merged.data$MAF.biom,N=merged.data$N.biom, sdY=unique(merged.data$sdY.biom))
                dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$BETA.eqtl, varbeta= (merged.data$SE.eqtl)^2,
                           N = as.numeric(merged.data$N.eqtl), type = "quant", MAF=merged.data$MAF.eqtl)
                #dataset.eqtl$MAF <-  maf.eqtl[match(merged.data$SNPID, maf.eqtl$snp ) ,"maf"]
         }
         (coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12))
         pp0       <- as.numeric(coloc.res$summary[2])
         pp1       <- as.numeric(coloc.res$summary[3])
         pp2       <- as.numeric(coloc.res$summary[4])
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])
         snp.biom <- merged.data[which.min(merged.data$PVAL.biom), "SNPID"]
         snp.eqtl <- merged.data[which.min(merged.data$PVAL.eqtl), "SNPID"]
         min.pval.biom <- min(merged.data$PVAL.biom)
         min.pval.eqtl <- min(merged.data$PVAL.eqtl)
         best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])
         ## Per locus likelihood
         # Take the logsum of the 4 models
         l1 = coloc.res$results$lABF.df1
         l2 = coloc.res$results$lABF.df2
         lsum <- coloc.res$results$internal.sum.lABF # lsum = l1 + l2
         lH0.abf <- 0
         lH1.abf <-  logsum(l1) - log(nsnps)
         lH2.abf <-  logsum(l2) - log(nsnps)
         lH3.abf <- logdiff(logsum(l1) + logsum(l2), logsum(lsum))  - log(nsnps^2)
         lH4.abf <- logsum(lsum) -log(nsnps)
         all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
         message(unique(merged.data$bed_region))
         res.temp = data.frame(ProbeID = ProbeID, Chr = chrom, pos.start=pos.start, pos.end=pos.end, nsnps = nsnps, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, PP0.coloc.priors=pp0, PP1.coloc.priors=pp1, PP2.coloc.priors=pp2, PP3.coloc.priors = pp3, PP4.coloc.priors=pp4, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, plotFiles=NA, files=NA, bed_region=unique(merged.data$bed_region))
         if (save.coloc.output) {
           coloc.out = paste(outfolder,"/",prefix,".coloc.output.perSNP/", sep="")
           if (!file.exists(coloc.out)) dir.create(coloc.out)
           if(all(is.na(merged.data$bed_region))){
           write.table(x=coloc.res$results, file=paste(coloc.out, ProbeID,'_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')
           }else{
           write.table(x=coloc.res$results, file=paste(coloc.out, ProbeID,"_",unique(merged.data$bed_region), '_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')
           res.temp$files= as.character(coloc.out)
           }
         }


         ############# PLOT
         # For now disable.
         if (plot & (pp4 > 0.2 | pp3 >=0.2) & nsnps > 2 & F) {
         # For plotting, input called chr has to be numeric!
         # Make sure this is the case otherwise no plot is produced (because the locusZoom script coerces this value to numeric and NA is produced if not numeric
         chrom = gsub("chr","", chrom)

                pvalue_BF_df = as.data.frame(coloc.res[2])
                #region_name <- paste(ProbeID,'.', biom.names[j], ".chr", chr.name, "_", pos.start, "_", pos.end, sep= '')
                region_name <- paste(prefix, ".", ProbeID, ".chr", chrom, "_", pos.start, "_", pos.end, sep= '')
                pvalue_BF_file <- paste(pval.fld, 'pval_', unique(merged.data$bed_region), '.txt', sep="")

                ### LocusZoom arguments:
                pvalue_BF_df$chr = chrom 
                pvalue_BF_df$pos = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "POS.biom"]
                pvalue_BF_df$locus_snp <- paste(pvalue_BF_df$chr, pvalue_BF_df$pos, sep=":")
                pvalue_BF_df$locus_snp <- paste("chr", pvalue_BF_df$locus_snp, sep="")
                # INDELS FORMAT FOR LOCUSZOOM: chr1:117930794:AAG_A (not rsid)
                pvalue_BF_df$locus_snp <- ifelse(grepl("*[:][:]*", pvalue_BF_df$results.snp), paste("chr", as.character(pvalue_BF_df$results.snp), sep=""), as.character(pvalue_BF_df$locus_snp))
                pvalue_BF_df$results.pvalues.df1 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.biom"]
                pvalue_BF_df$results.pvalues.df2 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.eqtl"]
                image.biom = paste(plot.fld, "/", region_name, '_df1.pdf', sep='')
                image.eqtl = paste(plot.fld, "/", region_name, '_df2.pdf', sep='')

                write.table(x =  pvalue_BF_df , file = pvalue_BF_file, row.names = FALSE, quote = FALSE, sep = '\t')

                message('Output pdf for biomarker: ', image.biom)
                pdf(image.biom, width = 9, height = 9)
                plotted = tryCatch(locuszoom.ms(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] ,
                  title = 'A',
                  pvalCol='results.pvalues.df1',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( biomarker P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = chrom,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.biom)

                 message('Output pdf for eQTL: ', image.eqtl)
                pdf(image.eqtl, width = 9, height = 9)
                plotted = tryCatch(locuszoom.ms(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] ,
                  title = 'B',
                  pvalCol='results.pvalues.df2',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( expression P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = chrom,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.eqtl)

                 res.temp$plotFiles= paste(as.character(paste0(plot.fld, region_name, "_df1.pdf", sep="")), as.character(paste(plot.fld, region_name, "_df2.pdf", sep="")),sep=",")
         }
         res.out = rbind(res.out,res.temp)
         #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
     }
    }
      #if(nrow(res.out)==0){
      #  return(NULL)
      #}
      #return(res.out)
      if(nrow(res.out)!=0){
        res.all = rbind(res.all, res.out)
        }
    }
	outfname = paste(outfolder, prefix, '_summary.tab', sep='')
    write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
    return(res.all)
}

sniff <- function(file=fname, eqtl = FALSE) {

    # check if there is a header
    if (all(sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)=="character")) header=TRUE else header=FALSE
    # Column called "F" interpreted as a logical, make sure this is not happening:
    if (!header) {
        checkIflogical=sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)
        if (length(checkIflogical[checkIflogical=="logical"])>0) checkIflogical[checkIflogical=="logical"]="character"
        if (all(checkIflogical=="character")) header=TRUE else header=FALSE
    }
    # For cc need: 
    # SNPID, CHR, POS, F, Z, NCASE, NCONTROL or
    # SNPID, CHR, POS, SE, Z
    # but ro compute Z need either pvalue and beta, or beta and SE so mucst import those too
    # Try to find correct columns to extract
    if (!header) stop("There is no header in the data")

    line = read.table(file, header = T, nrow = 100, stringsAsFactors=F)
    if (ncol(line)==1) {
        separator=","
        line = read.table(file, header = T, nrow = 100, stringsAsFactors=F, sep=separator)
    }
   SNPID = which(grepl("snp|SNP|MarkerName", names(line),  ignore.case = T))
    CHR = which(grepl("^CHROM$|^chr$|^CHR$|^Chromosome$|^chr_name$", names(line),  ignore.case = T))
    POS = which(grepl("^pos$|^bp$|^bp_hg19$|^position|^chrom_start$", names(line),  ignore.case = T))
    BETA = which(grepl("^beta$|^b$|^effect$|^or$|OR.A1|BMIadjMainEffects", names(line),  ignore.case = T))
    F = which(grepl("^F$|freq|FRQ|MAF", names(line),  ignore.case = T))[1] # sometimes have 2 (one for each allele), doesn't matter whcih you take for our applications (fGWAS and coloc)
    if (is.na(F)) F=integer(0)
    PVAL = which(grepl("^p$|pval|Pvalue|P.value|P.val|BMIadjMainP", names(line),  ignore.case = T))[1]
    SE = which(grepl("^se|^StdErr$|BMIadjMainSE", names(line),  ignore.case = T))

    if (length(PVAL)>0) {
        #Make sure this is not a heterogenous pvalue, if it is set to 0
        if (any(grepl("het", names(line)[PVAL], ignore.case =T))) {
            message("Pvalue found is a heterogenous p-value. Looking for another pvalue")
            PVAL = PVAL[!grepl("het", names(line)[PVAL], ignore.case =T)]
        }
        if (length(PVAL)==0) {
            PVAL = which(grepl("^p", names(line),  ignore.case = T))
            #print(names(line)[pval])
        }
    }
    print(names(line))
    A1 = which(grepl("A1|Allele1", names(line),ignore.case=T))
    A2 = which(grepl("A2|Allele2", names(line),ignore.case=T))
    if(length(A1) ==0 || length(A2)== 0){ 
        message("Could not find columns representing A1 and A2, continuing anyway")
    }
    if (eqtl) {
        ProbeID = which(grepl("ProbeID", names(line),  ignore.case = T))
        if (length(ProbeID)==0) stop("Column ProbeID is missing")
    }
    # If there is no "chr" column see if I can extract it from the SNP name (only if all SNP names contain chr info)
    # if (length(chr)==0) & all(grepl("chr",line[,snp])) {
    #  message("Try to extract chr info from file") # sapply(strsplit(as.character(basename(files)), "_",fixed = TRUE), "[", 1)

    # Necessary columns
    # if only chr, pos are missing could match with biomaRt to find rsid
    if (all(c(length(CHR), length(POS))==0) & length(SNPID)!=0) print("No chr, pos: must find from rsid")
    #if (any(c(length(snp), length(chr), length(pos), length(effect), length(pval))==0)) stop("Column names ", c("snp", "chr", "pos", "effect", "pval")[c(length(snp), length(chr), length(pos), length(effect), length(pval))==0], " not recognized")
    if (any(c(length(SNPID), length(PVAL))==0)) stop("Column names ", c("snp", "pval")[c(length(snp), length(pval))==0], " not recognized")
    if (any(c(length(SNPID), length(CHR), length(POS), length(BETA), length(SE), length(PVAL))>1)) stop("Some values match more than one column name: change complicated column names")

    # If have a sample size for each SNP, output this
    N = which(grepl("^N$|^TotalSampleSize$|^SAMPLE_SIZE$", names(line),  ignore.case = T))
    Ncases =  which(grepl("^N_CASES$|^N_CASE$|^Ncases$", names(line),  ignore.case = T)) 
    Ncontrols =  which(grepl("^N_CONTROLS$|^N_CONTROL$|^Ncontrols$", names(line),  ignore.case = T)) 

    info = which(grepl("INFO|RSQ", names(line),  ignore.case = T))
    Z = which(grepl("^Z$|zscore", names(line),  ignore.case = T))

    ### Sanity checks
    if (length(PVAL)>0) if ( sum(line[,PVAL] <0 | line[,PVAL]>1) >0 ) message("Invalid Pvalues")  
    if (length(SE)>0) if ( sum(line[,SE] <=0 | line[,SE]=="Inf" | line[,SE]>10) >0 ) message("Invalid SE")
    if (length(info)>0) if ( sum(line[,info] <0 | line[,info]>1) >0 ) message("Imputation information is less than 0 or greater than 1")

    print(line[1:2,])
    message("SNP column: ", names(line)[SNPID])
    message("CHR column: ", names(line)[CHR])
    message("POS column: ", names(line)[POS])
    message("EFFECT column: ", names(line)[BETA])
    message("PVAL column: ", names(line)[PVAL])
    message("Nsnp column: ", names(line)[N])
    message("INFO column: ", names(line)[info])
    message("SE column: ", names(line)[SE])
    message("Z column: ", names(line)[Z])
    message("FREQ column: ", names(line)[F])
    message("Ncases column: ", names(line)[Ncases])
    message("Ncontrols column: ", names(line)[Ncontrols])
    message("A1 column: ", names(line)[A1])
    message("A2 column: ",names(line)[A2])
    if (eqtl) message("ProbeID column: ", names(line)[ProbeID])

    output_cols = c(SNPID=SNPID, CHR=CHR, POS=POS, BETA=BETA, PVAL=PVAL, N=N, info=info, SE=SE, Z=Z, F=F, Ncases=Ncases, Ncontrols=Ncontrols,A1=A1,A2=A2)
    if (eqtl) output_cols = c(SNPID=SNPID, CHR=CHR, POS=POS, BETA=BETA, PVAL=PVAL, N=N, info=info, SE=SE, Z=Z, F=F, Ncases=Ncases, Ncontrols=Ncontrols, ProbeID=ProbeID,A1=A1,A2=A2)
    # Also output a list of the classes
    output_class = c(class(line[,SNPID]), class(line[,CHR]),class(line[,POS]),class(line[,BETA]),class(line[,PVAL]),class(line[,N]),class(line[,info]),class(line[,SE]),class(line[,Z]),class(line[,F]), class(line[,Ncases]), class(line[,Ncontrols]), class(line[,A1]), class(line[,A2]))
    if (eqtl) output_class = c(class(line[,SNPID]), class(line[,CHR]),class(line[,POS]),class(line[,BETA]),class(line[,PVAL]),class(line[,N]),class(line[,info]),class(line[,SE]),class(line[,Z]),class(line[,F]), class(line[,Ncases]), class(line[,Ncontrols]), class(line[,ProbeID]),class(line[,A1]),class(line[,A2]))
    # if data.frame make it NA for missing columns
    output_class[which(output_class=="data.frame")]=rep("NA", length(which(output_class=="data.frame")))
    # remove missing
    output_class = output_class[output_class!="NA"]
    # if integer make it numeric
      output_class[which(output_class=="integer")]=rep("numeric", length(which(output_class=="integer")))

      output=list(output_cols, output_class)
      return(output)
   }
  # check if there is a header

formatColoc <- function(fname = fname, type="cc", N=NA, Ncases=NA, info_filter=0.6, maf_filter=0.001, fread=T, eqtl=FALSE) { 
# read sniff and quickRead functions from here: 
 library(data.table)
 colsAll = sniff(file=fname, eqtl=eqtl)
 # filter by HetDf > 3 before importing?
 # if ("HetDf" %in% names(colsAll[[1]])) {
 #    colToBeFiltered = colsAll[[1]][which(names(colsAll[[1]])=="HetDf")]
 #    data = quickRead(file = fname, importCols = colsAll[[1]], colClasses = colsAll[[2]], addFilter=TRUE, colToBeFiltered=colToBeFiltered, valueToBeFiltered=3)
  #data = fread(fname, select = as.numeric(colsAll[[1]]), col.names=names(colsAll[[1]]), colClasses = colsAll[[2]])
  # fread select reads the columns in order so must re-order
  colsAll[[1]] = colsAll[[1]][order(colsAll[[1]])]
  #sel <- as.numeric(colsAll[[1]])
  #types of all columns
  data= fread(fname, select = colsAll[[1]], col.names=names(colsAll[[1]]), stringsAsFactors=FALSE)
  data = as.data.frame(data)
  if ( is.null(data ) )  {
     cat("Enter correct column with pvalues.\n")
  }
  if(type == "cc"){
  beta_idx = which(colnames(data) == "BETA")
  }

  if (type=="quant") {
  if (!("N" %in% names(colsAll[[1]]))) {
     #if (is.na(N)) stop("Please specify a sample size because it is not found within the data")
     if (!is.na(N)) {
     message("Sample size not found within the data; using ", N)
     data$N = N
     }
   }
  }
  if (type=="cc") {
  if (!all((c("Ncases", "Ncontrols") %in% names(colsAll[[1]])))) {
     #if ( is.na(Ncases) ) stop("The dataset is a case-control so must specify number of cases")
     if (!is.na(Ncases)) {
     message("Sample size not found within the data; using ", Ncases, " ", N)
     data$Ncases = Ncases
     data$N = N
     #data$Ncontrols = Ncontrols
     #data$N = data$Ncases + data$Ncontrols
     }
    }
  }

  if (length(which(names(data)=="info"))>0 & info_filter!=0 & !is.na(info_filter)) {
    data = data[data$info > info_filter,]
  } else {
  message("Cannot find info so cannot filter")
  }

  if (length(which(names(data)=="F"))>0 & maf_filter!=0 & !is.na(maf_filter)) {
    # could be either MAF or a frequency, but to filter find MAF
    data$F = as.numeric(as.character(data$F))
    data$MAF = ifelse(data$F<0.5, data$F, 1-data$F)
    data = data[data$MAF > maf_filter,]
   } else {
   message("Cannot find frequency so cannot filter")
  }
  
  # you need either effect and SE or pval for coloc
  if ( all(!(c("BETA", "SE") %in% names(data))) | !("PVAL" %in% names(data)) ) stop("Need either BETA and SE, or PVAL in the data")
  #outfile=paste(basename(fname), ".formatted.txt", sep="")
  #write.table(x = data, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  return(data)
}
