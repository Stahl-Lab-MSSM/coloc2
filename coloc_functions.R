#source("/u/project/eeskin/pasaniuc/clagiamb/scripts/moloc/moloc/R/functions_moloc.R")
#source("/u/project/eeskin/pasaniuc/clagiamb/scripts/moloc/moloc/R/align_alleles_betas.R")
#source("/u/project/eeskin/pasaniuc/clagiamb/scripts/moloc/moloc/R/remove_duplicate_snps.R")

# can either give fixed start and end here, or give start and stop based on any overlap with gwas range
prepare.eqtl.biom <- function(eqtl.df, biom.df, chromosome, start = 1, end = 300*10^6, anyOverlap=FALSE, useBETA=TRUE, outfolder, prefix= "NULL", bychrpos=TRUE, allele_merge=T, flip=TRUE, write=T){

require(data.table)
if (anyOverlap) require(GenomicRanges)


complement_snp <- function(x){
    as = x =="A"
    ts = x == "T"
    gs = x == "G"
    cs = x == "C"
    ins = x == "I"
    dels = x == "D"
    x[as] = "T"
    x[ts] = "A"
    x[gs] = "C"
    x[cs] = "G"
    x[ins] = "NA"
    x[dels] = "NA"
    return(x)
}

remove_dupl = function(data=data, snpcol = "SNP") {
    n_occur <- data.frame(table(data[,snpcol]))
    dupl = data[data[,snpcol] %in% n_occur$Var1[n_occur$Freq > 1],]
    rownames(dupl) = which(data[,snpcol] %in% n_occur$Var1[n_occur$Freq > 1])
         if (nrow(dupl)>0) {
          #removed_list <- rbind(removed_list, data.frame(Marker_removed = dupl$SNPID, reason = "Duplicated SNPs"))
          if (all(c("A1","A2") %in% names(data))) {
             dupl <- transform(dupl, n=nchar(as.character(dupl$A1)) + nchar(as.character(dupl$A2)))
             dupl=dupl[order(dupl$n, decreasing=T),]
          } else if ("MAF" %in% names(data))   {
             dupl=dupl[order(dupl$MAF, decreasing=T),]
          }
          dupl = dupl[ !duplicated(dupl[,snpcol]), ]
          toremove = rownames(data[as.numeric(rownames(dupl)),])
          #if (length(toremove)>0) {
            removed_list <- data.frame(Marker_removed = dupl[,snpcol], reason = "Duplicated SNPs")
          data = data[!(rownames(data) %in% toremove),]
          message("To remove: ", length(toremove), " duplicated SNP names")
          print(removed_list)
          }  else {
          removed_list <- data.frame(Marker_removed = NA, reason = "Duplicated SNPs")
          }
    #return(list(data, removed_list))
    return(data)
}

remove_dupl_biom.data.table <- function(data=data) {
     require(data.table)
     data = data.table(data)
     setkey(data, "SNP")
     key=data[, freq := .N > 1, by = key(data)]
     #unique=data[unique(data),freq:=.N>1]
     dupl = data[duplicated(data[, "SNP", with = FALSE])|duplicated(data[, "SNP", with = FALSE], fromLast=T),]
     if (nrow(dupl)>0) {
          if (all(c("A1","A2") %in% names(data))) {
               toremove = dupl[, crit := nchar(A1) + nchar(A2)][order(crit), .SD[-1L], by="SNP"]
               data=dupl[, crit := nchar(A1) + nchar(A2)][order(crit), .SD[1L], by="SNP"]
             } else {
               toremove = dupl[, crit := MAF][order(crit, decreasing=T), .SD[-1L], by="SNP"]
               data=dupl[, crit := MAF][order(crit, decreasing=T), .SD[1L], by="SNP"]
             }
          message("Removed ", length(toremove), " duplicated SNP names")
          removed_list <- data.frame(Marker_removed = dupl$SNP, reason = "Duplicated SNPs")
          }
     data = data.frame(data)
     data = data[, -which(names(data)%in%c("freq", "crit"))]
     #return(list(data, removed_list))
     return(data)
}


remove_dupl_eqtl.data.table <- function(data=data, byCols=c("ProbeID", "SNP")) {
     require(data.table)
     data = data.table(data)
     setkey(data, "ProbeID", "SNP")
     key=data[, freq := .N > 1, by = key(data)]
     # unique=data[unique(data),freq:=.N>1]
     dupl = data[duplicated(data[, "ProbeID", "SNP", with = FALSE])|duplicated(data[, "ProbeID", "SNP", with = FALSE], fromLast=T),]
          if (all(c("A1","A2") %in% names(data))) {
             data=dupl[, crit := nchar(A1) + nchar(A2)][order(crit), .SD[1L], by=c("ProbeID", "SNP")]
          } else {
             data=dupl[, crit := MAF][order(crit, decreasing=T), .SD[1L], by=c("ProbeID", "SNP")]
          }
          #message("Removed ", length(toremove), " duplicated SNP names")
     data = data.frame(data)
     data = data[, -which(names(data)%in%c("freq", "crit"))]
     return(data)
}

  #####################
  ## check all columns exist
  if ("Ncases" %in% names(biom.df)) cc=TRUE else cc=FALSE
  if (!cc) message("The GWAS is a quantitative trait")
  if (cc) message("The GWAS is a case/control")

  # Remove columns that are all NAs
  biom.df=Filter(function(x)!all(is.na(x)), biom.df)
  eqtl.df=Filter(function(x)!all(is.na(x)), eqtl.df)

  if (useBETA) {
   cols.eqtl = c("SNP", "CHR", "POS", "PVAL", "BETA", "SE", "ProbeID", "N") # Also keep PVAL to report min snp; We need the N only if we do the sdYest step
   cols.biom = c("SNP", "CHR", "POS", "PVAL", "BETA", "SE", "N")# Also keep PVAL to report min snp; We need the N only if we do the sdYest step
   if (cc) cols.biom = c(cols.biom, "Ncases")
  }
  if (!useBETA) {
   cols.eqtl = c("SNP", "CHR", "POS", "PVAL", "ProbeID", "N")
   cols.biom = c("SNP", "CHR", "POS", "PVAL", "N")
   if (cc) cols.biom = c(cols.biom, "Ncases")
  }

  if (!all(  cols.eqtl %in% names(eqtl.df))) stop("These columns are missing from the eQTL data: ", paste(cols.eqtl[!cols.eqtl %in% names(eqtl.df)], collapse=" , "))
  if (!all(  cols.biom %in% names(biom.df))) stop("These columns are missing from the biomarker data: ", paste(cols.biom[!cols.biom %in% names(biom.df)], collapse=" , "))

  if (class(eqtl.df$ProbeID)!="character") stop("When reading the data frame, make sure class of ProbeID in eQTL data is a character")

  if (all(c("A1", "A2") %in% names(eqtl.df)) & all(c("A1", "A2") %in% names(biom.df))) {
    cols.biom = c(cols.biom, "A1", "A2")
    cols.eqtl = c(cols.eqtl, "A1", "A2")
    allele_merge=TRUE 
    message("Matching data frames by alleles: non-matches will be removed")
    } else {
    allele_merge=FALSE
    message("Warning allele columns missing, will sometimes merge SNPs with indels")
  }

  if (!all(c("BETA", "SE") %in% names(biom.df)) & !all(c("BETA", "SE") %in% names(eqtl.df))) useBETA=FALSE


  #####################
  if (!file.exists(outfolder)) dir.create(outfolder)
  
  # now define the output file that contains the eQTL data
  region <- paste('chr', chromosome, '_', start, '_', end, sep = '')
  if (start <= 1 && end >= 300*10^6) region <- paste('chr', chromosome, sep = '')

  if (is.null(prefix)) prefix = region
  outfname = paste(outfolder, prefix, '_summary.tab', sep='')
  out_removed_snps= paste(outfolder, prefix, '_removed_snps.tab', sep='')

  removed_snp_list = data.frame()

  #####################
  # Limit by postion if provided
  # make sure chr, start, end are in numeric format
  chromosome = gsub("chr", "", chromosome)
  chromosome = as.numeric(chromosome)
  start = as.numeric(start)
  end = as.numeric(end)
  
          if (anyOverlap) {
            bed.ranges = GRanges(seqnames=chromosome,IRanges(start=start,end=end))
            eqtl.ranges = GRanges(seqnames=eqtl.df$CHR,IRanges(start=eqtl.df$POS,end=eqtl.df$POS))

            merged.overlaps = findOverlaps(eqtl.ranges,bed.ranges, type="any")
            s = min(eqtl.df[merged.overlaps@queryHits,"POS"])
            e = max(eqtl.df[merged.overlaps@queryHits,"POS"])
            # if there are any genes at the border include
            if (s<start) {
               message("Extending the start")
               start = s
            }
            if (e>end) {
               message("Extending the stop")
               end = e
            }
         }
  
  matches <- which(eqtl.df$CHR==chromosome & eqtl.df$POS > start & eqtl.df$POS < end )
  eqtl.df <- eqtl.df[matches, ]
  matches <- which(biom.df$CHR==chromosome & biom.df$POS > start & biom.df$POS < end )
  biom.df <- biom.df[matches, ]
  if ( nrow(eqtl.df) ==0 | nrow(biom.df) ==0 ) stop ("No SNPs found in the range specified")
  
  if (bychrpos) {
  names(eqtl.df)[which(names(eqtl.df)=="SNP")] = "input"
  names(biom.df)[which(names(biom.df)=="SNP")] = "input"
  biom.df$SNP = paste(biom.df$CHR, biom.df$POS, sep=":")
  eqtl.df$SNP = paste(eqtl.df$CHR, eqtl.df$POS, sep=":")
  }

  if (length(biom.df$BETA[biom.df$BETA<0])>0) log=TRUE  else log=FALSE # if there are negative value, then it is a logOR?
    if (!log) message("GWAS dataset seems to be a case-control and beta of OR needs to be logged: CHECK!!")
    # before taking the log must remove the SNPs with beta = 0
    #if (!log) (listData[[i]] = subset(listData[[i]], BETA != 0))
    #if (!log) (listData[[i]]$BETA = log(listData[[i]]$BETA))

  #####################
  # Set Filters 
  maf_filter = 0.001 # 0.05  #MAF filter applied to datasets
  rsq_filter = 0.3 #Imputation quality filter applied to datasets
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
  maf.eqtl = ifelse(any(c("MAF","F") %in% names(eqtl.df)), TRUE, FALSE) 
  maf.biom = ifelse(any(c("MAF","F") %in% names(biom.df)), TRUE, FALSE)
  if (!maf.eqtl & !maf.biom) {
    message("There is no MAF information in neither datasets, must use external") 
  }
  # First check if there is a MAF in eQTL data and use this, if not take the one in biom data
  # Filter by MAF
  if (maf.eqtl) {
   if ("F" %in% names(eqtl.df) & !("MAF" %in% names(eqtl.df))){
     eqtl.df$MAF = ifelse(eqtl.df$F<0.5, eqtl.df$F, 1-eqtl.df$F)
     }
    eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
    cols.eqtl = c(cols.eqtl, "MAF")
  }
  if (maf.biom) {
   if ("F" %in% names(biom.df) & !("MAF" %in% names(biom.df))){
     biom.df$MAF = ifelse(biom.df$F<0.5, biom.df$F, 1-biom.df$F)
    }
    biom.df = subset(biom.df, biom.df$MAF > maf_filter)
    cols.biom = c(cols.biom, "MAF")
  }
   
  #####################  
  # Remove missing data for SNP, PVAL, BETA, SE
  eqtl.df = eqtl.df[,cols.eqtl]
  biom.df = biom.df[,cols.biom]

  if (useBETA) {
  eqtl.df = eqtl.df[complete.cases(eqtl.df[,c("SNP", "BETA", "SE")]),]
  biom.df = biom.df[complete.cases(biom.df[,c("SNP", "BETA", "SE")]),]
  }
  if (!useBETA) {
  eqtl.df = eqtl.df[complete.cases(eqtl.df[,c("SNP", "PVAL")]),]
  biom.df = biom.df[complete.cases(biom.df[,c("SNP", "PVAL")]),]
  }

  # Remove duplicate SNPs in each data: This is too memory intensive!! Do this after when looping through each ProbeID, or remove first duplicated snp
  # biom.df = remove_dupl_biom.data.table(data=biom.df)
  # eqtl.df = remove_dupl_eqtl.data.table(data=eqtl.df)
  biom.df = biom.df[!duplicated(biom.df$SNP),]
  eqtl.df = eqtl.df[!duplicated(eqtl.df[c("SNP", "ProbeID")]),]

eqtl.df$s = 0.5
eqtl.df$type = "quant"
 if (cc) {
          biom.df$type= "cc"
          #  s = proportion of individuals that are cases (cases / N)
          if (!all(c("N", "Ncases") %in% names(biom.df))) stop("Case/control study: need ratio of cases/controls") 
          biom.df$s = biom.df$Ncases/biom.df$N
 }
 if (!cc) {
          biom.df$type = "quant"
          biom.df$s=rep(0.5, length(biom.df$SNP)) ## This will be ignored since the type is "quant"
 }

  if ( nrow(eqtl.df) ==0 | nrow(biom.df) ==0 ) stop ("No SNPs in the datasets after filters applied")

  merged.data <- merge(biom.df, eqtl.df, by = "SNP",  suffixes=c(".biom", ".eqtl"))

  if ( nrow(merged.data) ==0 ) stop ("No common SNPs in the datasets: make sure column SNP has matching SNP ids across datasets")

  #####################
  # SNP MATCHING
  # Check that the alleles match

      # Check that the alleles match
      if (allele_merge){
         
          # make sure all indels are in the I/D format
          merged.data$A2.biom[nchar(merged.data$A1.biom)>1] = "D"
          merged.data$A1.biom[nchar(merged.data$A1.biom)>1] = "I"
          merged.data$A1.biom[nchar(merged.data$A2.biom)>1] = "D"
          merged.data$A2.biom[nchar(merged.data$A2.biom)>1] = "I"

          merged.data$A2.eqtl[nchar(merged.data$A1.eqtl)>1] = "D"
          merged.data$A1.eqtl[nchar(merged.data$A1.eqtl)>1] = "I"
          merged.data$A1.eqtl[nchar(merged.data$A2.eqtl)>1] = "D"
          merged.data$A2.eqtl[nchar(merged.data$A2.eqtl)>1] = "I"

          match_correct = toupper(merged.data$A1.biom) == toupper(merged.data$A1.eqtl) & toupper(merged.data$A2.biom)== toupper(merged.data$A2.eqtl)
          match_flip = toupper(merged.data$A1.biom) == toupper(merged.data$A2.eqtl) & toupper(merged.data$A2.biom) == toupper(merged.data$A1.eqtl)
          match_comp_one = toupper(merged.data$A1.biom) == complement_snp(toupper(merged.data$A1.eqtl)) & toupper(merged.data$A2.biom)== complement_snp(toupper(merged.data$A2.eqtl))

          match_comp_two = toupper(merged.data$A1.biom) == complement_snp(toupper(merged.data$A2.eqtl)) & toupper(merged.data$A2.biom) == complement_snp(toupper(merged.data$A2.eqtl))

          if (flip) {
           if (any(which(match_flip)>0)) {
              merged.data <- transform(merged.data, A1.eqtl = ifelse(match_flip, A2.eqtl, A1.eqtl), A2.eqtl = ifelse(match_flip, A1.eqtl, A2.eqtl), BETA.eqtl = ifelse(match_flip, -BETA.eqtl, BETA.eqtl))
           }
          }

          snp_allele_match = match_flip | match_correct | match_comp_one | match_comp_two
          print(merged.data[!snp_allele_match,])
          message(sum(snp_allele_match), " SNPs out of ", length(snp_allele_match), " had the correct alleles, discarding SNPs without the correct alleles")
          if (nrow(merged.data[!snp_allele_match,])>0) {
            #removed_snp_list = rbind(removed_snp_list, data.frame(ProbeID = ProbeID, data="merged", Marker_removed=merged.data[!snp_allele_match,"SNP"], reason="Alleles do not match"))
            merged.data = merged.data[snp_allele_match,]
          }

        } #if (allele_merge)
 
  #####################
  # FINAL STEPS NECESSARY FOR COLOC (OLD) TO WORK
  # if only have one MAF per pair of dataset use that for both
         if( ( !maf.eqtl & maf.biom) | (!maf.biom & maf.eqtl) ){
            merged.data$MAF.eqtl = merged.data$MAF
            merged.data$MAF.biom = merged.data$MAF
         }

         # Remove the pvalues at or equal to zero, otherwise it gives an error!
         toremove=merged.data$PVAL.biom<=0 | merged.data$PVAL.eqtl<=0
         if ( sum(toremove) >0 ) {
             removed_snp_list = rbind(removed_snp_list, data.frame(ProbeID = ProbeID, data="merged", Marker_removed=merged.data$SNP[toremove], reason="PVALs zero"))
             merged.data = merged.data[!toremove,]
             # merged.data = merged.data[merged.data$PVAL.biom>0 & merged.data$PVAL.eqtl>0,]
         }

  #####################
if (write) write.table(x =  merged.data , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
return(merged.data)     
}

#######################################################
#######################################################
# provide merged data
coloc.eqtl.biom <- function(merged.data, list.probes=NULL, useBETA=TRUE, outfolder=".", prefix= "test", min_snps=50, p1=1e-4, p2=1e-4, p12=1e-6, moloc=FALSE, write=T){

  # original coloc
  # source("/u/project/eeskin/pasaniuc/clagiamb/scripts/coloc/original_coloc.R")
  library(coloc) 
  ##' This function calculates the log of the sum of the exponentiated
  ##' logs taking out the max, i.e. insuring that the sum is not Inf
  logsum <- function(x) {
   my.max <- max(x)                              ##take out the maximum value in log form
   my.res <- my.max + log(sum(exp(x - my.max )))
   return(my.res)
  }

  ##' This function calculates the log of the difference of the exponentiated
  ##' logs taking out the max, i.e. insuring that the difference is not negative
  logdiff <- function(x,y) {
   my.max <- max(x,y)                              ##take out the maximum value in log form
   my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
   return(my.res)
  }

  if ("Ncases.biom" %in% names(merged.data)) cc=TRUE else cc=FALSE

  if (useBETA) {
   cols.eqtl = c("SNP", "ProbeID","CHR", "POS", "BETA", "SE", "N", "type", "s") 
   cols.biom = c("SNP", "CHR", "POS", "BETA", "SE","N", "type", "s")
   if (cc) cols.biom = c(cols.biom, "Ncases")
  }
  if (!useBETA) {
   cols.eqtl = c("SNP", "ProbeID", "CHR", "POS", "PVAL", "N", "type", "s")
   cols.biom = c("SNP", "CHR", "POS", "PVAL", "N", "type", "s")
   if (cc) cols.biom = c(cols.biom, "Ncases")
  }
  cols.biom[-1] = paste(cols.biom[-1], ".biom", sep="")
  cols.eqtl[-c(1,2)] = paste(cols.eqtl[-c(1,2)], ".eqtl", sep="")
  if (!all(  cols.eqtl %in% names(merged.data))) stop("These columns are missing from the eQTL data: ", paste(cols.eqtl[!cols.eqtl %in% names(merged.data)], collapse=" , "))
  if (!all(  cols.eqtl %in% names(merged.data))) stop("These columns are missing from the eQTL data: ", paste(cols.eqtl[!cols.eqtl %in% names(merged.data)], collapse=" , "))

  #####################
  # Output
  if (!file.exists(outfolder)) dir.create(outfolder)
  outfname = paste(outfolder, prefix, '_summary.tab', sep='')

  res.all <- data.frame()

  #####################
  # Now go over all regions that overlap between eQTL table and input.data
   # list.probes = bed$ProbeID
   if (is.null(list.probes)) list.probes = names(table(merged.data$ProbeID))
   dfByProbe = split(seq(nrow(merged.data)), as.factor(merged.data$ProbeID))

   message("Looping through ", length(list.probes), " genes")

   for (i in 1:length(list.probes)) {

       ProbeID = as.character(list.probes[i]) ##the character bit is important for probe names that are numbers
       message("Looping through ProbeID: ", ProbeID)
       merged.region = merged.data[dfByProbe[[as.character(ProbeID)]],]

       nsnps = nrow(merged.region)
       if (nsnps <= min_snps ) {
             message("There are not enough common snps in the region")
             next()
       } else {

       # For now run with p-values (better for cc data)
       if (!useBETA) {
          dataset.biom = list(snp = merged.region$SNP, pvalues = merged.region$PVAL.biom,
                           N = merged.region$N.biom, s=merged.region$s.biom, type = unique(merged.region$type.biom), MAF=merged.data$MAF.biom)
           dataset.eqtl = list(snp = merged.region$SNP, pvalues = merged.region$PVAL.eqtl,
                           N = merged.region$N.eqtl, s=merged.region$s.eqtl, type = unique(merged.region$type.eqtl), MAF=merged.region$MAF.eqtl)
        } else {
           dataset.biom = list(snp = merged.region$SNP, beta = merged.region$BETA.biom, varbeta= (merged.region$SE.biom)^2, s=merged.region$s.biom, type = unique(merged.region$type.biom), MAF=merged.region$MAF.biom,N=merged.region$N.biom) #, sdY=unique(merged.data$sdY.biom))
           dataset.eqtl = list(snp = merged.region$SNP, beta = merged.region$BETA.eqtl, varbeta= (merged.region$SE.eqtl)^2, N = as.numeric(merged.region$N.eqtl), s=merged.region$s.eqtl, type = unique(merged.region$type.eqtl), MAF=merged.region$MAF.eqtl)
         }

         ### COLOC OLD
         capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12))
         pp0       <- as.numeric(coloc.res$summary[2])
         pp1       <- as.numeric(coloc.res$summary[3])
         pp2       <- as.numeric(coloc.res$summary[4])
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])

         min.pval.biom <- min(merged.region$PVAL.biom)
         min.pval.eqtl <- min(merged.region$PVAL.eqtl)
         min.pval.snp.biom <- merged.region[which.min(merged.region$PVAL.biom), "SNP"]
         min.pval.snp.eqtl <- merged.region[which.min(merged.region$PVAL.eqtl), "SNP"]
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
         coloc.lkl = c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
         coloc.ppa = c(pp0, pp1, pp2, pp3, pp4)

         res.temp = data.frame(ProbeID = ProbeID, locus=prefix, nsnps = nsnps, min.pval.biom = min.pval.biom, min.pval.eqtl=min.pval.eqtl, min.pval.snp.biom = min.pval.snp.biom, min.pval.snp.eqtl = min.pval.snp.eqtl, best.causal=best.causal, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, PP0=pp0, PP1=pp1, PP2=pp2, PP3=pp3, PP4=pp4)
         #res.temp = data.frame(ProbeID = ProbeID, locus=prefix, nsnps = nsnps, best.causal=best.causal, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, PP0=pp0, PP1=pp1, PP2=pp2, PP3=pp3, PP4=pp4)
        all = TRUE
        if (all) {
             res.SNP = data.frame(ProbeID = ProbeID, SNP = coloc.res$results$snp, lABF.df1 = l1, lABF.df2 = l2, internal.sum.lABF = lsum)
             res.temp = merge(res.temp, res.SNP, by="ProbeID")
         }
         
         ### MOLOC
         if (moloc) {
           library(moloc)
           merged.region.biom = merged.region[,cols.biom]
           names(merged.region.biom) = gsub(".biom", "", names(merged.region.biom))
           merged.region.eqtl = merged.region[,cols.eqtl]
           names(merged.region.eqtl) = gsub(".eqtl", "", names(merged.region.eqtl))
         
           moloc <- moloc_test(listData=list(merged.region.biom, merged.region.eqtl), overlap=FALSE, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06), from_p=FALSE)
           # Posteriors
           print(moloc[[1]])
           # Number of SNPs analyzed
           print(moloc[[2]])
           # Posterior of the most likely SNP co-localizing with another trait
           print(moloc[[3]])
           # put the PPA in same order as reported in coloc
           moloc.ppa = c(moloc[[1]]$PPA[5], moloc[[1]]$PPA[1], moloc[[1]]$PPA[3], moloc[[1]]$PPA[2], moloc[[1]]$PPA[4])
           moloc.ppa = paste(moloc.ppa, collapse=",")
           best.causal.moloc = as.character(moloc[[3]]$best.snp.coloc[3])
           res.temp = cbind.data.frame(res.temp, moloc.ppa, best.causal.moloc)
        }

        print(res.temp)
        res.all = rbind.data.frame(res.temp, res.all)
    }
    }
    if (write) {
       write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
    }
    return(res.all)

    #write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
}

#############################################################
#### POSTPROCESSING WITH ALL LOCI
#############################################################

add_fdr <- function(df, snpcol = "PP4") {
    df$PEP <- 1 - df[,snpcol]
    df <- df[order(df$PEP, decreasing=F),]
    df$FDR <- cumsum(df$PEP)/seq(df$PEP)
    message("There are ", nrow(df[df$FDR<=0.05,]), " results at FDR 5%")
    return(df)
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
  sumlkl = sumlkl + sum(log(c(dnorm(p[1], mean=2, sd=3), (dnorm(p[2:5], mean=-2, sd=3)))))
  #print(sumlkl)
  return(sumlkl)
}

logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}


logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}

# After optimization step to find the priors
# Apply this function to every row of data (every row is a locus)
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
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}

# Optimize to find the best parameters
est_lkl <- function(res.all, colnames.lkl = c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf"), bootstrap=TRUE, no_bootstraps=1000, cores=20) {

   lkl.frame = res.all[,colnames.lkl]
   lkl.frame <-as.matrix(sapply(lkl.frame, function(x) as.numeric(as.character(x))))
       
   alphas = optim(c(2, -2, -2, -2, -2), fn.pw.gwas, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
   optim.alphas.mle= exp(alphas$par)/ sum(exp(alphas$par))
   print(optim.alphas.mle)
   message("Model with 5 parameters: ", paste(optim.alphas.mle, collapse =" , "), sep="")

if(bootstrap){
library(foreach)
library(doParallel)
   bootstrap.all <-  foreach(i=1:no_bootstraps, .combine=rbind) %dopar% {
      llk.frame.temp = lkl.frame[sample(nrow(lkl.frame), size=nrow(lkl.frame), replace=T),]
      alphas = optim(c(2, -2, -2, -2, -2), fn.pw.gwas, data=llk.frame.temp, method = "Nelder-Mead", control=list(fnscale=-1),hessian=T)
      optim.alphas = exp(alphas$par)/ sum(exp(alphas$par))
      return(optim.alphas)
     }
     #boot_strap.out  <-  paste(outfolder, "boostrap_estimates.txt", sep="")
     #write.table(bootstrap.all, file="boostrap_estimates.txt", quote=F, row.names=F)
     cis = t(apply(bootstrap.all,2,function(x){ quantile(x, probs=c(0.025,0.975))}))
     ml_estimates = data.frame(low=cis[,1],mle=optim.alphas.mle,hi=cis[,2])
     #bootstrap.summary = paste(outfolder, "bootstrap_mle.txt", sep="")
     #write.table(ml_estimates,file="bootstrap_mle.txt", quote=F, row.names=F)
 }

  # compute posteriors using the already computed likelihoods per locus (lH1.abf etc) and the optimized parameters 
   new.coloc = apply(lkl.frame, 1, function(x) combine.abf.locus(x[1],x[2], x[3], x[4], x[5], a0 = optim.alphas.mle[1], a1 = optim.alphas.mle[2], a2 = optim.alphas.mle[3], a3 = optim.alphas.mle[4], a4 = optim.alphas.mle[5]))
   new.coloc=data.frame(t(new.coloc))
   # Add FDR for PP4
   new.coloc = add_fdr(df=new.coloc, snpcol = "PP.H4.abf")

   # Add priors estimation   
   res.all =cbind.data.frame(res.all, new.coloc, t(data.frame(optim.alphas.mle, row.names = c("prior0", "prior1", "prior2", "prior3", "prior4"))))

   if (bootstrap) {
        CI_low = t(ml_estimates$low)
        colnames(CI_low)=c("CI_low_p0","CI_low_p1","CI_low_p2","CI_low_p3","CI_low_p4")
        CI_hi = t(ml_estimates$hi)
        colnames(CI_hi)=c("CI_hi_p0","CI_hi_p1","CI_hi_p2","CI_hi_p3","CI_hi_p4")

       res.all = cbind.data.frame(res.all, CI_low, CI_hi)
   }

   return(res.all)

}
