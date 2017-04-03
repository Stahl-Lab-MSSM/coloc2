coloc_fld="/u/project/eeskin/pasaniuc/clagiamb/trans_coloc/"
outfname = paste(coloc_fld, "coloc_all", sep="")
traits = list.files(coloc_fld, full.names=T)
library(plyr)

all = data.frame()

for (f in traits) {
 files = list.files(f, full.names=T, pattern="coloc_")
 x = ldply(lapply(files, function(i){read.table(i, header=TRUE)}), data.frame);  
 # 
 library(dplyr)
 x %>% group_by(ProbeID) %>% summarise_each(funs(mean))
 df = data.frame( x[,c("ProbeID","PP0", "PP1","PP2", "PP3", "PP4")] %>% group_by(ProbeID) %>% summarise_each(funs(mean)))
 df$trait = basename(f)
 all = rbind(all, df)
}

all = all[order(all$PP4, decreasing=T),]
write.table(x =  all , file = , row.names = FALSE, quote = FALSE, sep = '\t')

