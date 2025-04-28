
library(dplyr)

alns <- read.table("hap2ref.table",header=F,stringsAsFactors=F)
colnames(alns) <- c("query","query_l","query_start",
                    "query_end","strand","target",
                    "target_l","target_start","starget_end",
                    "matches","aln_l","qual")
alns$diff <- abs(alns$query_l - alns$aln_l)
alns_filt <- alns %>%
  group_by(query) %>%
  slice(which.min(diff)) %>%
  as.data.frame()
write.table(alns_filt,"hap2ref_filt.txt",sep="\t",row.names=F,col.names=F,quote=F)