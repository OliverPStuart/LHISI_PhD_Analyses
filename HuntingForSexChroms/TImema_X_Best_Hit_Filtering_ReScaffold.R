
# Filtering blast of Timema X chromosome genes

library(dplyr)
library(patchwork)
library(ggplot2)

# Get data
blast_hits <- read.table("Timema_X_to_LHISI_ReScaffold.txt",
                         header=F,
                         stringsAsFactors = F)

# Keep best hit per gene
filtered_blast_hits <- blast_hits %>% 
  group_by(V1) %>%
  filter(V13==min(V13)) %>%
  as.data.frame()

# There are some duplicates in here
# As in, some genes blast with equally low score to more than one chromosome 
# We'll just leave them in, they form a very small minority

# Write a table of the best hits
write.table(data.frame(Gene=filtered_blast_hits$V1,
           Scaffold=filtered_blast_hits$V2,
           Start=filtered_blast_hits$V11,
           End=filtered_blast_hits$V12),
           "Timema_X_Best_Hits_ReScaffold.txt",col.names=T,row.names=F,quote=F,sep="\t")

# Also make a bed file of this table to pass to mosdepth
bed <- data.frame(Scaffold=filtered_blast_hits$V2,
           Start=filtered_blast_hits$V11,
           End=filtered_blast_hits$V12)

# Reorder columns since some hits will be on -ve strand
for(i in 1:nrow(bed)){
  
  if(bed$Start[i] > bed$End[i]){
    Start=min(bed$Start[i],bed$End[i])
    End=max(bed$Start[i],bed$End[i])
    bed$Start[i] <- Start
    bed$End[i] <- End
  }
  
}

write.table(bed,"Timema_X_Best_Hits_ReScaffold.bed",col.names=F,row.names=F,quote=F,sep="\t")
