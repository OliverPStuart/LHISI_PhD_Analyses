
library(dplyr)
library(tidyr)

# Finding reciprocal best hits in TimemaX to LHISISex genes

timema_search <- read.table("SexGenes_to_TimemaX.txt",header=F,stringsAsFactors = F) %>%
  group_by(V1) %>% filter(V11 == min(V11)) %>% as.data.frame

LHISI_search <- read.table("TimemaX_to_SexGenes.txt",header=F,stringsAsFactors = F) %>%
  group_by(V1) %>% filter(V11 == min(V11)) %>% as.data.frame

# Now we just match pairs
timema_pairs <- data.frame(Timema=timema_search[,2],LHISI=timema_search[,1])
LHISI_pairs <- data.frame(Timema=LHISI_search[,1],LHISI=LHISI_search[,2])

# And write the output to a textfile, these are to be excluded from the search for Y-linked genes
writeLines(gsub("-RA","",as.character(unique(inner_join(timema_pairs,LHISI_pairs)$LHISI))),
           "LHISI_TimemaXHomologs.txt")
