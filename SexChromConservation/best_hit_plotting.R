
### Plotting chord diagrams of best hits

# Environment
HOME_DIR="/Volumes/Alter/LHISI"
WORKING_DIR=paste0(HOME_DIR,"/Analyses/SexChromConservation")
setwd(WORKING_DIR)

library(circlize)
library(dplyr)
library(magrittr)

species <- c("Sseri", "Samer", "Scanc", "Sgreg", "Snite", "Spice")

# Raw data for species i
data_raw <- read.table(paste0("best_hits_",species[i],".txt"),header=F,stringsAsFactors=F)

# Filter by minimum evalue and max identity
# Anything left after that, just select random row
# Then tabulate by scaffold, V2
data_filt <- data_raw %>% group_by(V1) %>% filter(V11 == min(V11)) %>%
  filter(V3 == max(V3)) %>%
  sample_n(1) %>%
  ungroup() %>% 
  group_by(V2) %>% 
  summarise(Hits=n()) %>%
  mutate(Species=species[i])
colnames(data_filt)[1] <- c("Scaffold")

# Read in tables of assembly information
tmp <- read.table(paste0(species[i],"_seq_report.txt"),
                  header=F,comment.char="#",sep="\t")[,c(3,4,5)]
colnames(tmp) <- c("Molecule","Type","Scaffold")
tmp$Species <- species[i]
tmp[] <- lapply(tmp, function(x) gsub("na", "minor", x))
assign(species[i],tmp)
rm(tmp)
    
tmp <- merge(data_filt,eval(parse(text=species[i])))
tmp <- data.frame(to="Timema_X",
           from=tmp$Molecule,
           value=tmp$Hits)
chordDiagram(tmp, transparency = 0.3,
             col=colorRamp2(breaks=tmp$value,
                            colors=hcl.colors(n=nrow(tmp),
                                              palette="viridis")),
             grid.col = "grey")


