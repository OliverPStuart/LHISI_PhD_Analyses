
# Analysing coverage across the sex chromosome from RNAseq alignments

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

`%ni%` <- Negate(`%in%`)

# First we read in each individual separately, then merge all of them

F2 <- read.table("F2.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
F3 <- read.table("F3.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
F4 <- read.table("F4.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
M1 <- read.table("M1.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
M2 <- read.table("M2.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
M3 <- read.table("M3.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
M4 <- read.table("M4.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")

All <- cbind(F2[,c(1,2,3,4,5)],F3[,5],F4[,5],M1[,5],M2[,5],M3[,5],M4[,5])
colnames(All) <- c("Scaffold","Start","End","Annotation","F2","F3","F4","M1","M2","M3","M4")

# Let's try and identify locations where females have exceptionally low coverage and males have some

All$NoFemale <- rowMeans(All[,c(5:7)]) == 0
All$Male <- rowMeans(All[,c(8:11)]) > 0
PotentialRegions <- All[All$NoFemale == T & All$Male == T,]

# Let's also only pick those that fall into lower depth bins
coverage <- read.table("/Volumes/Alter/LHISI/Analyses/HuntingForSexChroms/OmniC_scaffoldGCandDepth.txt",
                       header=F,stringsAsFactors = F)
coverage <- coverage[coverage$V1 == "Scaffold_4",c(1,2,3,6)]
colnames(coverage) <- c("Scaffold","Start","End","Depth")

for(i in 1:nrow(PotentialRegions)){
  
  if(i == 1){
    depths <- c()
  }
  
  start.g=PotentialRegions$Start[i]
  end.g=PotentialRegions$End[i]
  scaf=PotentialRegions$Scaffold[i]
  
  # Get depth at the windows overlapping the start and end coordinates
  start.d <- coverage$Depth[coverage$Scaffold == scaf & coverage$Start < start.g & coverage$End > start.g]
  end.d <- coverage$Depth[coverage$Scaffold == scaf & coverage$Start < end.g & coverage$End > end.g]
  
  # Catch cases where the start/end is on the boundary
  if(length(start.d) == 0)(start.d = end.d)
  if(length(end.d) == 0)(end.d = start.d)
  
  # If the gene is spread across two bins, we take the average
  # If the gene is in one bin, this line changes nothing
  depths[i] <- mean(start.d,end.d)
  
  if(i == nrow(PotentialRegions)){
    PotentialRegions$Depth <- depths
  }
  
}

# Now bring in list of viral-like and transposon proteins and remove
trans <- readLines("ViralLikeOrTransposon.txt")
PotentialRegions$Name <- substr(PotentialRegions$Annotation,4,11)
PotentialRegions <- PotentialRegions[PotentialRegions$Name %ni% trans,]

# Now select only those with < 20 coverage and some similarity to a known protein
write.table(PotentialRegions[PotentialRegions$Depth < 18,c(1:4)],
            "PotentialYLinked_lessStringent.bed",col.names=F,row.names=F,quote=F,sep="\t")
