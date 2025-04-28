
library(ggplot2)
library(dplyr)
library(tidyr)

# Plotting the locations of genes wrt heterozygosity along the chromosomes

# Getting pi
het <- read.table("500kb_variant_counts.txt",header=F,stringsAsFactors = F)
colnames(het) <- c("Scaffold","Start","End","Hets")
het$pi <- het$Hets / 500000
het$Scaffold <- gsub("Scaffold_","",het$Scaffold)
het$Scaffold <- factor(het$Scaffold,
                       labels=1:17,
                       levels=1:17)

# Getting genes
genes <- read.table("protsOfInterest.gff",header=F,stringsAsFactors = F,sep="\t")[,c(1,3,4,5,9)]
genes <- genes[genes$V3 == "gene",c(1,3,4,5)]
colnames(genes) <- c("Scaffold","Start","End","Annotation")
genes$Scaffold <- gsub("Scaffold_","",genes$Scaffold)
genes <- genes[genes$Scaffold %in% c(1:17)[-4],]
genes$Scaffold <- factor(genes$Scaffold,
                       labels=1:17,
                       levels=1:17)

# Plotting the locations of genes along the chromosome wrt local diversity
ggplot(het) + 
  geom_point(alpha=0.2,aes(x=Start,y=pi)) + 
  geom_vline(data=genes,aes(xintercept=Start),
             linetype="dashed",
             colour="springgreen3",
             alpha=0.2) +
  facet_wrap(~Scaffold,scales="free_x") +
  scale_y_continuous(trans="log10") +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

# Okay, so that plot is not useful'
# Plotting the locations of genes in the distribution of pi values

# A loop to get, for each gene, the heterozygosity of the bin it falls into

for(i in 1:nrow(genes)){

  if(i == 1){
    pis <- c()
  }
    
  start.g=genes$Start[i]
  end.g=genes$End[i]
  scaf=genes$Scaffold[i]

  # Get pi at the windows overlaping the start and end coordinates
  start.p <- het$pi[het$Scaffold == scaf & het$Start < start.g & het$End > start.g]
  end.p <- het$pi[het$Scaffold == scaf & het$Start < end.g & het$End > end.g]
  
  # Catch cases where the start/end is on the boundary
  if(length(start.p) == 0)(start.p = end.p)
  if(length(end.p) == 0)(end.p = start.p)
  
  # If the gene is spread across two bins, we take the average
  # If the gene is in one bin, this line changes nothing
  pis[i] <- mean(start.p,end.p)
  
  if(i == nrow(genes)){
    genes$pi <- pis
  }
  
}

ggplot() + 
  geom_histogram(data=het,aes(x=pi)) + 
  geom_vline(data=genes,aes(xintercept=pi),
             linetype="dashed",
             colour="springgreen3",
             alpha=0.5) +
  facet_wrap(~Scaffold,scales="free_y") +
  scale_x_continuous(trans="log10") +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank())

# Also not super informative
# We could also pull out very specific genes if we want and see where they are
# We have a list of insect immune related genes

immune <- readLines("immuneGeneNames.txt")

# Make a column for gene IDs in the genes table
genes$ID <- gsub("^.*Similar to (.*):.*$","\\1",genes$Annotation)

for(i in 1:length(immune)){
  
  if(i == 1){
    genes$Match <- 0
    }
  
  genes$Match[grep(immune[i],genes$Annotation)] <- 1
  
}

ggplot() + 
  geom_histogram(data=het,aes(x=pi)) + 
  geom_vline(data=genes[genes$Match == 1,],aes(xintercept=pi),
             colour="springgreen3",
             lwd=1) +
  facet_wrap(~Scaffold,scales="free_y") +
  scale_x_continuous(trans="log10") +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank())

# Immune genes do not fall inside regions of elevated heterozygosity, unfortunately