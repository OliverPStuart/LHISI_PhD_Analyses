
### Plotting the raw variants from Freebayes on OmniC alignment

# Libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)

# Get the data
data <- read.table("VariantData.txt",header=F,stringsAsFactors=F)
colnames(data) <- c("Scaffold","Position","Site_Quality","Balance","Depth","Genotype","Type")

# Anything odd in there needs to go
Acceptable <- c("0/0","0/1","1/1")
data <- data[data$Genotype %in% Acceptable,]

# Order the scaffold factor
Scaf_Order <- paste0("Scaffold_",c(1:17))
data$Scaffold <- factor(data$Scaffold,
                        levels=Scaf_Order,
                        labels=Scaf_Order)

# Explore data at various resolutions
# Only care about heterozygous sites

Het_SNPs <- data %>% filter(Genotype == "0/1",Type == "snp")
Het_SNPs$Quality_Cutoffs <-cut(Het_SNPs$Site_Quality,breaks=c(0,30,50,100,200,300,1000,100000))

Depths <- Het_SNPs %>% 
  ggplot() + geom_histogram(aes(x=Depth,fill=Quality_Cutoffs)) + 
  scale_x_continuous(trans="log10") + facet_wrap(~Scaffold)

Quals <- Het_SNPs %>% 
  ggplot() + geom_histogram(aes(x=Balance,fill=Quality_Cutoffs)) + 
  facet_wrap(~Scaffold)

# Let's plot their distribution? How to do that
# Rolling window of SNPs per 100kb

# Read in file with contig lengths
Lengths <- read.table("majorScaffolds.bed",header=F,stringsAsFactors = F)[,c(1,3)]
colnames(Lengths) <- c("Scaffold","Length")
Lengths <- Lengths[order(Lengths$Length,decreasing = T),]

Q0 <- Het_SNPs %>% filter(Site_Quality > 0) %>% xtabs(~Scaffold,.) %>% as.data.frame()
Q30 <- Het_SNPs %>% filter(Site_Quality > 30) %>% xtabs(~Scaffold,.) %>% as.data.frame()
Q50 <- Het_SNPs %>% filter(Site_Quality > 50) %>% xtabs(~Scaffold,.) %>% as.data.frame()
Q100 <- Het_SNPs %>% filter(Site_Quality > 100) %>% xtabs(~Scaffold,.) %>% as.data.frame()
Q200 <- Het_SNPs %>% filter(Site_Quality > 200) %>% xtabs(~Scaffold,.) %>% as.data.frame()
Q300 <- Het_SNPs %>% filter(Site_Quality > 300) %>% xtabs(~Scaffold,.) %>% as.data.frame()

QCutoffs <- cbind(Q0,Q30[,2],Q50[,2],Q100[,2],Q200[,2],Q300[,2],Lengths[,2])
colnames(QCutoffs) <- c("Scaffold","0","30","50","100","200","300","Length")

Pi_Plot <- gather(QCutoffs,-c(Length,Scaffold),value="Variants",key="Quality_Cutoff") %>%
  ggplot(aes(x=as.numeric(Quality_Cutoff),y=Variants/Length,group=Scaffold,colour=Scaffold)) + 
  geom_line() + geom_point() + theme(legend.position="none")
Abs_Plot <- gather(QCutoffs,-c(Length,Scaffold),value="Variants",key="Quality_Cutoff") %>%
  ggplot(aes(x=as.numeric(Quality_Cutoff),y=Variants,group=Scaffold,colour=Scaffold)) + 
  geom_line() + geom_point()

Variants <- Pi_Plot + Abs_Plot


png("Raw_SNP_Depths.png",res=300,width=20,height=12,units="in")
plot(Depths)
dev.off()
png("Raw_SNP_Qualities.png",res=300,width=20,height=12,units="in")
plot(Quals)
dev.off()
png("Raw_SNP_Quality_Attrition.png",res=300,height=7,width=10,units="in")
plot(Variants)
dev.off()
