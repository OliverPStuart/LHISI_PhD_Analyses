
# Load libraries
library(dplyr) ; library(ggplot2) ; library(ggExtra)
library(pheatmap) ; library(grid) ; library(patchwork)


# Bring in file
d <- read.table("Scaffold_PacBio.QUAL_FILT.vcf.table",header=F,stringsAsFactors=F)

# Setup data.frame for plotting
colnames(d) <- c("Scaffold","Position","Ref","Alt","Qual","Depth","Ref_C","Alt_C","Amb_C","Geno","Hap","Original")
d$Bal <- d$Alt_C / (d$Ref_C + d$Alt_C)
d$Genotype <- ifelse(d$Geno == "1/1" | d$Geno == "1|1", "hom", "het")
d$In_Hap <- ifelse(d$Hap == ".","out","in")
Scaf_Order <- paste0("Scaffold_",c(1:17))
d$Scaffold <- factor(d$Scaffold,
                        levels=Scaf_Order,
                        labels=Scaf_Order)

# Plot of allele balance at het and hom sites
png("AllScaffold_Balances.png",res=300,height=10,width=10,units="in")
ggplot(d,aes(x=Genotype,y=Bal)) + 
  geom_boxplot() + 
  facet_wrap(~Scaffold,ncol=8) +
  theme_bw() +
  labs(x="In haplotype block?",y="Allele balance")
dev.off()

