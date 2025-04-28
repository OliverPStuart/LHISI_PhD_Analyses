
library(tidyr)
library(dplyr)
library(ggplot2)

# Plotting the number and type of variants per gene
# Using the table output by snpEff

data <- read.table(text = gsub("#", "", readLines("snpEff_basic/snpEff_genes.txt")[-1]),
                   header=T,stringsAsFactors = F)
colnames(data) <- colnames(data) %>% 
  gsub("variants_","",.) %>%
  gsub("effect_","",.) %>%
  gsub("_variant","",.)

# Remove splice and UTR variants
data <- data[,-c(grep("UTR",colnames(data)),grep("splice",colnames(data)))]

# Subset this by the known proteins, excluding proteins which might not do anything
interest <- readLines("protsOfInterest.txt")
data$OfInterest <- data$GeneName %in% interest

png("VariantTypes.png",res=300,width=14,height=10,units="in")
gather(data[c(1,5:ncol(data))],-c(GeneName,OfInterest),key=Type,value=Count) %>%
  ggplot() + geom_histogram(aes(x=Count,fill=OfInterest)) +
  facet_wrap(~Type,scales="free") +
  scale_y_continuous(trans="log10")
dev.off()

