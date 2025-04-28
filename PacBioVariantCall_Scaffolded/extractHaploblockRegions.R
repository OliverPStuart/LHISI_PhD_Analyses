
library(magrittr) ; library(dplyr) ; library(ggplot2)

file <- commandArgs(trailingOnly = T)
d <- read.table(file,header=F,stringsAsFactors=F)

# setup data.frame for plotting
colnames(d) <- c("chr","pos","ref","alt","qual","depth","ref_c","alt_c","amb_c","geno","hap","original")
d$genotype <- ifelse(d$geno == "1/1" | d$geno == "1|1", "hom", "het")
d <- d[d$ref_c > 60 & d$alt_c > 60,]

summ <- d[d$hap != ".",] %>% 
  group_by(chr,hap) %>% 
  summarise(het=sum(genotype == "het")/n(),
            len=max(pos)-as.numeric(hap),
            n=n()) %>% 
  filter(row_number()==1,len > 1000) %>% mutate(pi=n/len) %>%
  group_by(chr) %>% top_n(10,wt=len) %>% top_n(3,wt=pi) %>% 
  as.data.frame

write.table(summ,paste0(file,".candidateHaplotypes.txt"),
            quote=F,row.names=F,col.names=F,sep="\t")
