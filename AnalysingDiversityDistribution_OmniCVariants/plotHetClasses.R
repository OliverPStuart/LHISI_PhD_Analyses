
# This script plots the point estimates of pi per chromosome for a bunch of different sequence types
# Exons, introns, CDS, coding regions with and without repeats included

library(ggplot2)
library(tidyr)

data <- read.table("het_split_by_site_class.txt",header=F,stringsAsFactors = F)
colnames(data) <- c("Scaffold","Type","Length","Het_Sites")
data$Scaffold <- gsub("Scaffold_","",data$Scaffold)

data$Scaffold <- factor(data$Scaffold,
                        labels=c(1:17)[-4],
                        levels=c(1:17)[-4])

data$Pi <- data$Het_Sites / data$Length

png("Pi_PerClass_PerChrom.png",res=300,width=7,height=4,units="in")
data %>% ggplot(aes(x=Type,y=Pi,colour=Type,group=Type)) + 
  geom_point(position=position_dodge2(width=0.8),size=2.5) + 
  facet_wrap(~Scaffold,scales="free_x") +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  scale_y_continuous(limits=c(1e-6,7e-4))
dev.off()
