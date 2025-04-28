# Analysing coverage across the sex chromosome with ONT female alignments

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# Read in table

data <- read.table("ONT_Scaffold.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
colnames(data) <- c("Scaffold","Start","End","F_Depth")
data$Mid <- (data$Start + data$End) / 2

data$Scaffold <- factor(data$Scaffold,
                        levels=paste0("Scaffold_",1:17),
                        labels=paste0("Scaffold_",1:17))

# All scaffolds

png("ONT_F_depths.png",res=300,width=8,height=8,units="in")
ggplot(data,aes(x=Mid,y=F_Depth,colour=Scaffold=="Scaffold_4")) + 
  geom_point(alpha=0.01) +
  facet_wrap(~Scaffold,scales="free_x") + 
  coord_cartesian(ylim=c(0,110)) +
  scale_colour_manual(values=c("red3","steelblue3")) +
  theme_bw() +
  theme(legend.position="none",
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
dev.off()

ggplot(data,aes(x=F_Depth)) + 
  geom_histogram() +
  facet_wrap(~Scaffold,scales="free_y") + 
  scale_x_continuous(limits=c(0,110))

# Just sex chromosome

data %>% filter(Scaffold=="Scaffold_4") %>%
  ggplot(aes(x=Mid,y=F_Depth)) + 
  geom_point(alpha=0.1) +
  coord_cartesian(ylim=c(0,110))

data %>% filter(Scaffold=="Scaffold_4") %>%
  ggplot(aes(x=F_Depth)) + 
  geom_histogram() +
  scale_x_continuous(limits=c(0,110))

data %>% filter(Scaffold=="Scaffold_4") %>%
  ggplot(aes(x=Mid,y=F_Depth)) + 
  geom_point(alpha=0.2) +
  coord_cartesian(ylim=c(0,110),xlim=c(1e+08,2e+08))

data %>% filter(Scaffold=="Scaffold_4",F_Depth == 0)

