### File for plotting timema homolog histogram

HOME_DIR="/Volumes/Alter/LHISI"
WORKING_DIR=paste0(HOME_DIR,"/Analyses/HuntingForSexChroms")
setwd(WORKING_DIR)

library(dplyr)
library(ggplot2)
library(patchwork)
library(svglite)

# Read in data
hits <- read.table("Timema_X_Best_Hits.txt",header=T,stringsAsFactors = F)
scafOrder <- 1:17
hits$Scaffold <- factor(gsub("Scaffold_","",hits$Scaffold),labels=scafOrder,levels=scafOrder)

p <- hits %>% group_by(Scaffold) %>%
  summarise(Hits=n()) %>%
  ggplot(aes(x=Scaffold,y=Hits,fill=Scaffold=="4")) + geom_col(colour="black") +
  scale_y_continuous(limits=c(0,150),expand=c(0,0)) +
  scale_fill_manual(values=c("palegreen4", "palegreen3")) +
  theme(
    axis.ticks=element_blank(),
    axis.text=element_text(size=7),axis.title=element_text(size=10),
    legend.position = "none",
    panel.grid=element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour="black",fill=NA,size=1)
  ) +
  labs(y="Timema X\nhomologs")


svg("timema_hits.svg",height=4,width=3)
plot(p)
dev.off()

png("timema_hits.png",height=4,width=3,units='in',res=300)
plot(p)
dev.off()

