
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)

# Plotting the individual contig depths

# Read in data

data <- read.table("sexChromContigs_OmniC.regions.bed.gz",
                   header=F,stringsAsFactors = F)
colnames(data) <- c("Contig","Start","End","Depth")

# Plot histogram of depth

png("sexChromContigDepths_hist.png",res=300,width=15,height=15,units='in')
data %>% group_by(Contig) %>%
  filter(max(End) > 100000) %>%
ggplot() + 
  geom_histogram(aes(x=Depth)) +
  facet_wrap(~Contig,scales="free_y") +
  scale_x_continuous(limits=c(5,60)) +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank())
dev.off()

# Plot dot plot of depths

png("sexChromContigDepths_dotPlot.png",res=300,width=15,height=15,units='in')
data %>% group_by(Contig) %>%
  filter(max(End) > 100000) %>%
  ggplot() + 
  geom_point(aes(x=Start,y=Depth),alpha=0.2) +
  facet_wrap(~Contig,scales="free_x") +
  scale_y_continuous(limits=c(5,60)) +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank())
dev.off()
