
# Libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Read in data from command line argument
input <- commandArgs(trailingOnly=T)
raw.data <- read.table(input,header=F,stringsAsFactors=F)

# Select only best match per minor scaffold
# This is the scaffold with the highest number of matched bases over the alignment length

data <- raw.data %>% 
  group_by(V1) %>%
  slice_max(V10) %>%
  mutate(ident=V10/V11, matchLen=(V4-V3)/V2)
data$Scaffold <- gsub("(Scaffold_([0-9]+)).*","\\1",data$V6)

# Plot identity and length of matches
identPlot <- ggplot() +
  geom_point(data=data,aes(x=ident,y=matchLen,colour=Scaffold,size=V2)) +
  scale_size(trans="log10",name = "Minor\nscaffold\nsize") +
  labs(y="% of minor scaffold aligned",x="% matched bases") + 
  scale_x_continuous(limits=c(-0.02,1.02),expand=c(0,0)) +
  scale_y_continuous(limits=c(-0.02,1.02),expand=c(0,0)) +
  theme_bw() +
  scale_colour_discrete(guide=F)

# Plot histogram bins of scaffold length
# Coloured by the proportion in that bin mapped or unmapped
minor <- read.table("minorScaffoldLengths",header=F,stringsAsFactors = F)
minor <- minor %>% 
  mutate(mapped=V1 %in% raw.data$V1)
plot1 <- ggplot() + geom_boxplot(data=minor,aes(y=V2,fill=mapped,x=mapped)) + 
  scale_y_continuous(trans="log10") +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none") +
  labs(y="minor scaffold length")
plot2 <- ggplot() + geom_histogram(data=minor,aes(V2,fill=mapped)) +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(limits=c(0,250),expand=c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())

output <- (plot1 + plot2) / identPlot

png(paste0(input,"_plot.png"),res=300,width=8,height=8,units="in")
plot(output)
dev.off()
