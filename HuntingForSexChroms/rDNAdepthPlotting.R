
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggallin)

### PACBIO

# Get depth of rDNA positions
rdna <- read.table("PacBio_scaffold4GCandDepth.txt",header=F,stringsAsFactors = F)
colnames(rdna) <-  c("Scaffold","Start","End","AT","GC","Depth")

# Get depth of whole scaffold 4 alignment in 10k bins and save GC lookup table
genome <- read.table("PacBio_scaffoldGCandDepth.txt",header=F,stringsAsFactors = F)
colnames(genome) <- c("Scaffold","Start","End","AT","GC","Depth")
genome$GCbin <- cut(genome$GC,seq(from=0,to=1,by=0.02))
summaryGC <- genome %>% 
  group_by(GCbin) %>%
  summarise(meanDepth=mean(Depth),sdDepth=sd(Depth))
genome <- genome[rdna$Scaffold[1] == genome$Scaffold,]

zD <- c()
for(i in 1:nrow(genome)){
  
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == genome[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == genome[i,7]]
  zD[i] <- (genome[i,6] - meanGC) / sdGC
  
}
genome$depthZscore <- zD

rdna$GCbin <- cut(rdna$GC,seq(from=0,to=1,by=0.02))
zD <- c()
for(i in 1:nrow(rdna)){
  
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == rdna[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == rdna[i,7]]
  zD[i] <- (rdna[i,6] - meanGC) / sdGC
  
}
rdna$depthZscore <- zD

# Plot histogram of depths with indicators for depth of specified bins
depth_plot <- ggplot(genome,aes(x=Depth)) + 
  geom_histogram(bins = 100) + 
  scale_x_continuous(trans="log10",limits=c(20,300)) +
  scale_y_continuous(limits=c(0,2250),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  geom_vline(xintercept=rdna$Depth,colour="springgreen3",lwd=0.5,lty="dashed") +
  ggtitle("Depth at mapping\nlocations of rDNA") +
  coord_flip()

# Plot line plot with intervals of mapping locations highlighted
loc_plot <- ggplot() +
  geom_point(alpha=0.1,data=genome,aes(y=Depth,x=Start)) +
  scale_y_continuous(trans="log10",limits=c(20,300),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,genome$End[nrow(genome)]),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank()) +
  geom_vline(xintercept=rdna$Start,colour="springgreen3",lwd=0.5,lty="dashed") +
  labs(title = "Locations of rDNA mapping\non Scaffold 4",
       xlab="Position",y="Read depth\n(truncated from 20-300)")

png("PacBio_readDepth_rdnaMapping.png",res=300,height=6,width=10,units="in")
loc_plot + depth_plot
dev.off()

### And GC normalised depth

GCdepth_plot <- ggplot(genome,aes(x=depthZscore)) + 
  geom_histogram(bins = 100) + 
  scale_x_continuous(limits=c(-2.5,1.5)) +
  scale_y_continuous(limits=c(0,2100),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  geom_vline(xintercept=rdna$depthZscore,colour="springgreen3",lwd=0.5,lty="dashed") +
  ggtitle("GC normalised depth at\nmapping locations of rDNA") +
  coord_flip()

GCloc_plot <- ggplot() +
  geom_point(alpha=0.1,data=genome,aes(y=depthZscore,x=Start)) +
  scale_y_continuous(limits=c(-2.5,1.5),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,genome$End[nrow(genome)]),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank()) +
  geom_vline(xintercept=rdna$Start,colour="springgreen3",lwd=0.5,lty="dashed") +
  labs(title = "Locations of rDNA mapping\non Scaffold 4",
       xlab="Position",y="GC normalised read depth")

png("PacBio_GCreadDepth_rdnaMapping.png",res=300,height=6,width=10,units="in")
GCloc_plot + GCdepth_plot
dev.off()

### OmniC

# Get depth of rDNA positions
rdna <- read.table("OmniC_scaffold4GCandDepth.txt",header=F,stringsAsFactors = F)
colnames(rdna) <-  c("Scaffold","Start","End","AT","GC","Depth")

# Get depth of whole scaffold 4 alignment in 10k bins and save GC lookup table
genome <- read.table("OmniC_scaffoldGCandDepth.txt",header=F,stringsAsFactors = F)
colnames(genome) <- c("Scaffold","Start","End","AT","GC","Depth")
genome$GCbin <- cut(genome$GC,seq(from=0,to=1,by=0.02))
summaryGC <- genome %>% 
  group_by(GCbin) %>%
  summarise(meanDepth=mean(Depth),sdDepth=sd(Depth))
genome <- genome[rdna$Scaffold[1] == genome$Scaffold,]

zD <- c()
for(i in 1:nrow(genome)){
  
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == genome[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == genome[i,7]]
  zD[i] <- (genome[i,6] - meanGC) / sdGC
  
}
genome$depthZscore <- zD

rdna$GCbin <- cut(rdna$GC,seq(from=0,to=1,by=0.02))
zD <- c()
for(i in 1:nrow(rdna)){
  
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == rdna[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == rdna[i,7]]
  zD[i] <- (rdna[i,6] - meanGC) / sdGC
  
}
rdna$depthZscore <- zD

# Plot histogram of depths with indicators for depth of specified bins
depth_plot <- ggplot(genome,aes(x=Depth)) + 
  geom_histogram(bins = 150) + 
  scale_x_continuous(limits=c(5,250)) +
  scale_y_continuous(limits=c(0,7000),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  geom_vline(xintercept=rdna$Depth,colour="springgreen3",lwd=0.5,lty="dashed") +
  ggtitle("Depth at mapping\nlocations of rDNA") +
  coord_flip() +
  annotate(geom="text",y = 2000,x = 240,label="Additional mapping\nat depth = 19879.67")

# Plot line plot with intervals of mapping locations highlighted
loc_plot <- ggplot() +
  geom_point(alpha=0.1,data=genome,aes(y=Depth,x=Start)) +
  scale_y_continuous(limits=c(5,250)) +
  scale_x_continuous(limits=c(0,genome$End[nrow(genome)]),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank()) +
  geom_vline(xintercept=rdna$Start,colour="springgreen3",lwd=0.5,lty="dashed") +
  labs(title = "Locations of rDNA mapping\non Scaffold 4",
       xlab="Position",y="Read depth)")

png("OmniC_readDepth_rdnaMapping.png",res=300,height=6,width=10,units="in")
loc_plot + depth_plot
dev.off()

### And GC normalised depth

GCdepth_plot <- ggplot(genome,aes(x=depthZscore)) + 
  geom_histogram(bins=200) + 
  scale_x_continuous(limits=c(-1,1)) +
  scale_y_continuous(limits=c(0,4300),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()
        ) +
  geom_vline(xintercept=rdna$depthZscore,colour="springgreen3",lwd=0.5,lty="dashed") +
  ggtitle("GC normalised depth at\nmapping locations of rDNA") +
  coord_flip() + 
  annotate(geom="text",x = 0.8, y = 1000, label="Additional mappings\nat Z=14 and Z=464")

GCloc_plot <- ggplot() +
  geom_point(alpha=0.1,data=genome,aes(y=depthZscore,x=Start)) +
  scale_y_continuous(limits=c(-1,1)) +
  scale_x_continuous(limits=c(0,genome$End[nrow(genome)]),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank()) +
  geom_vline(xintercept=rdna$Start,colour="springgreen3",lwd=0.5,lty="dashed") +
  labs(title = "Locations of rDNA mapping\non Scaffold 4",
       xlab="Position",y="GC normalised read depth")

png("OmniC_GCreadDepth_rdnaMapping.png",res=300,height=6,width=10,units="in")
GCloc_plot + GCdepth_plot
dev.off()



### OmniC MAPQ60

# Get depth of rDNA positions
rdna <- read.table("OmniC_scaffold4GCandDepth_Q60.txt",header=F,stringsAsFactors = F)
colnames(rdna) <-  c("Scaffold","Start","End","AT","GC","Depth")

# Get depth of whole scaffold 4 alignment in 10k bins and save GC lookup table
genome <- read.table("OmniC_scaffoldGCandDepth.txt",header=F,stringsAsFactors = F)
colnames(genome) <- c("Scaffold","Start","End","AT","GC","Depth")
genome$GCbin <- cut(genome$GC,seq(from=0,to=1,by=0.02))
summaryGC <- genome %>% 
  group_by(GCbin) %>%
  summarise(meanDepth=mean(Depth),sdDepth=sd(Depth))
genome <- genome[rdna$Scaffold[1] == genome$Scaffold,]

zD <- c()
for(i in 1:nrow(genome)){
  
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == genome[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == genome[i,7]]
  zD[i] <- (genome[i,6] - meanGC) / sdGC
  
}
genome$depthZscore <- zD

rdna$GCbin <- cut(rdna$GC,seq(from=0,to=1,by=0.02))
zD <- c()
for(i in 1:nrow(rdna)){
  
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == rdna[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == rdna[i,7]]
  zD[i] <- (rdna[i,6] - meanGC) / sdGC
  
}
rdna$depthZscore <- zD

# Plot histogram of depths with indicators for depth of specified bins
depth_plot <- ggplot(genome,aes(x=Depth)) + 
  geom_histogram(bins = 150) + 
  scale_x_continuous(limits=c(5,70)) +
  scale_y_continuous(limits=c(0,1750),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()
        ) +
  geom_vline(xintercept=rdna$Depth,colour="springgreen3",lwd=0.5,lty="dashed") +
  ggtitle("Depth at mapping\nlocations of rDNA") +
  coord_flip() + 
  annotate(geom="text",x= 65, y = 300,label="Additional mapping\nat depth=8039")

# Plot line plot with intervals of mapping locations highlighted
loc_plot <- ggplot() +
  geom_point(alpha=0.1,data=genome,aes(y=Depth,x=Start)) +
  scale_y_continuous(limits=c(5,70)) +
  scale_x_continuous(limits=c(0,genome$End[nrow(genome)]),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank()) +
  geom_vline(xintercept=rdna$Start,colour="springgreen3",lwd=0.5,lty="dashed") +
  labs(title = "Locations of rDNA mapping\non Scaffold 4",
       xlab="Position",y="Read depth)")

png("OmniC_readDepth_rdnaMapping_Q60.png",res=300,height=6,width=10,units="in")
loc_plot + depth_plot
dev.off()

### And GC normalised depth

GCdepth_plot <- ggplot(genome,aes(x=depthZscore)) + 
  geom_histogram(bins=200) + 
  scale_x_continuous(limits=c(-0.8,0.3)) +
  scale_y_continuous(limits=c(0,4300),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()
  ) +
  geom_vline(xintercept=rdna$depthZscore,colour="springgreen3",lwd=0.5,lty="dashed") +
  ggtitle("GC normalised depth at\nmapping locations of rDNA") +
  coord_flip() + 
  annotate(geom="text",x = 0.1, y = 1000, label="Additional mappings\nat Z=62")

GCloc_plot <- ggplot() +
  geom_point(alpha=0.1,data=genome,aes(y=depthZscore,x=Start)) +
  scale_y_continuous(limits=c(-0.8,0.3)) +
  scale_x_continuous(limits=c(0,genome$End[nrow(genome)]),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank()) +
  geom_vline(xintercept=rdna$Start,colour="springgreen3",lwd=0.5,lty="dashed") +
  labs(title = "Locations of rDNA mapping\non Scaffold 4",
       xlab="Position",y="GC normalised read depth")

png("OmniC_GCreadDepth_rdnaMapping_Q60.png",res=300,height=6,width=10,units="in")
GCloc_plot + GCdepth_plot
dev.off()

