
# Computing GC normalised coverage

library(dplyr)
library(ggplot2)
library(patchwork)

#### PACBIO ALIGNMENT

# Read data
data <- read.table("PacBio_scaffoldGCandDepth.txt",header=F,stringsAsFactors=F)

# Add column names, make GC bins, clean up scaffold ID
colnames(data) <- c("Scaffold","WindowStart","WindowEnd","AT","GC","Depth")
data$GCbin <- cut(data$GC,seq(from=0,to=1,by=0.02))
Scaf_Order <- paste0("Scaffold_",c(1:17))
data$Scaffold <- factor(data$Scaffold,
                        levels=Scaf_Order,
                        labels=Scaf_Order)

# Create summary table for lookup
summaryGC <- data %>% 
  group_by(GCbin) %>%
  summarise(meanDepth=mean(Depth),sdDepth=sd(Depth))

# Empty vector for storage
zD <- c()

# For all rows
for(i in 1:nrow(data)){
  
  # Get the mean and sd of depth in that GC bin
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == data[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == data[i,7]]
  
  # Compute Z-score of depth
  zD[i] <- (data[i,6] - meanGC) / sdGC
  
}

data$depthZscore <- zD

# Plot Z-scores of GC normalisation
# This may help in future identifying duplications
depthPlot_GCnorm <- data %>% 
  ggplot(aes(y=depthZscore,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(-3,3)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="GC normalised depth")

# Plot depth, this is all that's necessary for getting the putative sex chromosome
depthPlot <- data %>% 
  ggplot(aes(y=Depth,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(0,500)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="depth")

# Plot GC
gcPlot <- data %>% 
  ggplot(aes(y=GC,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(0.2,0.6)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="GC %")

# Plot relationship between GC and Depth overall

colours.vec=c("mediumpurple1","purple4")

gcDens <- ggplot(data,aes(x=GC,group=Scaffold,
                          colour=Scaffold=="Scaffold_4",
                          fill=Scaffold=="Scaffold_4")) + 
  geom_density(alpha=0.3) + 
  scale_x_continuous(limits=c(0.2,0.6),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,20),expand=c(0,0)) +
  scale_colour_manual(values=colours.vec) + 
  scale_fill_manual(values=colours.vec) + 
  coord_flip() +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())

gcDepthPlot <- ggplot(data,aes(x=Depth,y=GC,
                               colour=Scaffold=="Scaffold_4")) + 
  geom_point(alpha=0.01) +
  scale_y_continuous(limits=c(0.2,0.6),expand=c(0,0)) +
  scale_x_continuous(limits=c(9,300),expand=c(0,0)) +
  scale_colour_manual(values=colours.vec) + 
  theme_bw() +
  theme(legend.position="none",
        axis.ticks=element_blank())

depthDens <- ggplot(data,aes(x=Depth,group=Scaffold,
                             colour=Scaffold=="Scaffold_4",
                             fill=Scaffold=="Scaffold_4")) +
  geom_density(alpha=0.3) +
  scale_x_continuous(limits=c(9,300),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.035),expand=c(0,0)) +
  scale_colour_manual(values=colours.vec) + 
  scale_fill_manual(values=colours.vec) + 
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        legend.position="none")

fullPlot <- depthDens + plot_spacer() + 
  gcDepthPlot + gcDens + 
  plot_layout(ncol=2,widths=c(3,1),nrow=2,heights=c(1,3))

# Lastly, just plotting depths histograms
PacBio_Depth_Hist <- ggplot(data,aes(x=Depth)) + geom_histogram(bins=50) + facet_wrap(~Scaffold) +
  scale_x_continuous(limits=c(20,300))

# Now save plots

png("PacBio_scaffoldDepths.png",res=300,width=10,height=10,units="in")
plot(depthPlot)
dev.off()

png("scaffoldGC.png",res=300,width=10,height=10,units="in")
plot(gcPlot)
dev.off()

png("PacBio_scaffoldDepths_GCnorm.png",res=300,width=10,height=10,units="in")
plot(depthPlot_GCnorm)
dev.off()

png("PacBio_gcVdepth.png",res=300,height=10,width=10,units="in")
plot(fullPlot)
dev.off()

png("PacBio_Depth_Hist.png",res=300,width=8,height=8,units="in")
plot(PacBio_Depth_Hist)
dev.off()


#### OMNIC ALIGNMENT

# Read data
data <- read.table("OmniC_scaffoldGCandDepth.txt",header=F,stringsAsFactors=F)

# Add column names, make GC bins, clean up scaffold ID
colnames(data) <- c("Scaffold","WindowStart","WindowEnd","AT","GC","Depth")
data$GCbin <- cut(data$GC,seq(from=0,to=1,by=0.02))
Scaf_Order <- paste0("Scaffold_",c(1:17))
data$Scaffold <- factor(data$Scaffold,
                        levels=Scaf_Order,
                        labels=Scaf_Order)


# Create summary table for lookup
summaryGC <- data %>% 
  group_by(GCbin) %>%
  summarise(meanDepth=mean(Depth),sdDepth=sd(Depth))

# Empty vector for storage
zD <- c()

# For all rows
for(i in 1:nrow(data)){
  
  # Get the mean and sd of depth in that GC bin
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == data[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == data[i,7]]
  
  # Compute Z-score of depth
  zD[i] <- (data[i,6] - meanGC) / sdGC
  
}

data$depthZscore <- zD


# Plot Z-scores of GC normalisation
# This may help in future identifying duplications
depthPlot_GCnorm <- data %>% 
  ggplot(aes(y=depthZscore,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(-3,3)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="depth")


# Plot depth, this is all that's necessary for getting the putative sex chromosome
depthPlot <- data %>% 
  ggplot(aes(y=Depth,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(0,155)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="depth")

# Plot relationship betwen GC and Depth overall
gcDens <- ggplot(data,aes(x=GC,group=Scaffold,
                          colour=Scaffold=="Scaffold_4",
                          fill=Scaffold=="Scaffold_4")) + 
  geom_density(alpha=0.3) + 
  scale_x_continuous(limits=c(0.2,0.6),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,20),expand=c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())

gcDepthPlot <- ggplot(data,aes(x=Depth,y=GC,
                               colour=Scaffold=="Scaffold_4")) + 
  geom_point(alpha=0.01) +
  scale_y_continuous(limits=c(0.2,0.6),expand=c(0,0)) +
  scale_x_continuous(limits=c(5,70),expand=c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        axis.ticks=element_blank())

depthDens <- ggplot(data,aes(x=Depth,group=Scaffold,
                             colour=Scaffold=="Scaffold_4",
                             fill=Scaffold=="Scaffold_4")) +
  geom_density(alpha=0.3) +
  scale_x_continuous(limits=c(5,70),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.2),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        legend.position="none")

fullPlot <- depthDens + plot_spacer() + 
  gcDepthPlot + gcDens + 
  plot_layout(ncol=2,widths=c(3,1),nrow=2,heights=c(1,3))

# Lastly, plot simple histograms of depth
OmniC_Depth_Hist <- ggplot(data,aes(x=Depth)) + geom_histogram() + facet_wrap(~Scaffold) +
  scale_x_continuous(limits=c(5,70))

png("OmniC_Depth_Hist.png",res=300,width=8,height=8,units="in")
plot(OmniC_Depth_Hist)
dev.off()

png("OmniC_scaffoldDepths.png",res=300,width=10,height=10,units="in")
plot(depthPlot)
dev.off()

png("OmniC_scaffoldDepths_GCnorm.png",res=300,width=10,height=10,units="in")
plot(depthPlot_GCnorm)
dev.off()

png("OmniC_gcVdepth.png",res=300,height=10,width=10,units="in")
plot(fullPlot)
dev.off()


#### OMNIC, MAPQ60 ALIGNMENT

# Read data
data <- read.table("OmniC_scaffoldGCandDepth_Q60.txt",header=F,stringsAsFactors=F)

# Add column names, make GC bins, clean up scaffold ID
colnames(data) <- c("Scaffold","WindowStart","WindowEnd","AT","GC","Depth")
data$GCbin <- cut(data$GC,seq(from=0,to=1,by=0.02))
Scaf_Order <- paste0("Scaffold_",c(1:17))
data$Scaffold <- factor(data$Scaffold,
                        levels=Scaf_Order,
                        labels=Scaf_Order)


# Create summary table for lookup
summaryGC <- data %>% 
  group_by(GCbin) %>%
  summarise(meanDepth=mean(Depth),sdDepth=sd(Depth))

# Empty vector for storage
zD <- c()

# For all rows
for(i in 1:nrow(data)){
  
  # Get the mean and sd of depth in that GC bin
  meanGC <- summaryGC$meanDepth[summaryGC$GCbin == data[i,7]]
  sdGC <- summaryGC$sdDepth[summaryGC$GCbin == data[i,7]]
  
  # Compute Z-score of depth
  zD[i] <- (data[i,6] - meanGC) / sdGC
  
}

data$depthZscore <- zD

# Plot Z-scores of GC normalisation
# This may help in future identifying duplications
depthPlot_GCnorm <- data %>% 
  ggplot(aes(y=depthZscore,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(-3,3)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="depth")


# Plot depth, this is all that's necessary for getting the putative sex chromosome
depthPlot <- data %>% 
  ggplot(aes(y=Depth,x=WindowStart)) + 
  geom_point(alpha=0.01) + 
  facet_wrap(~Scaffold,nrow=4,ncol=5,scales="free") +
  scale_y_continuous(limits=c(0,155)) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  labs(y="depth")

# Plot relationship betwen GC and Depth overall
gcDens <- ggplot(data,aes(x=GC,group=Scaffold,
                          colour=Scaffold=="Scaffold_4",
                          fill=Scaffold=="Scaffold_4")) + 
  geom_density(alpha=0.3) + 
  scale_x_continuous(limits=c(0.2,0.6),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,20),expand=c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())

gcDepthPlot <- ggplot(data,aes(x=Depth,y=GC,
                               colour=Scaffold=="Scaffold_4")) + 
  geom_point(alpha=0.01) +
  scale_y_continuous(limits=c(0.2,0.6),expand=c(0,0)) +
  scale_x_continuous(limits=c(5,70),expand=c(0,0)) +
  theme_bw() +
  theme(legend.position="none",
        axis.ticks=element_blank())

depthDens <- ggplot(data,aes(x=Depth,group=Scaffold,
                             colour=Scaffold=="Scaffold_4",
                             fill=Scaffold=="Scaffold_4")) +
  geom_density(alpha=0.3) +
  scale_x_continuous(limits=c(5,70),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.2),expand=c(0,0)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        legend.position="none")

fullPlot <- depthDens + plot_spacer() + 
  gcDepthPlot + gcDens + 
  plot_layout(ncol=2,widths=c(3,1),nrow=2,heights=c(1,3))

# Lastly, plot simple histograms of depth

OmniC_Depth_Hist <- ggplot(data,aes(x=Depth)) + geom_histogram() + facet_wrap(~Scaffold) +
  scale_x_continuous(limits=c(5,70))

png("OmniC_Depth_Hist_Q60.png",res=300,width=8,height=8,units="in")
plot(OmniC_Depth_Hist)
dev.off()

png("OmniC_scaffoldDepths_Q60.png",res=300,width=10,height=10,units="in")
plot(depthPlot)
dev.off()

png("OmniC_scaffoldDepths_GCnorm_Q60.png",res=300,width=10,height=10,units="in")
plot(depthPlot_GCnorm)
dev.off()

png("OmniC_gcVdepth_Q60.png",res=300,height=10,width=10,units="in")
plot(fullPlot)
dev.off()

