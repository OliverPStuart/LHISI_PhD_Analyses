
library(dplyr)
library(ggplot2)
library(patchwork)

# Plotting the locations of the Timema X chromosome gene blast hits

# First, bring in PacBio and OmniC depth data
pacbio_depth_all <- read.table("PacBio_scaffoldGCandDepth.txt",header=F,stringsAsFactors=F)
omnic_depth_all <- read.table("OmniC_scaffoldGCandDepth_Q60.txt",header=F,stringsAsFactors=F)
colnames(pacbio_depth_all) <- c("Scaffold","Start","End","AT","GC","Depth")
colnames(omnic_depth_all) <- c("Scaffold","Start","End","AT","GC","Depth")
scafOrder <- paste0("Scaffold_",c(1:17))
pacbio_depth_all$Scaffold <- factor(pacbio_depth_all$Scaffold,labels=scafOrder,levels=scafOrder)
omnic_depth_all$Scaffold <- factor(omnic_depth_all$Scaffold,labels=scafOrder,levels=scafOrder)

# Now bring in table of locations
hits <- read.table("Timema_X_Best_Hits.txt",header=T,stringsAsFactors = F)
hits$Scaffold <- factor(hits$Scaffold,labels=scafOrder,levels=scafOrder)

# Now plot the locations of the hits along all scaffolds

all_scafs_pacbio_hits <- ggplot() + 
  geom_point(data=pacbio_depth_all,aes(x=Start,y=Depth),alpha=0.01) + 
  geom_vline(data=hits,aes(xintercept=Start),colour="springgreen3", lty="dashed", lwd=0.4) +
  facet_wrap(~Scaffold,scales="free_x") +
  scale_y_continuous(limits=c(0,300)) +
  ggtitle("PacBio, all scaffolds")

all_scafs_omnic_hits <- ggplot() + 
  geom_point(data=omnic_depth_all,aes(x=Start,y=Depth),alpha=0.01) + 
  geom_vline(data=hits,aes(xintercept=Start),colour="springgreen3", lty="dashed", lwd=0.4) +
  facet_wrap(~Scaffold,scales="free_x") +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("OmniC (filtered), all scaffolds")

png("Timema_X_Hits_All_Scaffolds.png",res=300,width=20,height=10,units='in')
plot(all_scafs_pacbio_hits + all_scafs_omnic_hits)
dev.off()


# Now bring in the hits with their depth information for pacbio and omnic
hits_4_omnic <- read.table("OmniCQ60_TimemaXtoScaffold.regions.bed.gz",
                           header=F,
                           stringsAsFactors=F) %>% filter(V1 == "Scaffold_4")
hits_4_pacbio <- read.table("PacBio_TimemaXtoScaffold.regions.bed.gz",
                           header=F,
                           stringsAsFactors=F) %>% filter(V1 == "Scaffold_4")
colnames(hits_4_omnic) <- c("Scaffold","Start","End","Depth")
colnames(hits_4_pacbio) <- c("Scaffold","Start","End","Depth")
hits_4_omnic$Scaffold <- factor(hits_4_omnic$Scaffold,labels=scafOrder,levels=scafOrder)
hits_4_pacbio$Scaffold <- factor(hits_4_pacbio$Scaffold,labels=scafOrder,levels=scafOrder)
pacbio_depth_4 <- pacbio_depth_all %>% filter(Scaffold=="Scaffold_4")
omnic_depth_4 <- omnic_depth_all %>% filter(Scaffold=="Scaffold_4")


# Now plot their locations on the chromosome as well as their locations in the depth histogram
# Pacbio
# Locations
scaf_4_pacbio_hits <- ggplot() + 
  geom_point(data=pacbio_depth_4,aes(x=Start,y=Depth),alpha=0.1) +
  geom_vline(data=hits_4_pacbio,aes(xintercept=Start),colour="springgreen3", lty="dashed", lwd=0.4) +
  scale_y_continuous(limits=c(0,150),expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank())
# Histogram
scaf_4_pacbio_hist <- ggplot() +
  geom_histogram(data=pacbio_depth_4,aes(x=Depth)) +
  geom_vline(data=hits_4_pacbio,aes(xintercept=Depth),colour="springgreen3", lty="dashed", lwd=0.4) +
  scale_x_continuous(limits=c(0,150),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,6000),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()) +
  coord_flip()

png("Timema_X_Best_Hits_PacBio_Scaf_4.png",res=300,width=10,height=8, units='in')
plot(scaf_4_pacbio_hits + scaf_4_pacbio_hist + plot_layout(widths=c(7,3)))
dev.off()

# And OmniC
# Locations
scaf_4_omnic_hits <- ggplot() + 
  geom_point(data=omnic_depth_4,aes(x=Start,y=Depth),alpha=0.1) +
  geom_vline(data=hits_4_omnic,aes(xintercept=Start),colour="springgreen3", lty="dashed", lwd=0.4) +
  scale_y_continuous(limits=c(0,50),expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank())
# Histogram
scaf_4_omnic_hist <- ggplot() +
  geom_histogram(data=omnic_depth_4,aes(x=Depth)) +
  geom_vline(data=hits_4_omnic,aes(xintercept=Depth),colour="springgreen3", lty="dashed", lwd=0.4) +
  scale_x_continuous(limits=c(0,50),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,6500),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()) +
  coord_flip()

png("Timema_X_Best_Hits_OmniC_Scaf_4.png",res=300,width=10,height=8, units='in')
plot(scaf_4_omnic_hits + scaf_4_omnic_hist + plot_layout(widths=c(7,3)))
dev.off()

# Now, just for omnic, do the same thing, but using only the depths in 10000 windows
# i.e. the depth of the X gene hit is just the depth of the window it falls into

hits_4 <- hits %>% filter(Scaffold=="Scaffold_4")

# Reorder columns since some hits will be on -ve strand
for(i in 1:nrow(hits_4)){
  
  if(hits_4$Start[i] > hits_4$End[i]){
    Start=min(hits_4$Start[i],hits_4$End[i])
    End=max(hits_4$Start[i],hits_4$End[i])
    hits_4$Start[i] <- Start
    hits_4$End[i] <- End
  }
  
}

new_depths <- c()
# For every gene's coordinates
for(i in 1:nrow(hits_4)){
  
  # Find the bin with start coordinate that is less than 10'000 bp before the hit's start
  start_i <- omnic_depth_4[hits_4[i,3] - omnic_depth_4$Start < 10000 & 
                             hits_4[i,3] - omnic_depth_4$Start > 0,]
  
  # Find the bins with end coordinate that is less than 10'000 bp after the hit's end
  end_i <- omnic_depth_4[omnic_depth_4$End - hits_4[i,4] < 10000 & 
                           omnic_depth_4$End - hits_4[i,4] > 0,]
  
  # If they are the same, no problem
  if(sum(start_i == end_i) == ncol(start_i)){
    new_depths[i] <- start_i$Depth
  # If they are different and the hits falls into two bins, just take the mean
  # Thought about taking a weighted mean but that's too much effort for a simple task
  } else {
    new_depths[i] <- mean(start_i$Depth,end_i$Depth)
  }
  
  if(i == nrow(hits_4)){
    hits_4$New_Depth <- new_depths
  }
  
}

# Locations
scaf_4_omnic_hits <- ggplot() + 
  geom_point(data=omnic_depth_4,aes(x=Start,y=Depth),alpha=0.1) +
  geom_vline(data=hits_4,aes(xintercept=Start),colour="springgreen3", lty="dashed", lwd=0.4) +
  scale_y_continuous(limits=c(0,50),expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank())
# Histogram
scaf_4_omnic_hist <- ggplot() +
  geom_histogram(data=omnic_depth_4,aes(x=Depth)) +
  geom_vline(data=hits_4,aes(xintercept=New_Depth),colour="springgreen3", lty="dashed", lwd=0.4) +
  scale_x_continuous(limits=c(0,50),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,6500),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()) +
  coord_flip()

png("Timema_X_Best_Hits_OmniC_Scaf_4_BinDepths.png",res=300,width=10,height=8, units='in')
plot(scaf_4_omnic_hits + scaf_4_omnic_hist + plot_layout(widths=c(7,3)))
dev.off()

# Now a simple histogram of the bin depths
# Both the depths of the 10'000 bins they fall into, and the depths of the mapped regions themselves

bin_depths
p1 <- ggplot(data=hits_4,aes(x=New_Depth)) + 
  geom_histogram() +
  theme_bw() +
  scale_y_continuous(limits=c(0,35),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,55),expand=c(0,0)) +
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank()) +
  ggtitle("Average depth in 10'000 bins that Timema X hits fall inside")

p2 <- ggplot(data=hits_4_omnic,aes(x=Depth)) + 
  geom_histogram() +
  theme_bw() +
  scale_y_continuous(limits=c(0,35),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,55),expand=c(0,0)) +
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  ggtitle("Average depth in mapping regions of Timema X hits")

png("Depths_Timema_X_Hits.png",res=300,width=8,height=6,units='in')
plot(p1 + p2 + plot_layout(nrow=2))
dev.off()

