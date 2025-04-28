
library(ggplot2)
library(patchwork)
lib

# Read in data
d_o_mn <- read.table("sexChromContigs_OmniC.regions.bed.gz",header=F,stringsAsFactors=F)
d_p_mn <- read.table("sexChromContigs_PacBio.regions.bed.gz",header=F,stringsAsFactors=F)


# Plot both with and without depth filters

plot_function <- function(x) {
  
  tempmean <- mean(x[,5])
  tempsd <- sd(x[,5])
  p <- (ggplot(x,aes(x=V5)) + geom_histogram() + 
          scale_x_continuous(limits=c(tempmean-tempsd,tempmean+tempsd),trans="log10") + 
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank()) + 
    labs(x="depth of coverage") + ggtitle("expected depth range")) +
  (ggplot(x,aes(x=V5)) + geom_histogram() + scale_x_continuous(trans="log10") + 
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank()) + 
    labs(x="depth of coverage") + ggtitle("overall depth range")) +
    (ggplot(x,aes(x=V3-V2,y=V5)) + geom_point() + scale_y_continuous(trans="log10") +
       scale_x_continuous(trans="log10") +
       labs(x="contig length",y="contig depth"))
  
  plot(p)
  
}

png("OmniC_SexChromDepths_AlignmentToScaffold.png",res=300,width=10,height=5,units='in')
plot_function(d_o_mn)
dev.off()
png("PacBio_SexChromDepths_AlignmentToScaffold.png",res=300,width=10,height=5,units='in')
plot_function(d_p_mn)
dev.off()

# One final plot, correlation between OmniC and PacBio read depth
png("ReadDepthCorrelation.png",res=300,width=5,height=5,units='in')
ggplot(data.frame(OmniC=d_o_mn$V5,PacBio=d_p_mn$V5),
       aes(x=OmniC,y=PacBio)) + geom_point() +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10")
dev.off()
