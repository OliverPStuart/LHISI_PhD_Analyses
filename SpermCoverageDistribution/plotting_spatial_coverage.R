
### Plotting sperm coverage along scaffolds
### Here, we plot the spatial distribution

library(ggplot2)
library(dplyr)

# Read data, remove minor scaffolds

data <- read.table("spatial_coverage_sperm.txt",header=T,stringsAsFactors=F,sep="\t")
major <- paste0("Scaffold_",c(1:17))
data_major <- data[data$Scaffold %in% major,]

data_major$Scaffold <- factor(data_major$Scaffold,levels=paste0("Scaffold_",1:17))

samples <- unique(data$Sample)

for(i in 1:length(samples)){
  
  sample_i <- samples[i]
  
  p <- data_major %>% filter(Sample==sample_i) %>%
    ggplot(aes(x=Start,y=Coverage)) + geom_point(alpha=0.1) +
    facet_wrap(~Scaffold,scales="free_x") + 
    scale_y_continuous(trans="log10") + 
    geom_smooth() + ggtitle(sample_i) +
    theme(axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank())
  
  png(paste0(sample_i,"_spatial_coverage.png"),res=300,width=7,height=7,units='in')
  plot(p)
  dev.off()
  
}
