# Analysing coverage across the minor scaffolds with male/female data

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# Read in table

ont <- read.table("ONT_minorScaffold.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
pb <- read.table("PacBio_minorScaffold.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
both <- merge(ont,pb,by=c("V1","V2","V3"))
colnames(both) <- c("Scaffold","Start","End","F_Depth","M_Depth")

# Select places where F_Depth == 0, over the whole contig
# Select for minimum mean depth in male to remove anything that looks odd

F_Missing <- both %>% group_by(Scaffold) %>%
  summarise(Fe = mean(F_Depth),Ma=mean(M_Depth)) %>%
  ungroup %>% filter(Fe == 0 & Ma > 10) %>% select(Scaffold)

# Now for all of these, we make a plot of M and F depth
# First truncate upper tail of distributions for visualisation

both$M_Depth[both$M_Depth > 250] <- 250
both$F_Depth[both$F_Depth > 250] <- 250
both$Mid <- (both$Start + both$End)/2

for(i in 1:nrow(F_Missing)){
  
  scaf <- as.character(F_Missing[i,])
  temp <- gather(both[both$Scaffold == scaf,],key="Sex",value="Depth",-Scaffold,-Start,-End,-Mid)
  
  line_plot <- ggplot(temp,aes(x=Mid,y=Depth)) + geom_line() + 
    facet_wrap(~Sex,nrow=2) + scale_y_continuous(limits=c(0,300))
  hist_plot <- ggplot(temp,aes(x=Depth)) + geom_histogram() + 
    facet_wrap(~Sex,nrow=2) + scale_x_continuous(limits=c(0,300))
  
  png(paste0("minorScaffoldPlots/",scaf,".png"),res=300,width=7,height=7,units='in')
  plot(line_plot + hist_plot)
  dev.off()
  
}

# And from these plots, we can identify by eye the contigs to focus on

# Quickly check that there aren't any reciprocal cases
# i.e. where there is male missing but high female depth

M_Missing <- both %>% group_by(Scaffold) %>%
  summarise(Fe = mean(F_Depth),Ma=mean(M_Depth)) %>%
  ungroup %>% filter(Fe > 0 & Ma == 0) %>% select(Scaffold)

M_Missing_Depths <- both[both$Scaffold %in% M_Missing$Scaffold,]
View(M_Missing_Depths)

# Mostly very low, so we can safely ignore these, and focus on F_Missing contigs with higher depth