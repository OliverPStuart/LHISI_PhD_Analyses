### File for plotting depth across D. australis scaffolds

HOME_DIR="/Volumes/Alter/LHISI"
WORKING_DIR=paste0(HOME_DIR,"/Analyses/HuntingForSexChroms")
setwd(WORKING_DIR)

library(dplyr)
library(ggplot2)
library(patchwork)
library(svglite)

### Script

# Read in scaffold information
scaffolds <- read.table(paste0(HOME_DIR,"/References/major_scaffold_regions.bed"))[,c(1,3)]
scaffoldsV3 <- scaffolds$V3 + 1
colnames(scaffolds) <- c("Scaffold","Length")

# Get start and end in global coordinates
scaffolds$CumEnd <- cumsum(as.numeric(scaffolds$Length))
scaffolds$CumStart <- c(1,lag(scaffolds$CumEnd)[-1]+1)

# Read data
data <- read.table("PacBio_scaffoldGCandDepth.txt",header=F,stringsAsFactors=F)

#### For now, subsample to make it easier to plot and test
data <- data[seq(from=1,nrow(data),by=50),]

# Add column names, make midpoint variable
colnames(data) <- c("Scaffold","WindowStart","WindowEnd","AT","GC","Depth")

# Give every position a cumulative position
data <- merge(scaffolds,data) %>% mutate(CumPosition=WindowStart+CumStart-1)

# Get positions of axis labels
axis_set <- scaffolds %>% 
  group_by(Scaffold) %>% 
  summarize(center = (CumStart+CumEnd)/2)
axis_set$Scaffold <- gsub("Scaffold_","",axis_set$Scaffold)

# Set y-axis minimum and maximum
ymin=0
ymax <- data %>% 
  mutate(ylim = mean(Depth) + 200) %>% 
  summarise(ylim=first(ylim)) %>%
  pull(ylim)

# Data.frame to make panel background rectangles
rect_breaks <- seq(from=ymin,to=ymax,by=(ymax-ymin)/5)
rect <- data.frame(ymin=rect_breaks[c(2,4)],
                   ymax=rect_breaks[c(3,5)],
                   xmax=Inf,xmin=-Inf)

# Make y_labels from this, max and min slightly nudged inward 
y_labels <- rect_breaks
y_labels[1] <- rect_breaks[1] + (rect_breaks[6]/100)
y_labels[6] <- rect_breaks[6] - (rect_breaks[6]/100)

p <- ggplot(data,aes(x=CumPosition,y=Depth,colour=as.factor(Scaffold))) + 
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray95",fill="gray95") +
  geom_point(alpha=0.2) + 
  scale_x_continuous(label = axis_set$Scaffold, breaks = axis_set$center,
                     limits=c(1-(max(scaffolds$CumEnd)*0.01),
                              max(scaffolds$CumEnd)+(max(scaffolds$CumEnd)*0.01)),
                     expand=c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, ymax)) +
  scale_color_manual(values = rep(c("palegreen4", "palegreen3"),
                                  unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Scaffold",
       y = "Depth") + 
  theme( 
    legend.position = "none",
    panel.grid=element_blank(),
    panel.background = element_blank(),
    axis.ticks=element_blank(),
    panel.border = element_rect(colour="black",fill=NA,size=1),
    axis.text=element_text(size=7), axis.title=element_text(size=10)
  )


svg("continuous_depth.svg",height=4,width=8)
plot(p)
dev.off()

png("continuous_depth.png",height=4,width=8,units='in',res=300)
plot(p)
dev.off()

