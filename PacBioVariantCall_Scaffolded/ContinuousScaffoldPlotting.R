
library(ggplot2)
library(dplyr)
library(patchwork)

# Variant counts
variants <- read.table("counts",header=F,stringsAsFactors = F)
colnames(variants) <- c("Scaffold","Start","End","Count")
variants$Scaffold <- as.numeric(gsub("Scaffold_","",variants$Scaffold))

# Windows
# Sum over windows to get per window masked
windows <- read.table("MaskedWindows.bed",header=F,stringsAsFactors = F)[,c(1,2,3,7)]
colnames(windows) <- c("Scaffold","Start","End","Mask")
windows <- windows %>% group_by(Scaffold,Start,End) %>%
  summarise(Start=min(Start),End=min(End),Mask=sum(Mask)) %>% 
  as.data.frame()
windows$Scaffold <- as.numeric(gsub("Scaffold_","",windows$Scaffold))

# Bring the two together
data <- merge(variants,windows,by=c("Scaffold","Start","End"))
data$IncludedLength <- 1e6 - data$Mask

# Order data.frame by scaffold order
data <- data[with(data, order(Scaffold, Start)),]

# Make bins for plotting
# Make continuous variable for plotting
data$Bin <- c(1:nrow(data))

# Get median and select corresponding bin number to get location of axis label
for(i in 1:length(unique(data$Scaffold))){
  
  if(i == 1){
    Locations <- c()
  }
  
  Locations <- c(Locations,
                 data %>% 
                   subset(Scaffold==unique(data$Scaffold)[i]) %>%
                   summarise(Bin=median(Bin)) %>% unlist)
  
  if(i == length(unique(data$Scaffold))){
    Locations <- data.frame(Locations=Locations,Scaffold=unique(data$Scaffold))
  }

}

# Some formatting
data$Scaffold <- factor(data$Scaffold,levels=c(1:17),labels=c(1:17))

# Vector of colours
colours <- rep(c("springgreen4","lightgreen"),9)[-18]

# Set y-axis minimum and maximum
ymin=0
ymax=100

# Data.frame to make panel background rectangles
rect_breaks <- seq(from=ymin,to=ymax,by=(ymax-ymin)/5)
rect <- data.frame(ymin=rect_breaks[c(2,4)],
                          ymax=rect_breaks[c(3,5)],
                          xmax=Inf,xmin=-Inf)

# Make y_labels from this, max and min slightly nudged inward 
y_labels <- rect_breaks
y_labels[1] <- rect_breaks[1] + (rect_breaks[6]/100)
y_labels[6] <- rect_breaks[6] - (rect_breaks[6]/100)

masked <- ggplot(data,aes(x=Bin,group=Scaffold,y=Mask/1e4,colour=Scaffold)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),panel.grid=element_blank(),
        legend.position="none") +
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray95",fill="gray95") +
  geom_line() +
  scale_x_continuous(expand=c(0,0),lim=c(1,nrow(data)),breaks=Locations$Locations,labels=Locations$Scaffold) +
  scale_y_continuous(breaks=y_labels,labels=rect_breaks) + 
  coord_cartesian(ylim=c(ymin,ymax),expand=F) + 
  scale_colour_manual(breaks=c(1:17),values=colours) +
  labs(y="Percent masked",x="Scaffold")


# Now formatting for heterozygosity
# We scale up the counts using the proportion of each window masked by repeats
# Then divide each bin to get a per kilobase value
#data$H_perKb <- data$Count * (1e6/(1e6-data$IncludedLength)) / 1000

# Alternative, just divide count by included length
data$H_perKb <- data$Count / data$IncludedLength #* 1000


# Set y-axis minimum and maximum
ymin=0
ymax=0.005

# Data.frame to make panel background rectangles
rect_breaks <- seq(from=ymin,to=ymax,by=(ymax-ymin)/5)
rect <- data.frame(ymin=rect_breaks[c(2,4)],
                   ymax=rect_breaks[c(3,5)],
                   xmax=Inf,xmin=-Inf)

# Make y_labels from this, max and min slightly nudged inward 
y_labels <- rect_breaks
y_labels[1] <- rect_breaks[1] + (rect_breaks[6]/100)
y_labels[6] <- rect_breaks[6] - (rect_breaks[6]/100)

het <- ggplot(data,aes(x=Bin,group=Scaffold,y=H_perKb,colour=Scaffold)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),panel.grid=element_blank(),
        legend.position="none") +
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray95",fill="gray95") +
  geom_line() +
  scale_x_continuous(expand=c(0,0),lim=c(1,nrow(data)),breaks=Locations$Locations,labels=Locations$Scaffold) +
  scale_y_continuous(breaks=y_labels,labels=rect_breaks) + 
  coord_cartesian(ylim=c(ymin,ymax),expand=F) + 
  scale_colour_manual(breaks=c(1:17),values=colours) +
  labs(y="Heterozygosity",x="Scaffold")

#(masked + theme(axis.title.x=element_blank(),axis.text.x=element_blank())) / het

png("het_per_kb_plot.png",res=300,width=7,height=4,units='in')
het
dev.off()

png("mask_per_kb_plot.png",res=300,width=7,height=4,units='in')
masked
dev.off()

# Now we are also going to bring in depth data and examine heterozygosity in these windows in relation to the mean depth
depths <- read.table("PacBio_Scaffold_1M.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
depths$V1 <- gsub("Scaffold_","",depths$V1)
depths$V1<- factor(depths$V1,levels=c(1:17),labels=c(1:17))
colnames(depths) <- c("Scaffold","Start","End","Depth")
data <- merge(data,depths,by=c("Scaffold","Start","End"))

ggplot(data,aes(x=Depth,y=H_perKb)) + geom_point() + 
  facet_wrap(~Scaffold) + 
  coord_cartesian(xlim=c(50,300),ylim=c(0,1.5))
# No relationship between long-read depth and heterozygosity
# Probably only decipherable with genotype data from multiple individuals





### Plotting for a black background
data$H <- data$H_perKb/1000

# Set y-axis minimum and maximum
ymin=0
ymax=5/1000

# Data.frame to make panel background rectangles
rect_breaks <- seq(from=ymin,to=ymax,by=(ymax-ymin)/5)
rect <- data.frame(ymin=rect_breaks[c(2,4)],
                   ymax=rect_breaks[c(3,5)],
                   xmax=Inf,xmin=-Inf)

# Make y_labels from this, max and min slightly nudged inward 
y_labels <- rect_breaks
y_labels[1] <- rect_breaks[1] + (rect_breaks[6]/100)
y_labels[6] <- rect_breaks[6] - (rect_breaks[6]/100)


het <- ggplot(data,aes(x=Bin,group=Scaffold,y=H,colour=Scaffold)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),panel.grid=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill="black"),
        panel.border = element_blank(),
        plot.background = element_rect(fill="black"),
        axis.title = element_text(colour="grey95"),
        axis.text=element_text(colour="grey95")) +
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray25",fill="gray25") +
  geom_line() +
  scale_x_continuous(expand=c(0,0),lim=c(1,nrow(data)),breaks=Locations$Locations,labels=Locations$Scaffold) +
  scale_y_continuous(breaks=y_labels,labels=c("0","1e-3","2e-3","3e-3","4e-3","5e-3")) + 
  coord_cartesian(ylim=c(ymin,ymax),expand=F) + 
  scale_colour_manual(breaks=c(1:17),values=colours) +
  labs(y="Heterozygosity",x="Scaffold")

png("het_blackbackground.png",res=300,width=7,height=3,units='in')
plot(het)
dev.off()
