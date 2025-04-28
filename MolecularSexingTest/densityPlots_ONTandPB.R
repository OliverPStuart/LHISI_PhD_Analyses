library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# Read in table

ont <- read.table("ONT_Scaffold.regions.bed.gz",header=F,stringsAsFactors = F,sep="\t")
pb <- read.table("../HuntingForSexChroms/PacBio_scaffoldGCandDepth.txt",header=F,stringsAsFactors = F,sep="\t")

ont$norm <- ont$V4/mean(ont$V4)
ggplot(ont,aes(x=V2,y=norm)) + geom_point(alpha=0.01) + facet_wrap(~V1,scales = "free_x") + scale_y_continuous(limits=c(0,2))

pb$norm <- pb$V6/mean(pb$V6[pb$V1 != "Scaffold_4"])
ggplot(pb,aes(x=V2,y=norm)) + geom_point(alpha=0.01) + facet_wrap(~V1,scales = "free_x") + scale_y_continuous(limits=c(0,2))

combined <- data.frame(Scaffold=ont$V1,
                       Position=(ont$V2+ont$V3)/2,
                       Male=pb$norm,
                       Female=ont$norm)

colours.vec=c("mediumpurple1","springgreen4")

dots <- combined %>% filter(Scaffold=="Scaffold_4") %>%
  gather(key=Sex,value=Normalised_Depth,-Scaffold,-Position) %>%
  ggplot(aes(x=Position,y=Normalised_Depth,colour=Sex)) + 
  geom_point(alpha=0.15) + 
  scale_y_continuous(limits=c(0,1.75),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,326569942),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),axis.title.x=element_blank(),
        legend.position="none",
        axis.text.x=element_text(hjust=-0.1)) + 
  scale_colour_manual(values=colours.vec)


denss <- combined %>% filter(Scaffold=="Scaffold_4") %>%
  gather(key=Sex,value=Normalised_Depth,-Scaffold,-Position) %>%
  ggplot(aes(x=Normalised_Depth,fill=Sex)) + 
  geom_density() + 
  scale_x_continuous(limits=c(0,1.75),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,5),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),axis.title=element_blank(),
        axis.text=element_blank()) + 
  coord_flip() + 
  scale_fill_manual(values=colours.vec)

png("Depths.png",res=300,width=10,height=8,units="in")
plot(dots+denss+plot_layout(widths=c(1,0.5)))
dev.off()
  
