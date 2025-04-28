
### Comparing heterozygosity per scaffold with different call sets

library(ggplot2)
library(tidyverse)
library(lemon)
### Read in data

lengths <- read.table("scaffoldLengths.txt",header=F,stringsAsFactors = F)
colnames(lengths) <- c("Scaffold","Length")
lengths$Scaffold <- as.numeric(unlist(lapply(str_split(lengths$Scaffold,pattern = "_"),`[[`,2)))

OmniC <- read.table("OmniC_Het_Table.txt",header=T,stringsAsFactors=F)
PacBio <- read.table("PacBio_Het_Table.txt",header=T,stringsAsFactors=F)
OmniC$Scaffold <- as.numeric(unlist(lapply(str_split(OmniC$Scaffold,pattern = "_"),`[[`,2)))
PacBio$Scaffold <- as.numeric(unlist(lapply(str_split(PacBio$Scaffold,pattern = "_"),`[[`,2)))

# Re-order OmniC... should really go back and fix this so it's not an issue
OmniC <- OmniC[order(OmniC$Scaffold,decreasing=F),]

# Combine
comb <- cbind(OmniC,PacBio[,-1]) %>% gather(key="Dataset",-Scaffold,value="Variants")
comb <- merge(comb,lengths,by="Scaffold")
comb$Ho <- comb$Variants/comb$Length

scafOrder <- c(1:3,5:17)
comb$Scaffold <- factor(comb$Scaffold,labels=scafOrder,levels=scafOrder)

# Plot
overlay <- comb %>% 
  group_by(Scaffold) %>% 
  summarise(average=mean(Ho),
            max=max(Ho),
            min=min(Ho)) 

p <- ggplot() + 
  geom_point(overlay,mapping=aes(x=Scaffold,y=average),
             alpha=0.3) +
  geom_linerange(overlay,mapping=aes(x=Scaffold,ymin=min,ymax=max),
                alpha=0.3) +
  geom_point(comb,mapping=aes(x=Scaffold,y=Ho,colour=Dataset),
             size=3, inherit.aes=F,
             position=position_dodge2(width=0.8,reverse = F)) +
  scale_y_continuous(trans="log10") +
  theme(axis.ticks=element_blank()) +
  theme_bw() +
  labs(y="Average heterozygosity")

png("comparingDatasetHo.png",res=300,width=8,height=8,units='in')
plot(p)
dev.off()
