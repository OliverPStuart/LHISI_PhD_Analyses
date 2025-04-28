# Processing slimulation outputs for half-sib heterozygosity

# Environment

HOME_DIR="/Volumes/Alter/LHISI"
REF_DIR=paste0(HOME_DIR,"/References")

# Libraries

library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
options(dplyr.summarise.inform = FALSE)

# Environment

WORKING_DIR=paste0(HOME_DIR,"/Analyses/SLiMulation")
setwd(WORKING_DIR)

# Read in data

het <- read.table("06012023_1131_half_sib_het_results.txt",
                  sep=",",stringsAsFactors=F,header=T)

# Make difference variable

het$diff <- abs(het$ind1 - het$ind2)

# Now get data for the two real individuals

real_hets <- read.table("../SingleSampleH/mean_heterozygosity_autosome.txt",
                        sep="\t",stringsAsFactors=F,header=T)


real_diffs <-real_hets %>% filter(Pop == "WILD") %>% 
  select(Sample,Depth,mean_het) %>%
  group_by(Depth) %>% 
  pivot_wider(names_from = Sample,values_from = mean_het) %>%
  ungroup %>%
  mutate(diff = abs(PAUL4 - VAN2)) %>% pull(diff)

# Plot this as histogram

p <- ggplot(het,aes(x=diff)) + 
  geom_histogram(fill="grey94",colour="black") + 
  geom_vline(xintercept=real_diffs,
             colour=c("springgreen3","cornflowerblue"),
             size=1.5,
             linetype="dashed") + 
  theme_bw() + 
#  scale_y_continuous(limits=c(0,90),expand=c(0,0)) +
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) + 
  labs(x="Difference in H\nbetween half-siblings")


png("half_sub_het_k300.png",res=300,width=6.5,height=6.5,units='in')
plot(p)
dev.off()


