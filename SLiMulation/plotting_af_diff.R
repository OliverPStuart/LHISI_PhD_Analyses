
# Processing slimulation outputs for allele frequency change dynamics

# Environment

HOME_DIR="/Volumes/Alter/LHISI"
REF_DIR=paste0(HOME_DIR,"/References")

# Libraries

library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
options(dplyr.summarise.inform = FALSE)

# First get actual data
# Now bring in actual data

# Environment

WORKING_DIR=paste0(HOME_DIR,"/Analyses/DeleteriousMutations")
setwd(WORKING_DIR)

# Making data.frame

# Function for complement

`%ni%` <- Negate(`%in%`)

# Read data in, rename for ease of use

Captive <- read.table("../AlleleFrequencies/lhisi.mafs.gz",
                      header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
Hybrid <- read.table("../AlleleFrequencies/lhip.mafs.gz",
                     header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
Wild <- read.table("../AlleleFrequencies/Wild.mafs.gz",header=T,
                   stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
colnames(Captive) <- c("Scaffold","Position","MAF_Captive","N_Captive")
colnames(Hybrid) <- c("Scaffold","Position","MAF_Hybrid","N_Hybrid")
colnames(Wild) <- c("Scaffold","Position","MAF_Wild","N_Wild")

# Read in effect data

effects <- read.table("vars_classified.txt",
                      header=T,stringsAsFactors=F,sep="\t")[,c(1,2,3,4,5,6,7,8,10)]

# Get site information

scafs <- c(1:17)[-4]
sites <- data.frame(Scaffold=character(),
                    Position=character(),
                    Ref=character(),
                    Alt=character())
for(i in 1:length(scafs)){
  tmp <- eval(parse(text=paste0("read.table('../AlleleFrequencies/scaffold",
                                scafs[i],
                                "_sites',header=F,stringsAsFactors=F)")))
  sites <- rbind(sites,tmp)
}
rm(tmp) ; colnames(sites) <- c("Scaffold","Position","Ref","Alt")

# Give all sites their alleles

effects <- merge(effects,sites)
rm(sites)

# For each population, score every site by that population's frequency for the allele

classify_maf <- function(x){
  if(x <= 0.05){return("absent")} else
    if (x & x < 0.95 ){return("segregating")} else
      if (x >= 0.95 ){return("fixed")}
}

Wild$fate_Wild <- sapply(Wild$MAF_Wild,classify_maf)
Captive$fate_Captive <- sapply(Captive$MAF_Captive,classify_maf)
Hybrid$fate_Hybrid <- sapply(Hybrid$MAF_Hybrid,classify_maf)

# Now get all of these into the effects dataset

effects <- merge(effects,Wild,all.x=T)
effects <- merge(effects,Captive,all.x=T)
effects <- merge(effects,Hybrid,all.x=T)

# Make an ID column for later

effects$ID <- paste0(effects$Scaffold,"_",effects$Position)

rm(Captive,Hybrid,Wild)

# Some refactorisation to order specific variables

effects$Variant_Effect <- factor(effects$Variant_Effect,
                                 levels=c("MODIFIER","LOW","MODERATE","HIGH"),
                                 labels=c("MODIFIER","LOW","MODERATE","HIGH"))

# Get wild alleles

Wild <- effects %>% filter(fate_Wild != "absent" & !is.na(fate_Wild))

# Calculate difference in frequency

Wild$Captive_change <- Wild$MAF_Captive - Wild$MAF_Wild
Wild$Hybrid_change <- Wild$MAF_Hybrid - Wild$MAF_Wild
Wild <- Wild %>% drop_na(Captive_change,Hybrid_change)

Wild <- Wild %>% select(Variant_Effect,fate_Captive,fate_Hybrid,
                        Hybrid_change,Captive_change)

# Now recode factors and rename variables to agree with the simulated dataset
colnames(Wild) <- c("effect","fate_lhisi","fate_lhip","lhip_change","lhisi_change")

Wild$effect <- recode(Wild$effect,
                      MODIFIER = "neutral",
                      LOW = "weak",
                      MODERATE = "mod",
                      HIGH = "strong")

# Now get the simulated data

# Environment 

WORKING_DIR=paste0(HOME_DIR,"/Analyses/SLiMulation/10112022_1456")
setwd(WORKING_DIR)

# Change into run specific directory
#setwd("")

# Define number of replicates to extract data from and empty data.frame to hold results
reps=1:50
full_data <- data.frame()

# Loop over simulation outputs to get a data.frame with all outputs and replicate tags
for(rep in 1:length(reps)){
  
  print(reps[rep])
  
  # If all files do not exist, then go to next iteration
  if(file.exists(paste0("MUTS_",reps[rep],".txt")) +
     file.exists(paste0("WILD_MUTS_",reps[rep],".txt")) + 
     file.exists(paste0("LHISI_MUTS_",reps[rep],".txt")) + 
     file.exists(paste0("LHIP_MUTS_",reps[rep],".txt")) != 4){
       next()
  }
  # Files might not exist if simulation ended early due to extinction

  # Read in mutation table
  mutations <- read.table(paste0("MUTS_",reps[rep],".txt"),header=T)
  
  # Read in all mutation lists and format into frequencies
  wild <- as.numeric(readLines(paste0("WILD_MUTS_",reps[rep],".txt")))
  wild <- data.frame(mutation=names(xtabs(~wild)),
                     freq_wild=as.vector(xtabs(~wild)/4))
  lhisi <- as.numeric(readLines(paste0("LHISI_MUTS_",reps[rep],".txt")))
  lhisi <- data.frame(mutation=names(xtabs(~lhisi)),
                      freq_lhisi=as.vector(xtabs(~lhisi)/18))
  lhip <- as.numeric(readLines(paste0("LHIP_MUTS_",reps[rep],".txt")))
  lhip <- data.frame(mutation=names(xtabs(~lhip)),
                     freq_lhip=as.vector(xtabs(~lhip)/16))
  
  # Write fate columns
  classify_maf <- function(x){
    if(x <= 0.05){return("absent")} else
      if (x > 0.05 & x < 0.95 ){return("segregating")} else
        if (x >= 0.95 ){return("fixed")}
  }
  
  wild$fate_wild <- sapply(wild$freq_wild,classify_maf)
  lhisi$fate_lhisi <- sapply(lhisi$freq_lhisi,classify_maf)
  lhip$fate_lhip <- sapply(lhip$freq_lhip,classify_maf)
  
  # Merge all data together
  mutations <- merge(mutations,wild,all.x=T)
  mutations <- merge(mutations,lhisi,all.x=T)
  mutations <- merge(mutations,lhip,all.x=T)
  
  # Remove all mutations not present in sample
  missing_mut <- function(x) {
    
    # accessing elements from first column
    if(is.na(x[3]) & is.na(x[5]) & is.na(x[7])){
      return(F)
    } else {
      return(T)
    }
  }
  
  mutations <- mutations[apply(X=mutations,FUN=missing_mut,MARGIN=1),]
  
  # Remove anything fixed in all samples, these are not interesting
  for(i in 1:nrow(mutations)){
    
    if(i == 1){vec <- c()}
    
    tmp <- as.vector(mutations[i,c(3,5,7)])
    if(sum(is.na(tmp)) == 0 & sum(tmp) == 3) { 
      vec[i] <- F } else { vec[i] <- T 
      }
    
  }
  mutations <- mutations[vec,]
  
  # Now code to generate categories of mutations
  # classify_var <- function(x){
  #   if(x <= -0.1){return("v_strong")} else
  #     if(x > -0.1 & x <= -0.01 ){return("strong")} else
  #       if(x > -0.01 & x <= -0.001 ){return("mod")} else {
  #         if(x > -0.001 & x < 0 ){return("weak")} else {
  #           if(x == 0){return("neutral")}
  #         }
  #     }
  # }
  
  # Alternatively, collapse strong and v strong
  classify_var <- function(x){
    if(x <= -0.01 ){return("strong")} else
      if(x > -0.01 & x <= -0.001 ){return("mod")} else {
        if(x > -0.001 & x < 0 ){return("weak")} else {
          if(x == 0){return("neutral")}
        }
      }
  }
  
  
  mutations$effect <- sapply(mutations$s,classify_var)
  
  # Turn into factor for plotting order
  mutations$effect <- factor(mutations$effect,
                             levels=c("neutral","weak","mod","strong"),
                             labels=c("neutral","weak","mod","strong"))
  
  # Now calculate metrics
  mutations$lhisi_change <- mutations$freq_lhisi - mutations$freq_wild
  mutations$lhip_change <- mutations$freq_lhip - mutations$freq_wild
  
  # Mean +- sd plot
  # mutations %>%
  #   filter(fate_lhisi == "segregating",effect != "v_strong") %>%
  #   group_by(effect) %>%
  #   filter(!is.na(lhisi_change)) %>%
  #   dplyr::summarise(mean=mean(lhisi_change),
  #             min=mean(lhisi_change)-sd(lhisi_change),
  #             max=mean(lhisi_change)+sd(lhisi_change),
  #             effect = first(effect)) %>%
  #   ggplot(aes(x=effect,y=mean,ymin=min,ymax=max)) +
  #   geom_errorbar(width=0) +
  #   geom_point(size=3,pch=21, fill = "cornflowerblue")
  
  # And pull to add to combined data.frame
  
  mutations <- mutations %>% 
    select(effect,lhisi_change,lhip_change,fate_lhisi,fate_lhip) %>%
    mutate(replicate=reps[rep])
  
  full_data <- rbind(full_data,
                     mutations)
  
}

rm(mutations, lhisi, lhip, wild, tmp)

# Get colour scheme

# Bring in 2-dimensional colour scheme and make legend figures

new_colours <- read.delim(paste0(REF_DIR,"/two_pop_2d_colour_scale.txt"),
                          header=F,stringsAsFactors = F)
colnames(new_colours) <- c("effect","LHISI","LHIP","Wild")

new_colours$effect <- factor(new_colours$effect,levels=c("neutral","weak","mod","strong"),
                             labels=c("Non-coding","Weak","Moderate","Strong"))


# Now plot
# Plot all replicates together and then plot the actual data on top of it

full_data$effect <- factor(full_data$effect,levels=c("neutral","weak","mod","strong"),
                           labels=c("Non-coding","Weak","Moderate","Strong"))

Wild$effect <- factor(Wild$effect,levels=c("neutral","weak","mod","strong"),
                      labels=c("Non-coding","Weak","Moderate","Strong"))

p1 <- full_data %>%
  filter(fate_lhisi == "segregating") %>%
  group_by(effect) %>%
  filter(!is.na(lhisi_change)) %>%
  ggplot(aes(x=lhisi_change,group=replicate)) +
  stat_ecdf(alpha=0.2) + facet_grid(effect~.) + 
  scale_x_continuous(limits=c(-1,1)) + 
  scale_y_continuous(limits=c(-0.05,1.05),breaks=c(0,0.5,1),labels=c(0,0.5,1),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        strip.text.y=element_text(angle=0),
        axis.text.y=element_blank(),
        legend.position="none") + 
  stat_ecdf(data=Wild[Wild$fate_lhisi == "segregating",],
            inherit.aes=F,aes(x=lhisi_change),colour="black",size=2.2) +
  stat_ecdf(data=Wild[Wild$fate_lhisi == "segregating",],
            inherit.aes=F,aes(x=lhisi_change,colour=effect),size=1.5) + 
  labs(x="Allele frequency change in captivity",y="Cumulative density") + 
  scale_colour_manual(values=new_colours$LHIP)

png("../lhisi_change_sim_100100_extra_generations_in_captivity.png",
    res=300,width=6,height=6,units='in')
plot(p1)
dev.off()

full_data %>%
  filter(fate_lhip == "segregating") %>%
  group_by(effect) %>%
  filter(!is.na(lhip_change)) %>%
  ggplot(aes(x=lhip_change,group=replicate)) +
  stat_ecdf(alpha=0.2) + facet_grid(effect~.) + 
  scale_x_continuous(limits=c(-1,1)) + 
  scale_y_continuous(limits=c(-0.05,1.05),breaks=c(0,0.5,1),labels=c(0,0.5,1),expand=c(0,0)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        strip.text.y=element_text(angle=0),
        axis.text.y=element_blank(),
        legend.position="none") + 
  stat_ecdf(data=Wild[Wild$fate_lhip == "segregating",],
            inherit.aes=F,aes(x=lhip_change),colour="black",size=2.2) +
  stat_ecdf(data=Wild[Wild$fate_lhip == "segregating",],
            inherit.aes=F,aes(x=lhip_change,colour=effect),size=1.5) + 
  labs(x="Allele frequency change in captivity",y="Cumulative density") + 
  scale_colour_manual(values=new_colours$LHISI)


