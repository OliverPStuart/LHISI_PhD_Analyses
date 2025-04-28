
### Plotting sperm coverage across scaffolds
### Starting with just per-scaffold # reads

library(ggplot2)
library(dplyr)

# Read data

data_sperm <- read.table("coverage_stats_sperm.txt",header=T,stringsAsFactors=F,sep="\t")

# Start by splitting into major and minor
# We want to aggregate all minor scaffolds together

`%ni%` <- Negate(`%in%`)
major <- paste0("Scaffold_",c(1:17))
data_sperm_major <- data_sperm[data_sperm$Scaffold %in% major,]
data_sperm_minor <- data_sperm[data_sperm$Scaffold %ni% major,]

data_sperm_minor <- data_sperm_minor %>% group_by(Sample) %>%
  summarise(Scaffold="Minor",
            Length=sum(Length),
            Reads=sum(Reads),
            Sample=first(Sample)) %>%
  select(Scaffold,Length,Reads,Sample)

data_sperm <- rbind(data_sperm_major,data_sperm_minor)

# Now we just plot length versus number of reads
# As well as the reads/length ratio per scaffold, multiplied by read length to get estimated mean coverage

length_reads_sperm <- ggplot(data_sperm,aes(x=Length,y=Reads,colour=Sample,shape=Scaffold=="Minor")) + 
  geom_point() + geom_line()

mean_reads_sperm <- ggplot(data_sperm,aes(x=Scaffold,y=Reads/Length*150,colour=Sample)) + 
  geom_point() + scale_y_continuous(trans="log10") + 
  scale_x_discrete(labels=c("Minor",1:17)) + ylab("Estimated mean coverage")

# So clearly, the minor scaffolds have more reads than thet should for their combined length
# This makes sense, since they're very repeat rich

# Let's also do this for the adult sequencing as well
# The question is: do the sperm samples have a different minor/major coverage ratio than the whole genome stuff

data_adult <- read.table("coverage_stats_adult.txt",header=T,stringsAsFactors=F,sep="\t")
data_adult_major <- data_adult[data_adult$Scaffold %in% major,]
data_adult_minor <- data_adult[data_adult$Scaffold %ni% major,]
data_adult_minor <- data_adult_minor %>% group_by(Sample) %>%
  summarise(Scaffold="Minor",
            Length=sum(Length),
            Reads=sum(Reads),
            Sample=first(Sample)) %>%
  select(Scaffold,Length,Reads,Sample)
data_adult <- rbind(data_adult_major,data_adult_minor)

length_reads_adult <- ggplot(data_adult,aes(x=Length,y=Reads,colour=Sample,shape=Scaffold=="Minor")) + 
  geom_point() + geom_line()
mean_reads_adult <- ggplot(data_adult,aes(x=Scaffold,y=Reads/Length*150,colour=Sample)) + 
  geom_point() + scale_y_continuous(trans="log10") + 
  scale_x_discrete(labels=c("Minor",1:17)) + ylab("Estimated mean coverage")

# Save all png files

png("length_v_reads_sperm.png",res=300,width=6,height=6,units='in')
plot(length_reads_sperm)
dev.off()
png("length_v_reads_adult.png",res=300,width=6,height=6,units='in')
plot(length_reads_adult)
dev.off()
png("mean_reads_sperm.png",res=300,width=6,height=6,units='in')
plot(mean_reads_sperm)
dev.off()
png("mean_reads_adult.png",res=300,width=6,height=6,units='in')
plot(mean_reads_adult)
dev.off()

# Very quickly, calculate per sample the ratio of major v minor reads/length

# Sperm

data_sperm_major_ratio <- data_sperm_major %>% mutate(Ratio=Reads/Length) %>%
  group_by(Sample) %>%
  summarise(Sample=first(Sample),
            Mean_Ratio_Major=mean(Ratio))
data_sperm_minor_ratio <- data_sperm_minor %>% mutate(Ratio=Reads/Length) %>%
  group_by(Sample) %>%
  summarise(Sample=first(Sample),
            Mean_Ratio_Minor=mean(Ratio))
ratios_sperm <- merge(data_sperm_major_ratio,data_sperm_minor_ratio) %>% 
  mutate(Minor_Major_Mean_Ratio = Mean_Ratio_Minor/Mean_Ratio_Major,Sample="Sperm")

# Adult

data_adult_major_ratio <- data_adult_major %>% mutate(Ratio=Reads/Length) %>%
  group_by(Sample) %>%
  summarise(Sample=first(Sample),
            Mean_Ratio_Major=mean(Ratio))
data_adult_minor_ratio <- data_adult_minor %>% mutate(Ratio=Reads/Length) %>%
  group_by(Sample) %>%
  summarise(Sample=first(Sample),
            Mean_Ratio_Minor=mean(Ratio))
ratios_adult <- merge(data_adult_major_ratio,data_adult_minor_ratio) %>% 
  mutate(Minor_Major_Mean_Ratio = Mean_Ratio_Minor/Mean_Ratio_Major,Sample="Adult")

ratio_plot <- rbind(ratios_sperm,ratios_adult) %>% ggplot(aes(x=Sample,y=Minor_Major_Mean_Ratio)) + geom_point() + 
  ylab("Ratio of mean reads/length in minor versus major scaffolds per sample")

png("ratio_plot.png",res=300,width=6,height=6,units='in')
plot(ratio_plot)
dev.off()
