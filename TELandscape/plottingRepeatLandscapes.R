#!/usr/bin/env R

# Author: Oliver Stuart
# Date: 17/12/2021

# This script processes a divergence summary file from RepeatMasker
# It plots repeat landscapes  in a few different forms

##################################
######## Initial Checks ##########
##################################

#Are all required packages installed?
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
packages <- c("dplyr","tidyr","ggplot2")
if(sum(is.installed(packages)) < length(packages)){
  print(paste0("Something is missing, check that all required packages (",paste(packages,collapse=", "),") are installed."))
  stop()
}
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggplot2"))

#Bring in the inputs
file <- commandArgs(trailingOnly=T)

#Do the inputs actually exist?
if (length(file)==0){
  stop("Provide a divergence summary file from RepeatMasker.", call.=FALSE)
}

##################################
######### Preparing data #########
##################################

# Read file
lines <- readLines(file)

# Get prefix
prefix <- strsplit(file,"\\.")[[1]][1]

# Process text into a table
data <- read.table(text = lines[(grep("Coverage for each repeat class and divergence",lines)+1):length(lines)],
                        header=T,stringsAsFactors=F,sep=" ")
# Remove mystery NA column
data <- data[,-ncol(data)]

#Gather
data.gath <- data %>% gather(key="Class",value="Divergence",-Div)

# Plot without aggregating by class
full.plot <- data.gath %>%
  ggplot(aes(x=Div,y=Divergence,fill=Class)) + geom_col()

# Aggregate by class
data.gath$Class.Agg <- unlist(lapply(strsplit(data.gath$Class,"\\."),`[[`,1))
agg.plot <- data.gath %>%
  ggplot(aes(x=Div,y=Divergence,fill=Class.Agg)) + geom_col()

# Facet by type
agg.plot.facet <- data.gath %>%
  ggplot(aes(x=Div,y=Divergence,fill=Class.Agg)) + geom_col() + 
  facet_wrap(~Class.Agg,scales="free_y")

png(paste0(prefix,"_fullPlot.png"),res=300,width=15,height=10,units='in')
plot(full.plot)
dev.off()

png(paste0(prefix,"_aggPlot.png"),res=300,width=15,height=10,units='in')
plot(agg.plot)
dev.off()

png(paste0(prefix,"_aggPlotFacet.png"),res=300,width=15,height=10,units='in')
plot(agg.plot.facet)
dev.off()
