library(tikzDevice)
library(ggplot2)
library(tidyverse)
library(dplyr)
require(data.table)
library(colorspace)
require(ggh4x)


# Setting the source file directory as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data = list.files(path = "../../Code/data",pattern="^score_(.*)__(.*)_(.*).csv", full.names = TRUE)

# We import the data
df = do.call(rbind, lapply(data, fread))
df <- as.data.frame(unclass(df))


df$density <- as.double(as.integer(df$density*100+0.1))/100.

dfclean <- df %>% filter(subchromatique <= 23) %>% select(seed, type, mode, score, density)

totbyseed = as.integer((dfclean %>% group_by(seed, type, mode) %>% summarize(nbtot=n()) %>% select(seed, type, nbtot))[1, 3])

dftotnb <- dfclean %>% group_by(type, mode, density) %>% summarize(nbtot=n()/totbyseed)

dfclean <- merge(dfclean, dftotnb, by=c("type", "mode", "density"))

dfgroup <- dfclean %>% mutate(groupe=as.integer(as.integer(dfclean$score*100 + 9)/10)) %>%
  group_by(seed, type, mode, density, groupe, nbtot) %>%
  summarise(nb=n())

# Calculates mean, sd, se and IC (Standard Deviation, Standard Error, Confidence Interval)
dffinal <- dfgroup %>%
  group_by(type, mode, density, groupe, nbtot) %>%
  summarise(
    sum=sum(nb),
    wrongn=n(),
    wrongmean=mean(nb),
    wrongsd=sd(nb)
  ) %>%
  mutate(mean=sum/nbtot) %>%
  mutate(sd=sqrt((wrongsd*wrongsd+wrongmean*wrongmean)*wrongn/nbtot-mean*mean)) %>%
  # mutate(se=sd/sqrt(nbtot))  %>%
  # mutate(ic=se * qt((1-0.05)/2 + .5, nbtot-1)) %>%
  mutate(groupe=groupe/10)



  for (j in 1:2) {
    dftemp2 <- dffinal %>% filter(mode == c("GCA", "SCA")[j])
    
    type = c("GCA", "SCA")[j]
    
    
    plot <- ggplot(dftemp2) +
      geom_bar( aes(x=groupe, y=mean/168), stat="identity", fill="forestgreen", alpha=0.5) +
      geom_errorbar( aes(x=groupe, ymin=(mean-sd)/168, ymax=(mean+sd)/168), width=0.1, colour="orange", alpha=0.9, size=0.8) +
      ggtitle(paste("Histogram CiscoMeraki", type, sep=" ")) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(-0.05, 1.05)) +
      labs(x = "Pv", y = "Nb")

        pdf(paste("../Figure/histo_meraki_", type ,".pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/histo_meraki_", type ,".tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
    
  }


















