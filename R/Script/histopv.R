library(tikzDevice)
library(ggplot2)
library(tidyverse)
library(dplyr)
require(data.table)
library(colorspace)
require(ggh4x)

library(ggpubr) # For get_legend

# Setting the source file directory as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data = list.files(path = "../../Code/data",pattern="^score_(.*)_(.*)_(.*)_(.*).csv", full.names = TRUE)

# We import the data
df = do.call(rbind, lapply(data, fread))
df <- as.data.frame(unclass(df))


df$density <- as.double(as.integer(df$density*100+0.1))/100.

dfclean <- df %>% filter(type != "ST", density == 0.08, subchromatique <= 23) %>% select(seed, type, mode, score)

totnb = as.integer((dfclean %>% group_by(type, mode) %>% summarize(nbtot=n()) %>% select(type, nbtot))[1, 2])/as.integer((dfclean %>% group_by(seed, type, mode) %>% summarize(nbtot=n()) %>% select(seed, type, nbtot))[1, 3])

dfgroup <- dfclean %>% mutate(groupe=as.integer(as.integer(dfclean$score*100 + 9)/10)) %>%
  group_by(seed, type, mode, groupe) %>%
  summarise(nb=n())



# Calculates mean, sd, se and IC (Standard Deviation, Standard Error, Confidence Interval)
dffinal <- dfgroup %>%
  group_by(type, mode, groupe) %>%
  summarise(
    sum=sum(nb),
    wrongn=n(),
    wrongmean=mean(nb),
    wrongsd=sd(nb)
  ) %>%
  mutate(mean=(sum/totnb)) %>%
  mutate(sd=sqrt((wrongsd*wrongsd+wrongmean*wrongmean)*wrongn/totnb-mean*mean)) %>%
  mutate(se=sd/sqrt(totnb))  %>%
  mutate(ic=se * qt((1-0.05)/2 + .5, totnb-1)) %>%
  mutate(groupe=groupe/10) %>% select(groupe, mean, sd, type, mode)

dfnormal <- dffinal %>% mutate(mean = mean/1000) %>% mutate(sd = sd/1000)
dfnormal$groupe[dfnormal$mode == "SCA"] <- dfnormal$groupe[dfnormal$mode == "SCA"]+0.04
dfnormal$groupe <- dfnormal$groupe - 0.02

# dfadd <- dffinal %>% filter(mode == "SCA")
# dfadd$groupe = 0.0
# dfadd$mean = 0.0
# dfadd$sd = 0.0
# dfadd <- distinct(dfadd)
# 
# dffinal2 <- merge(dffinal, dfadd, all=TRUE)

color1 <- "forestgreen"
color2 <- "red"


for (i in 1:4) {
  dftemp1 <- dfnormal %>% filter(type == c("QUDG", "QUBG", "ER", "SBM")[i])
  colors <- c()
  for (j in 0:10) {
    if (Reduce("|", dftemp1$groupe==toString(-0.02+j*0.1))) {
      colors <- c(colors, color1)
    }
  }
  for (j in 0:10) {
    if (Reduce("|", dftemp1$groupe==toString(0.02+j*0.1))) {
      colors <- c(colors, color2)
    }
  }
  # for (j in 1:2) {
    # dftemp2 <- dftemp1 %>% filter(mode == c("GCA", "SCA")[j]) 
    
    name = c("Quasi-Unit Disk Graphs", "Quasi-Unit Balls Graphs", "Erdös-Rényi", "Stochastic Block Models")[i]
    # type = c("GCA", "SCA")[j]
    
    
    plot <- ggplot(dftemp1) +
      geom_bar( aes(x=groupe, y=mean), stat="identity", fill=colors, alpha=0.5) +
      geom_errorbar( aes(x=groupe, ymin=mean-sd, ymax=mean+sd), width=0.04, colour="black", alpha=0.9, linewidth=0.6) +
      ggtitle(paste("Histogram", name, sep=" ")) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(-0.075, 1.075)) +
      theme(legend.position="top", legend.box="vertical", legend.margin=margin(), legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7), panel.background = element_rect(fill="white"), panel.grid = element_blank()) +
      labs(x = "Pv", y = "Proportion", legends="")
    # geom_point(cex = 3) +
    # geom_errorbar(aes(xmin=avggp-sdgp, xmax=avggp+sdgp), linewidth = 0.1, width=.1) +
    # geom_errorbar(aes(ymin=avgyield-sdyield, ymax=avgyield+sdyield), linewidth = 0.1, width=.1) +
    # theme(axis.text.x = element_text(size=7)) +
    # theme(axis.text.y = element_text(size=7)) +
    # theme(axis.title.x = element_text(size=10)) +
    # theme(axis.title.y = element_text(size=10)) +
    # guides(colour=guide_legend(nrow=2)) +
    # theme(legend.title=element_blank()) +
    # xlim(0.7, 1) +

    pdf(paste("../Figure/histo_", name, ".pdf", sep=""), width = 5, height = 3.5)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/histo_", name, ".tex", sep="")), width = 5, height = 3.5)
    print(plot)
    dev.off()
    
  }
# }

dflegend <- dfnormal %>% group_by(mode) %>% summarise(mean=mean(mean), groupe=mean(groupe))
thelegendasaggplot <- ggplot(dflegend, aes(x=groupe, fill = mode)) + 
  geom_histogram(alpha = 0.5) + 
  theme(legend.position="top", legend.box="vertical", legend.margin=margin(), legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7)) +
  scale_fill_manual(name="", values=c("forestgreen","red"),labels=c("GCA","SCA"))

thelegendasagtable <- get_legend(thelegendasaggplot)
#Convert to a ggplot and print
thelegendasaggplot <- as_ggplot(thelegendasagtable)
pdf(paste("../Figure/histo_legend.pdf", sep=""), width = 8, height = 0.5)
print(thelegendasaggplot)
dev.off()
tikz(file = sprintf(paste("../Figure/histo_legend.tex", sep="")), width = 5, height = 0.5)
print(thelegendasaggplot)
dev.off()
















