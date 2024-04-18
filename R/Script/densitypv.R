library(tikzDevice)
library(ggplot2)
library(tidyverse)
library(dplyr)
require(data.table)
library(colorspace)
require(ggh4x)

# Setting the source file directory as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data = list.files(path = "../../Code/data",pattern="^score_(.*)_(.*)_(.*)_(.*).csv", full.names = TRUE)

# We import the data
df = do.call(rbind, lapply(data, fread))
df <- as.data.frame(unclass(df))


df$density <- as.double(as.integer(df$density*100+0.1))/100.

dfclean <- df %>% filter(type != "ST", subchromatique <= 23) %>% select(seed, type, mode, score, density)


# Calculates mean, sd, se and IC (Standard Deviation, Standard Error, Confidence Interval)
dffinal <- dfclean %>%
  group_by(type, mode, density) %>%
  summarise(
    firstd = quantile(score, probs = 0.1),
    lastd = quantile(score, probs = 0.9),
    firstq = quantile(score, probs = 0.25),
    lastq = quantile(score, probs = 0.75),
    med = quantile(score, probs = 0.5),
    n=n(),
    mean=mean(score)
    # sd=sd(score)
  )
# %>%
  # mutate(se=sd/sqrt(n))  %>%
  # mutate(ic=se * qt((1-0.05)/2 + .5, n-1))



for (i in 1:2) {
  dftemp1 <- dffinal %>% filter(type == c("QUDG", "QUBG", "ER", "SBM")[i])
  # for (j in 1:2) {
  #   dftemp2 <- dftemp1 %>% filter(mode == c("GCA", "SCA")[j])
    
  name = c("Quasi-Unit Disk Graphs", "Quasi-Unit Balls Graphs", "Erdös-Rényi", "Stochastic Block Models")[i]
  # type = c("GCA", "SCA")[j]
    
    
    plot <- ggplot(dftemp1) +
      geom_line(aes(x=density, y=med, group=mode, color=mode), linewidth = 0.8, linetype = 1) +
      geom_line(aes(x=density, y=firstq, group=mode, color=mode), linewidth = 0.6, linetype = 2) +
      geom_line(aes(x=density, y=lastq, group=mode, color=mode), linewidth = 0.6, linetype = 2) +
      geom_line(aes(x=density, y=firstd, group=mode, color=mode), linewidth = 0.4, linetype = 3) +
      geom_line(aes(x=density, y=lastd, group=mode, color=mode), linewidth = 0.4, linetype = 3)
    
    plot <- plot +
      scale_color_manual(name="", values=c("forestgreen", "red"),labels=c("GCA","SCA")) +
      # scale_shape_manual(name = "Lines", breaks = c("linetype", "linetype", "linetype"), values=c(1, 2, 3),labels=c("Mediane","Quartile", "Decile")) +
      theme(legend.position="top", legend.box="vertical", panel.background = element_rect(fill="white"), panel.grid = element_blank()) +
      labs(x = "Density", y = "Pv")
    
      
    plot <- plot + ggtitle(paste("Density", name, sep=" "))

    pdf(paste("../Figure/density_", name, ".pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/density_", name, ".tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
    

}


















