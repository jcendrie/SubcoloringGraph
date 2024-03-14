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



for (i in 1:4) {
  dftemp1 <- dffinal %>% filter(type == c("QUDG", "QUBG", "ER", "SBM")[i])
  for (j in 1:2) {
    dftemp2 <- dftemp1 %>% filter(mode == c("GCA", "SCA")[j])
    
    name = c("QUDG", "QUBG", "ER", "SBM")[i]
    type = c("GCA", "SCA")[j]
    
    
    plot <- ggplot(dftemp2) +
      geom_line(aes(x=density, y=med), linewidth = 0.8, linetype = 1) +
      geom_line(aes(x=density, y=firstq), linewidth = 0.6, linetype = 2) +
      geom_line(aes(x=density, y=lastq), linewidth = 0.6, linetype = 2) +
      geom_line(aes(x=density, y=firstd), linewidth = 0.4, linetype = 3) +
      geom_line(aes(x=density, y=lastd), linewidth = 0.4, linetype = 3) +
      # geom_errorbar( aes(x=density, ymin=firstd, ymax=lastd), width=0.1, colour="orange", alpha=0.9, size=0.8) +
      ggtitle(paste("Density", name, type, sep=" ")) +
      labs(x = "density", y = "Pv")
    # geom_point(cex = 3) +
    # geom_errorbar(aes(xmin=avggp-sdgp, xmax=avggp+sdgp), linewidth = 0.1, width=.1) +
    # geom_errorbar(aes(ymin=avgyield-sdyield, ymax=avgyield+sdyield), linewidth = 0.1, width=.1) +
    # theme(legend.position="top", legend.box="vertical", legend.margin=margin(), legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7)) +
    # theme(axis.text.x = element_text(size=7)) +
    # theme(axis.text.y = element_text(size=7)) +
    # theme(axis.title.x = element_text(size=10)) +
    # theme(axis.title.y = element_text(size=10)) +
    # guides(colour=guide_legend(nrow=2)) +
    # theme(legend.title=element_blank()) +
    # xlim(0.7, 1) +
    
    pdf(paste("../Figure/density_", name, "_", type ,".pdf", sep=""), width = 4.5, height = 2.5)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/density_", name, "_", type ,".tex", sep="")), width = 4.5, height = 2.5)
    print(plot)
    dev.off()
    
  }
}


















