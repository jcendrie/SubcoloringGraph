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

dfclean <- df %>% filter(mode == "SCA") %>% filter(type != "ST") %>% select(seed, type, density, subchromatique)

dffinal <- dfclean %>%
  group_by(type, density) %>%
  summarise(
    mean=mean(subchromatique),
  )

for (i in 1:2) {
  
  plot <- ggplot(dffinal) +
      geom_line( aes(x=density, y=mean, group=type, colour=type), linewidth = 0.8, linetype = 1) +
    geom_segment(aes(x = 0, y = 23, xend = 1, yend = 23), linetype = 2, color="purple", size = 0.8) +
      scale_color_discrete(name="", breaks=c("QUDG", "QUBG", "ER", "SBM"), labels=c("Quasi-Unit Disk Graphs", "Quasi-Unit Balls Graphs", "Erdös-Rényi", "Stochastic Block Models")) +
      ggtitle("Subchromatic number based on density") +
      labs(x = "Density", y = "Subchromatic number", legends= "Types of graph generation:") +
      guides(color=guide_legend(nrow=2)) +
      theme(legend.position="top", legend.box="vertical", panel.background = element_rect(fill="white"), panel.grid = element_blank())
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
  
    if (i == 1) {
      plot <- plot + scale_y_continuous(breaks = c(23, 50, 75, 100))
      pdf(paste("../Figure/subchrom.pdf", sep=""), width = 5, height = 4)
      print(plot)
      dev.off()
      
      tikz(file = sprintf(paste("../Figure/subchrom.tex", sep="")), width = 5, height = 4)
      print(plot)
      dev.off()
    } else {
      plot <- plot + xlim(0, 0.3) + scale_y_continuous(breaks = c(23, 30, 60), limits = c(0,65))
      pdf(paste("../Figure/subchrom_frag.pdf", sep=""), width = 5, height = 4)
      print(plot)
      dev.off()
      
      tikz(file = sprintf(paste("../Figure/subchrom_frag.tex", sep="")), width = 5, height = 4)
      print(plot)
      dev.off()
    }
}


















