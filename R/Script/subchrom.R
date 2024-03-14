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

dfclean <- df %>% filter(mode == "SCA") %>% select(seed, type, density, subchromatique)

dffinal <- dfclean %>%
  group_by(seed, type, density) %>%
  summarise(
    mean=mean(subchromatique),
  )

  plot <- ggplot(dffinal) +
      geom_line( aes(x=density, y=mean, group=type, colour=type), linewidth = 0.8, linetype = 1) +
      geom_hline(yintercept=23, linetype = 2, color="purple", size = 0.8) +
      scale_y_continuous(breaks = c(23, 30, 60, 90)) +
      ggtitle("Subchromatique") +
      labs(x = "Density", y = "Nb Subchromatique")
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
    
    pdf(paste("../Figure/subchrom.pdf", sep=""), width = 4.5, height = 2.5)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchrom.tex", sep="")), width = 4.5, height = 2.5)
    print(plot)
    dev.off()
  

















