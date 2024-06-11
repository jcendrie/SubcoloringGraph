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

dfudg <- dffinal %>% filter(type == "UDG") %>% mutate(UDG = mean) %>% ungroup() %>% select(density, UDG)
dfubg <- dffinal %>% filter(type == "UBG") %>% mutate(UBG = mean) %>% ungroup() %>% select(density, UBG)
dfqudg <- dffinal %>% filter(type == "QUDG") %>% mutate(QUDG = mean) %>% ungroup() %>% select(density, QUDG)
dfqubg <- dffinal %>% filter(type == "QUBG") %>% mutate(QUBG = mean) %>% ungroup() %>% select(density, QUBG)
dfer <- dffinal %>% filter(type == "ER") %>% mutate(ER = mean) %>% ungroup() %>% select(density, ER)
dfsbm <- dffinal %>% filter(type == "SBM") %>% mutate(SBM = mean) %>% ungroup() %>% select(density, SBM)

dftemp1 <- merge(dfudg, dfubg, by=c("density"))
dftemp2 <- merge(dftemp1, dfqudg, by=c("density"))
dftemp3 <- merge(dftemp2, dfqubg, by=c("density"))
dftemp4 <- merge(dftemp3, dfer, by=c("density"))
dffinal2 <- merge(dftemp4, dfsbm, by=c("density"))

for (i in 1:2) {
  
  plot <- ggplot(dffinal) +
    geom_line( aes(x=density, y=mean, group=type, colour=type), linewidth = 0.8, linetype = 1) +
    geom_segment(aes(x = 0, y = 23, xend = 1, yend = 23), linetype = 2, color="purple", size = 0.8) +
    scale_color_discrete(name="", breaks=c("UDG", "UBG", "QUDG", "QUBG", "ER", "SBM"), labels=c("Unit Disk Graphs", "Unit Balls Graphs", "Quasi-Unit Disk Graphs", "Quasi-Unit Balls Graphs", "Erdös-Rényi", "Stochastic Block Models")) +
    # ggtitle("Subchromatic number based on density") +
    ggtitle("") +
    labs(x = "Density", y = "Subchromatic number", legends= "Types of graph generation:") +
    guides(color=guide_legend(nrow=2)) +
    theme(legend.position="top", legend.box="vertical", panel.background = element_rect(fill="white"), legend.key = element_rect(fill="white"), legend.key.size = unit(0.3, 'cm'), panel.grid = element_blank(), 
          panel.grid.major = element_line(color = "grey", size = 0.25, linetype = 1),
          axis.line = element_line(arrow = arrow(angle = 30,
                                                 length = unit(0.2, "cm"),
                                                 ends = "last", 
                                                 type = "closed"), color = "black", size = 0.5, linetype = 1)
    )
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
    plot <- plot + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
    pdf(paste("../Figure/subchrom.pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchrom.tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
  } else {
    plot <- plot + xlim(0, 0.3) + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60), limits = c(0,65))
    pdf(paste("../Figure/subchrom_frag.pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchrom_frag.tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
  }
}

for (i in 1:2) {
  
  plot <- ggplot(dffinal2) +
    geom_line( aes(x=density, y=UDG, color="1"), linewidth = 0.8, linetype = 1) +
    geom_line( aes(x=density, y=UBG, color="2"), linewidth = 0.8, linetype = 1) +
    geom_line( aes(x=density, y=QUDG, color="3"), linewidth = 0.8, linetype = 2) +
    geom_line( aes(x=density, y=QUBG, color="4"), linewidth = 0.8, linetype = 2) +
    geom_line( aes(x=density, y=ER, color="5"), linewidth = 0.8, linetype = 3) +
    geom_line( aes(x=density, y=SBM, color="6"), linewidth = 0.8, linetype = 3) +
    geom_segment(aes(x = 0, y = 23, xend = 1, yend = 23), linetype = 2, color="purple", size = 0.6) +
    scale_color_manual(name="", values=c("orange", "green", "red", "forestgreen", "brown", "cyan"), labels=c("Unit Disk Graphs", "Unit Balls Graphs", "Quasi-Unit Disk Graphs", "Quasi-Unit Balls Graphs", "Erdös-Rényi", "Stochastic Block Models")) +
    # ggtitle("Subchromatic number based on density") +
    ggtitle("") +
    labs(x = "Density", y = "Subchromatic number", legends= "Types of graph generation:") +
    guides(color=guide_legend(nrow=2)) +
    theme(legend.position="top", legend.box="vertical", panel.background = element_rect(fill="white"), legend.key = element_rect(fill="white"), legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size=7), axis.title = element_text(size=15), axis.text = element_text(size=10), panel.grid = element_blank(), 
          panel.grid.major = element_line(color = "grey", size = 0.25, linetype = 1),
          axis.line = element_line(arrow = arrow(angle = 30,
                                                 length = unit(0.2, "cm"),
                                                 ends = "last", 
                                                 type = "closed"), color = "black", size = 0.5, linetype = 1)
    )
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
    plot <- plot + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
    pdf(paste("../Figure/subchromV2.pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchromV2.tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
  } else {
    plot <- plot + xlim(0, 0.3) + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60), limits = c(0,65))
    pdf(paste("../Figure/subchrom_fragV2.pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchrom_fragV2.tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
  }
}

for (i in 1:2) {
  
  plot <- ggplot(dffinal2) +
    geom_line( aes(x=density, y=UDG, color="1"), linewidth = 0.8, linetype = 1) +
    geom_line( aes(x=density, y=UBG, color="2"), linewidth = 0.8, linetype = 1) +
    geom_line( aes(x=density, y=ER, color="5"), linewidth = 0.8, linetype = 3) +
    geom_line( aes(x=density, y=SBM, color="6"), linewidth = 0.8, linetype = 3) +
    geom_segment(aes(x = 0, y = 23, xend = 1, yend = 23), linetype = 2, color="purple", size = 0.6) +
    scale_color_manual(name="", values=c("orange", "green", "brown", "cyan"), labels=c("Unit Disk Graphs", "Unit Balls Graphs", "Erdös-Rényi", "Stochastic Block Models")) +
    # ggtitle("Subchromatic number based on density") +
    ggtitle("") +
    labs(x = "Density", y = "Subchromatic number", legends= "Types of graph generation:") +
    guides(color=guide_legend(nrow=2)) +
    theme(legend.position="top", legend.box="vertical", panel.background = element_rect(fill="white"), legend.key = element_rect(fill="white"), legend.key.size = unit(0.8, 'cm'), legend.text = element_text(size=7), axis.title = element_text(size=15), axis.text = element_text(size=10), panel.grid = element_blank(), 
          panel.grid.major = element_line(color = "grey", size = 0.25, linetype = 1),
          axis.line = element_line(arrow = arrow(angle = 30,
                                                 length = unit(0.2, "cm"),
                                                 ends = "last", 
                                                 type = "closed"), color = "black", size = 0.5, linetype = 1)
    )
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
    plot <- plot + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
    pdf(paste("../Figure/subchromV3.pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchromV3.tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
  } else {
    plot <- plot + xlim(0, 0.3) + scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60), limits = c(0,65))
    pdf(paste("../Figure/subchrom_fragV3.pdf", sep=""), width = 5, height = 4)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/subchrom_fragV3.tex", sep="")), width = 5, height = 4)
    print(plot)
    dev.off()
  }
}

















