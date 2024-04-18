  library(tikzDevice)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  require(data.table)
  library(colorspace)
  require(ggh4x)
  
  # Setting the source file directory as the working directory
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  namefilein = c("CompGFull.csv", "CompGPart.csv", "CompSFull.csv", "CompSPart.csv")
  namefileout = c("GCA_Full", "GCA_Part", "SCA_Full", "SCA_Part")
  
  for (i in 1:4) {
    
    data = list.files(path = "../..",pattern=namefilein[i], full.names = TRUE)
    
    # We import the data
    df = do.call(rbind, lapply(data, fread))
    df <- as.data.frame(unclass(df))
    
    color1 <- "darkblue"
    if (i <= 2) {
      color2 <- "forestgreen"
    } else {
      color2 <- "red"
    }
    
    df$nb <- 0
    
    df$nb <- which(df$nb==0)
    
    df1 <- df %>% select(nb, ns3) %>% rename(score = ns3)
    df1$type <- "ns3"
    df2 <- df %>% select(nb, theo) %>% rename(score = theo)
    if (i <= 2) {
      df2$type <- "GCA"
    } else {
      df2$type <- "SCA"
    }
  
    dffinal <- merge(df1, df2, all=TRUE)
    
    if (i <= 2) {
      dffinal$type <- factor(dffinal$type, c("ns3", "GCA"))
    } else {
      dffinal$type <- factor(dffinal$type, c("ns3", "SCA"))
    }
    
    
    plot <- ggplot(dffinal) +
      geom_bar( aes(x=nb, y=score, fill=type), stat="identity", position=position_dodge(width = 0.9), alpha=0.5) +
      scale_fill_manual(values=c(color1, color2)) +
      theme(legend.position="top", legend.box="vertical", legend.margin=margin(), legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7), panel.background = element_rect(fill="white"), panel.grid = element_blank()) +
      scale_x_continuous(breaks = seq(1, c(168, 18, 168, 14)[i], by = c(3, 1, 3, 1)[i])) +
        labs(x = "Nodes", y = "Pv", legends="")
    
    pdf(paste("../Figure/", namefileout[i], "v2.pdf", sep=""), width = c(5, 5, 5, 5)[i], height = 3.5)
    print(plot)
    dev.off()
    
    tikz(file = sprintf(paste("../Figure/", namefileout[i], "v2.tex", sep="")), width = c(5, 5, 5, 5)[i], height = 3.5)
    print(plot)
    dev.off()
  
  }
  
  
  
  
