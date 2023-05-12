# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation

# Section: 5.8.2.2. Taxonomic assignment of RNA reads

library(readr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)

#load("E:/TFM/RNA/tax_analysis/Raw_reads/tax_analysis_rawRna.RData")

# Read taxonomy assignment files from bracken
read_brackens <- function(bracken_path, vars, output = "Abundance") {
  brack_df <- NULL
  for (file in list.files(bracken_path)) {
    
    brack_df <- rbind.data.frame(brack_df,
                                 cbind(read.table(paste("E:/bracken_output", file, sep = "/"), sep = "\t", header = TRUE), 
                                       Sample_ID = gsub(".bracken", "", file)))
    
  }
  
  brack_df <- reshape(brack_df[,vars], timevar = vars[1], idvar = vars[2], direction = "wide", v.names = vars[3])
  brack_df[is.na(brack_df)] <- 0
  row.names(brack_df) <- brack_df[,1]
  brack_df <- brack_df[,-1]
  colnames(brack_df) <- gsub(paste(vars[3], ".", sep = ""), "", colnames(brack_df))
  
  if (output == "Abundance") {
    brack_df <- as.data.frame(apply(brack_df, 2, function(x) x / sum(x)))
  } else  if (output == "Raw") {
    brack_df <- brack_df
  }
  
  return(brack_df)
}


br_df <- read_brackens("E:/bracken_output/", vars = c("Sample_ID", "name", "new_est_reads"), output = "Raw")
br_TotalCounts <- apply(br_df, 2, function(x) x/sum(x)) 

#### TAXONOMIC EXPLORATION ####

br_TotalCounts.melted <- melt(br_TotalCounts)
colnames(br_TotalCounts.melted) <- c("Genus", "Sample_ID", "value") 
br_TotalCounts.melted$Sample_ID <- gsub("output_", "", br_TotalCounts.melted$Sample_ID)

br_plot <- aggregate(br_TotalCounts.melted$value, list(br_TotalCounts.melted$Sample_ID, br_TotalCounts.melted$Genus), sum)
colnames(br_plot) <- c("Sample_ID", "Genus", "value")

br_plot$Genus <- factor(br_plot$Genus)
levels(br_plot$Genus) <- c(levels(br_plot$Genus), "Other")
br_plot$Genus[br_plot$value < 0.05] <- "Other"

orderG <- levels(factor(br_plot$Genus)) 
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

sample_df <- read_delim("E:/TFM/ITS/SGM/sample_SGM.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

br_plot <- merge(br_plot, sample_df[,1:4], by = "Sample_ID") 
br_plot$Origin <- factor(br_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
br_plot$Sample_name <- paste(br_plot$Farming, br_plot$Condition, sep = "-")
br_plot$Sample_name <- factor(br_plot$Sample_name, levels = c("CONV-Control", "CONV-18C", "CONV-NH4", "CONV-SO2",
                                                                            "ECO-Control", "ECO-18C", "ECO-NH4", "ECO-SO2"))

unique(br_plot$Genus)
genus_plot <- ggplot(br_plot, 
                     aes(x = Sample_name, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") +
  annotate(geom = "text", x = 2.5, y = 1.1, label = "CONV", size=4) +
  annotate(geom = "text", x = 6.5, y = 1.1, label = "ECO", size=4) +
  scale_fill_manual(name= "Genus", values = c("Kluyveromyces" = "#B452CD", "Lachancea" = "#B0C4DE", 
                                              "Saccharomyces" = "#336666", "Citeromyces" = "#D1EEEE", 
                                              "Fusarium" = "#308344", "Torulaspora" = "#FFFF5B",
                                              "Hanseniaspora" = "#EE2C2C", "Other" = "#990066")) +
  
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        aspect.ratio = 2,
        axis.text.x = element_text(angle = 90, size = 9, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(vjust = -1, size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("") + ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3)) + facet_wrap(~Origin, nrow = 1)


genus_plot

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/RNA/Reads_taxa/genus_RNA.png", genus_plot, bg = "white",
       width = 11, height = 5)


#save.image("E:/TFM/RNA/tax_analysis/Raw_reads/tax_analysis_rawRna.RData")

