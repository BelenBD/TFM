# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.5.1. TAXONOMIC EXPLORATION OF ITS READS (SGM)

library(readr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)


#load("E:/TFM/ITS/SGM/tax_analysis/tax_analysis_sgm.RData")


# Load output files from dada2_SGM.R script
ASV_sgm <- readRDS("E:/TFM/ITS/SGM/dada2/output/ASV_SGM.rds") 
ASV_sgm_TotalCounts <- apply(ASV_sgm, 1, function(x) x/sum(x)) # Gives ASV proportion for each sample

tax_sgm <- readRDS("E:/TFM/ITS/SGM/dada2/output/tax_SGM.rds") # Taxonomic assignment file from dada2 (ITS reads)
tax_sgm <- cbind.data.frame(tax_sgm, Id = row.names(tax_sgm)) 
tax_sgm_matrix <- as.matrix(tax_sgm) 
tax_sgm_matrix[is.na(tax_sgm_matrix)] <- "Unidentified" 

# Sample data 
sample_df <- read_delim("E:/TFM/ITS/SGM/sample_SGM.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#### TAXONOMIC EXPLORATION ####

ASV_sgm_TotalCounts.p <- melt(ASV_sgm_TotalCounts) 
colnames(ASV_sgm_TotalCounts.p) <- c("Id", "Sample_ID", "value") 
ASV_sgm_TotalCounts.p$Sample_ID <- gsub("NGS025-21-RUN-2-", "", ASV_sgm_TotalCounts.p$Sample_ID)

# Merge taxonomy and ASVs by genus
ASV_sgm_TotalCounts_genus <- merge(ASV_sgm_TotalCounts.p, tax_sgm_matrix[,c(6,8)], by = "Id")
ASV_sgm_TotalCounts_genus <- merge(ASV_sgm_TotalCounts_genus, sample_df[c(1:6)], by = "Sample_ID")


# Genus plot

#Hacemos una tabla con los valores de porporción de frecuencia de ASVs (value -> x), etiqueta (condiciones) de cada muestra (Sample_ID), género y etapa temporal 
ASV_genus_plot <- aggregate(ASV_sgm_TotalCounts_genus$value, list(ASV_sgm_TotalCounts_genus$Sample_ID, ASV_sgm_TotalCounts_genus$Genus), sum)
colnames(ASV_genus_plot) <- c("Sample_ID", "Genus", "value")

ASV_genus_plot$Genus[ASV_genus_plot$value < 0.05] <- "Other"
ASV_genus_plot <- aggregate(ASV_genus_plot$value, list(ASV_genus_plot$Sample_ID, ASV_genus_plot$Genus), sum)
colnames(ASV_genus_plot) <- c("Sample_ID", "Genus", "value")
ASV_genus_plot$Genus <- gsub("g__", "", ASV_genus_plot$Genus)

orderG <- levels(factor(ASV_genus_plot$Genus)) 
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

ASV_genus_plot <- merge(ASV_genus_plot, sample_df[,1:4], by = "Sample_ID") 
ASV_genus_plot$Origin <- factor(ASV_genus_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))

ASV_genus_plot$Sample_name <- paste(ASV_genus_plot$Farming, ASV_genus_plot$Condition, sep = "-")
ASV_genus_plot$Sample_name <- factor(ASV_genus_plot$Sample_name, levels = c("CONV-Control", "CONV-18C", "CONV-NH4", "CONV-SO2",
                                                                            "ECO-Control", "ECO-18C", "ECO-NH4", "ECO-SO2"))

# set.seed(22) # Note that we first generated the graph with random colours, assessed the present genera, and then customised specific colours for them
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
unique(ASV_genus_plot$Genus)

genus_plot <- ggplot(ASV_genus_plot, 
                     aes(x = Sample_name, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") +
  annotate(geom = "text", x = 2.5, y = 1.1, label = "CONV", size=4) +
  annotate(geom = "text", x = 6.5, y = 1.1, label = "ECO", size=4) +
  scale_fill_manual(name= "Genus", values = c("Kluyveromyces" = "#B452CD", "Lachancea" = "#B0C4DE", 
                                              "Saccharomyces" = "#336666", "Citeromyces" = "#D1EEEE", 
                                              "Wickerhamomyces" = "#9ACD32", "Meyerozyma" = "#EECFA1", 
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

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/ITS/SGM/Taxa_genus_sgm.png", genus_plot, bg = "white",
       width = 11, height = 5)

#save.image("E:/TFM/ITS/SGM/tax_analysis/tax_analysis_sgm.RData")
