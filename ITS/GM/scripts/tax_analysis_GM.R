# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.5.1. TAXONOMIC EXPLORATION OF ITS READS (GM)


library(readr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# load("E:/TFM/ITS/GM/tax_analysis/tax_analysis_gm.RData")

# Load output files from dada2_GM.R script
ASV_GM <- readRDS("E:/Results/GM/ASV_GM.rds")
ASV_GM_TotalCounts <- apply(ASV_GM, 1, function(x) x/sum(x)) # % of each sequence for every sample (ITS)
                                                              
tax_GM <- readRDS("E:/Results/GM/tax_GM.rds")
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM)) 
tax_GM_matrix <- as.matrix(tax_GM) 
tax_GM_matrix[is.na(tax_GM_matrix)] <- "Unidentified" 

# Sample data 
sample_df <- read_csv("C:/Users/Windows/OneDrive/Escritorio/sample_GM.txt")
row.names(sample_df) <- sample_df$Seq_ID 
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2")) 
sample_df <- cbind.data.frame(Sample_ID = paste(sample_df$Origin, sample_df$Farming, sample_df$Condition, sep = "-"),
                              sample_df) 

sample_df$Tax_rep <- paste(sample_df$Farming, sample_df$Tax_rep, sep = "-")

#### TAXONOMIC EXPLORATION ####

ASV_GM_TotalCounts.p <- melt(ASV_GM_TotalCounts) # Organise all ASVs for each ITS
colnames(ASV_GM_TotalCounts.p) <- c("Id", "Seq_ID", "value") 

# Merge taxonomy and ASVs by genus
ASV_GM_TotalCounts_genus <- merge(ASV_GM_TotalCounts.p, tax_GM_matrix[,c(6,8)], by = "Id")
ASV_GM_TotalCounts_genus <- merge(ASV_GM_TotalCounts_genus, sample_df[c(1:7)], by = "Seq_ID")


# Genus plot
ASV_genus_plot <- aggregate(ASV_GM_TotalCounts_genus$value, list(ASV_GM_TotalCounts_genus$Sample_ID, ASV_GM_TotalCounts_genus$Genus, ASV_GM_TotalCounts_genus$Stage), sum)
colnames(ASV_genus_plot) <- c("Sample_ID", "Genus", "Stage", "value")
ASV_genus_plot$Genus[ASV_genus_plot$value < 0.05] <- "Other"
ASV_genus_plot <- aggregate(ASV_genus_plot$value, list(ASV_genus_plot$Sample_ID, ASV_genus_plot$Genus, ASV_genus_plot$Stage), sum)
colnames(ASV_genus_plot) <- c("Sample_ID", "Genus", "Stage", "value")
ASV_genus_plot$Genus <- gsub("g__", "", ASV_genus_plot$Genus)

orderG <- levels(factor(ASV_genus_plot$Genus)) 
orderG <- orderG[! orderG %in% c("Other", "Unidentified")]
orderG <- append(orderG, c("Other", "Unidentified"))

ASV_genus_plot <- merge(ASV_genus_plot, sample_df[,1:7], by = c("Sample_ID", "Stage")) 
ASV_genus_plot$Origin <- factor(ASV_genus_plot$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", "R3B", "R3C"))
ASV_genus_plot$Sample_name <- paste(ASV_genus_plot$Farming, ASV_genus_plot$Condition, sep = "-")
ASV_genus_plot$Sample_name <- factor(ASV_genus_plot$Sample_name, levels = c("CONV-Control", "CONV-18C", "CONV-NH4", "CONV-SO2",
                                                                            "ECO-Control", "ECO-18C", "ECO-NH4", "ECO-SO2"))
# Select plot data for the initial time of fermentation
ASV_genus_plot_t0<- subset(ASV_genus_plot, ASV_genus_plot$Stage=='0_initial')
unique(ASV_genus_plot_t0$Genus)

# set.seed(42) # Note that we first generated the graph with random colours, assessed the present genera, and then customised specific colours for them
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
genus_plot<- ggplot(ASV_genus_plot_t0, 
                      aes(x = Sample_name, y = value, fill = factor(Genus, levels = orderG))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_x_discrete(labels = c("I", "II", "III", "IV", "I", "II", "III", "IV"))+ #As we filtered by stage == '0_initial' (i.e., time=0), the experimental conditions
                                                                                # have not been applied yet. Thereby, all samples from the same origin and farming are considered replicates
  annotate(geom = "text", x = 2.5, y = 1.05, label = "CONV", size=4) +
  annotate(geom = "text", x = 6.5, y = 1.05, label = "ECO", size=4) +
  scale_fill_manual(name= "Genus", values = c("Acremonium"="#918485", "Alternaria" = "#00C5CD", "Aureobasidium"="#CD96CD", 
                                               "Botrytis"="#663300", "Cladosporium" = "#FFD94A", 
                                               "Diplodia" = "#7A67EE", "Erysiphe" = "#1e74eb", 
                                               "Filobasidium" = "#cc9999", "Lachancea" = "#B0C4DE", "Quambalaria" = "#003333", 
                                               "Saccharomyces" = "#336666", "Stemphylium" = "#9ACD32", "Cytospora" = "lightsalmon2",
                                               "Other" = "#990066")) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        aspect.ratio = 2,
        axis.text.x = element_text(angle = 360, hjust = 0.5, vjust = 3, size = 9, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(vjust = -1, size = 15, color = "black"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  xlab("")+ ylab("Abundance") +
  guides(fill = guide_legend(nrow = 3)) + facet_wrap(~Origin, nrow = 1)

genus_plot

ggsave("E:/TFM/ITS/GM/tax_analysis/Taxa_genus_gm.png", genus_plot, bg = "white",
       width = 11, height = 7)


save.image("E:/TFM/ITS/GM/tax_analysis/tax_analysis_gm.RData")
