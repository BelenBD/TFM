# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.6: DIVERSITY ANALYSIS

library(ggplot2)
library(reshape2)
library(vegan)
library(hillR)
library(cowplot)
library(ggpubr)
library(data.table)
library(agricolae)
library(sf)
library(geosphere)
library(msa)
library(phangorn)
library(DECIPHER)
library(phyloseq)


#load("E:/R/Diversity.RData")


# Set taxonomic names as "Unidentified" in those cases where a specific taxonomic name has not been provided
set_unid <- function(tax_df) {
  
  for (i in 1:ncol(tax_df)) {
    if (i != 7) {
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  paste("Unidentified", tax_df[,i-1], sep = " ")),
                           tax_df[,i])
    } else{
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  "sp."),
                           tax_df[,i])
      
    }
    
  }
  
  return(tax_df)
  
}


# Sample data

sample_df <- read.table("C:/Users/Windows/OneDrive/Escritorio/sample_GM.txt", sep = "\t", header = TRUE)
row.names(sample_df) <- sample_df$Seq_ID
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("RdG", "VLP", "LM", "M", "R1", "R2", "R3A", 
                                                        "R3B", "R3C"))
sample_df$Region <- factor(sample_df$Region, levels = c("Ribera del Guadiana", "Valdepenas", "La Mancha",
                                                        "Madrid", "Rioja"))

sample_df <- cbind.data.frame(Sample_ID = paste(sample_df$Origin, sample_df$Farming, sample_df$Condition, sep = "-"),
                              sample_df)

sample_df <- subset(sample_df, Stage == "0_initial")

## Community data
asv_GM <- readRDS("E:/R/ASV_GM.rds")
asv_GM <- asv_GM[row.names(asv_GM) %in% sample_df$Seq_ID, ]
asv_GM <- asv_GM[, colSums(asv_GM != 0) > 0]

asv.t_GM <- apply(asv_GM, 1, function(x) x/sum(x)) #MARGIN=1 es que la funcion se aplica sobre las filas. Lo que no pillo es por que se le da la vuelta.

tax_GM <- as.matrix(readRDS("E:/R/tax_GM.rds"))
tax_GM <- cbind.data.frame(tax_GM, Id = row.names(tax_GM))
tax_GM <- tax_GM[row.names(tax_GM) %in% colnames(asv_GM), ]

tax_GM[is.na(tax_GM)] <- "Unidentified"


#### ALPHA DIVERSITY ANALYSIS ####

alpha_GM <- cbind.data.frame(    
  # Hill based taxonomic alpha diversity
  t.q0 = hill_taxa(t(asv.t_GM), q = 0),
  t.q1 = hill_taxa(t(asv.t_GM), q = 1),
  t.q2 = hill_taxa(t(asv.t_GM), q = 2))

alpha_GM$Seq_ID <- row.names(alpha_GM)
alpha_GM <- merge(alpha_GM, sample_df[, 1:6], by = "Seq_ID")

alpha_GM.plot <- melt(alpha_GM)

## Origin

summary(aov(t.q0~Origin, alpha_GM))
summary(aov(t.q1~Origin, alpha_GM))
summary(aov(t.q2~Origin, alpha_GM))

t.q0 <- as.data.frame(LSD.test(aov(t.q0~Origin, alpha_GM), "Origin")$group)
t.q1 <- as.data.frame(LSD.test(aov(t.q1~Origin, alpha_GM), "Origin")$group)
t.q2 <- as.data.frame(LSD.test(aov(t.q2~Origin, alpha_GM), "Origin")$group)

alpha_GM.plot$LSD.test.groups <- 0

for(i in 1:length(alpha_GM.plot$Origin)) {                           
  
  df <- as.character(alpha_GM.plot$variable[i])
  origin <- as.character(alpha_GM.plot$Origin[i])
  if (df=='t.q0'){
    alpha_GM.plot$LSD.test.groups[i] <- t.q0[origin, 2]
  } else if (df=='t.q1'){
    alpha_GM.plot$LSD.test.groups[i] <- t.q1[origin, 2]
  } else {
    alpha_GM.plot$LSD.test.groups[i] <- t.q2[origin, 2]
  }
  
}

alpha_GM.plot$Alpha_div <- ifelse(alpha_GM.plot$variable == "t.q0", "Richness", 
                                  ifelse(alpha_GM.plot$variable == "t.q1", "Hill-Shannon", "Hill-Simpson"))

alpha_GM.plot$Alpha_div <- reorder(alpha_GM.plot$Alpha_div, 
                                   as.numeric(factor(alpha_GM.plot$Alpha_div, levels = c("Richness", "Hill-Shannon", "Hill-Simpson"))))

alpha_origin <- ggplot(data = alpha_GM.plot, aes(x = Origin, y = value)) +
  geom_boxplot(fatten=1.5) +
  facet_wrap(~Alpha_div, scales = "free_y", nrow = 3) +
  geom_text(data = alpha_GM.plot, aes(x= Origin, y= 0.2, label = LSD.test.groups), show.legend = F, size=4, vjust=0)+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, color = "black"))

alpha_origin

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/alpha_div_origin.png", alpha_origin, bg = "white",
       width = 7, height = 7, dpi = 600)

## Farming
alpha_farming <- ggplot(data = alpha_GM.plot, aes(x = Farming, y = value)) +
  geom_boxplot() +
  facet_wrap(~Alpha_div, scales = "free_y", nrow = 1)+ geom_signif(comparisons = list(c("CONV", "ECO")),test='t.test', textsize=7, vjust=0.5, map_signif_level=TRUE) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, color = "black"))

alpha_farming

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/alpha_div_farming.png", alpha_farming, bg = "white",
       width = 11, height = 5, dpi = 600)

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/alpha_div.png", 
       plot_grid(alpha_origin, alpha_farming, align = "v", labels = c("a", "b"), label_size = 20, nrow = 2), 
       width = 7, height = 10, dpi = 300, bg = "white")


#ANOVA

summary(aov(alpha_GM$t.q0~alpha_GM$Farming * alpha_GM$Origin))
summary(aov(alpha_GM$t.q1~alpha_GM$Farming * alpha_GM$Origin))
summary(aov(alpha_GM$t.q2~alpha_GM$Farming * alpha_GM$Origin))


#### BETA DIVERSITY ANALYSIS ####

# Agglomerating Genus

tax_GM.un <- set_unid(tax_GM)
tax_GM[is.na(tax_GM)] <- "Unidentified"
asv.t_gen <- aggregate(asv.t_GM, list(tax_GM.un[,"Genus"]), sum)
row.names(asv.t_gen) <- asv.t_gen[,1]
asv.t_gen <- asv.t_gen[,-1]

## Calculate Bray-Curtis dissimilarity and NMDS
bray_gen <- vegdist(t(asv.t_gen), method = "bray", binary = TRUE)

set.seed(1)
nMDS_gen <- metaMDS(bray_gen)
nMDS_gen$stress

nMDS_gen.plot <- as.data.frame(nMDS_gen[["points"]])
nMDS_gen.plot$Seq_ID <- rownames(nMDS_gen.plot)
nMDS_gen.plot <- merge(nMDS_gen.plot, sample_df, by = "Seq_ID")

gg.nmds.gen_bray <- ggplot(nMDS_gen.plot) + 
  geom_point(aes(x = MDS1, y = MDS2, color = Origin, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#9fffd1", "#dba54d", "#ae0e36", 
                                         "#4e2469", "#2227a3", "#e8a6ed", "#993f60", "#19b7c2")) +
                                           annotate(geom = "text", x = Inf, y = Inf, label = round(nMDS_gen$stress, 3), size = 6, vjust = 2, hjust = 1.5) +
  theme_bw() +
  theme(aspect.ratio = 0.66,
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"))

gg.nmds.gen_bray

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/beta_div.png", gg.nmds.gen_bray, bg = "white",
       width = 11, height = 7, dpi = 600)

set.seed(1)
adonis2(bray_gen ~ Origin*Farming, data = sample_df)


#save.image("E:/R/Diversity.RData")
