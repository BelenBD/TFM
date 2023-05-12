# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.7. METABOLITE PROFILING OF FERMENTED WINES


library(vegan)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(igraph)
library(grid)
library(ggplotify)
library(goseq)
library(cowplot)


#load("E:/R/Metabolomics_final.RData")


# KEGG orthology table
ko.n_df <- readRDS("E:/Wineteractions/ko.n_df.rds")
ko.deg_df <- readRDS("E:/Wineteractions/ko.deg_df.rds")

# Sample data
sample_df <- read.table("C:/Users/Windows/OneDrive/Escritorio/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Dom.Genus <- ifelse(sample_df$Genus %in% c("Citeromyces", "Kluyveromyces"), "Other", sample_df$Genus)
sample_df$Dom.Genus <- factor(sample_df$Dom.Genus, levels = c("Hanseniaspora", "Lachancea", "Saccharomyces", "Other"))



col.cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
col.gen <- c("#cc3939", "#8da0cb", "#1b9e77", "gray70")


sgm_group <- data.frame(Group = c("Alcohols", rep("Acidity", 5), "Sugars", "Alcohols", rep("Volatiles", 6)),
                        variable = c("Ethanol", "Acetic_acid", "Lactic_acid", "Tartaric_acid", 
                                     "Citric_acid", "Succinic_acid", "Sugars", "Glycerol",
                                     "Ethyl.acetate", "Fusel.alcohol.acetates", "Fusel.alcohols", "EEFA", "SCFA", "MCFA"),
                        cols = c(1, rep(3, 5), 2, 1, rep(4, 6)))
#
#### PCA ####

sgm.pca_df <- sample_df[,c(5:16,19:24)]
row.names(sgm.pca_df) <- sample_df[,1]

#
## MICE for filling missing data
mice_pca <- mice::mice(sgm.pca_df, maxit = 999, method = "pmm", seed = 1)
sgm.pca_df2 <- mice::complete(mice_pca, 1)

## DENSITY PLOTS

# Paso 1: extraer los datos originales y los imputados
originales <- sgm.pca_df
imputados <- sgm.pca_df2

if (!is.numeric(originales) || !is.numeric(imputados)) {
  # Convertir los dataframes a valores numéricos
  originales <- sapply(originales, as.numeric)
  imputados <- sapply(imputados, as.numeric)
  
}

dens_plot_all <- function(datos_orig, datos_imp, varnames) {
  par(mfrow = c(3, 2), cex.main = 3, cex.lab = 2, cex.axis = 2, mar = c(10, 10, 5, 1), mgp = c(6, 2, 1), las = 1)
  for (i in 1:length(varnames)) {
    dens_orig <- density(datos_orig[, i], na.rm = TRUE)
    dens_imp <- density(datos_imp[, i], na.rm = TRUE)
    plot(dens_orig, col = "blue", main = varnames[i])
    lines(dens_imp, col = "red")
    legend("topright", legend = c("Original", "Imputed"), col = c("blue", "red"), lty = 1, cex = 2)
  }
}

varnames <- colnames(sgm.pca_df)
dens_plot_all(originales, imputados, varnames)

png("E:/Wineteractions/Results_GM_tf/Density_plots.png", width = 2000, height = 2000, res = 120)
dens_plot_all(originales, imputados, varnames)
dev.off()


## Calculate PCA with filled missing data
sgm_pca <- prcomp(~ ., data = sgm.pca_df2[,c(1,4:8,11:18)], scale. = TRUE)

# Points
sgm_pca.plot <- as.data.frame(scores(sgm_pca))
sgm_pca.plot <- merge(sample_df[,c(1:4,25)], sgm_pca.plot, by.x = "Sample_ID", by.y = "row.names")

# Vectors
sgm_pca.vplot <- scores(sgm_pca, display = "species")
sgm_pca.vplot <- cbind.data.frame(variable = row.names(sgm_pca.vplot), sgm_pca.vplot)
sgm_pca.vplot <- merge(sgm_group, sgm_pca.vplot[,1:3], by = "variable")
sgm_pca.vplot$Group <- factor(sgm_pca.vplot$Group, levels = c("Alcohols", "Sugars", "Acidity", "Volatiles"))

## PCA Plot (all points)
gg.sgm_pca.gen <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Dom.Genus), shape = 21, size = 3) +
  scale_fill_manual(values = col.gen) +
  geom_segment(data = sgm_pca.vplot, aes(x = 0, y = 0, xend = PC1*11, yend = PC2*11, color = Group), 
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) + 
  geom_text(data = sgm_pca.vplot, aes(x = PC1*10.5, y = PC2*11.5, label = variable, color = Group), 
            size = 5, show.legend = FALSE) +
  scale_color_manual(values = c("#bbcc06", "#81c236", "#30b3bf", "#82107a")) +
  xlab(paste("PC1: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21), order = 1), 
         shape = guide_legend(override.aes = list(fill = "black"), order = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        aspect.ratio = -1,
        axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 25, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        legend.text = element_text(size = 19, color = "black", hjust = 0),
        legend.title = element_text(size = 25, color = "black", hjust = 0),
        axis.text.x = element_text(size = 17, color = "black"))


gg.sgm_pca.gen


ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/Metabolites/met_GM_gen.png", gg.sgm_pca.gen, 
       width = 11, height = 10)


gg.sgm_pca.cond <- ggplot() + 
  geom_point(data = sgm_pca.plot, aes(x = PC1, y = PC2, fill = Condition), shape = 21, size = 3) +
  scale_fill_manual(values = col.cond) +
  geom_segment(data = sgm_pca.vplot, aes(x = 0, y = 0, xend = PC1*11, yend = PC2*11, color = Group), 
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) + 
  geom_text(data = sgm_pca.vplot, aes(x = PC1*10.5, y = PC2*11.5, label = variable, color = Group), 
            size = 5, show.legend = FALSE) +
  scale_color_manual(values = c("#bbcc06", "#81c236", "#30b3bf", "#82107a")) +
  xlab(paste("PC1: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((sgm_pca$sdev)^2 / sum((sgm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21), order = 1), 
         shape = guide_legend(override.aes = list(fill = "black"), order = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        aspect.ratio = -1,
        axis.text.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 25, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        legend.text = element_text(size = 19, color = "black", hjust = 0),
        legend.title = element_text(size = 25, color = "black", hjust = 0),
        axis.text.x = element_text(size = 17, color = "black"))


gg.sgm_pca.cond

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/Metabolites/met_GM_cond.png", gg.sgm_pca.cond, 
       width = 11, height = 10)

summary(aov(sgm_pca.plot$PC1 ~ sgm_pca.plot$Dom.Genus*sgm_pca.plot$Condition))
summary(aov(sgm_pca.plot$PC2 ~ sgm_pca.plot$Dom.Genus*sgm_pca.plot$Condition))

sgm_sc.plot <- sgm_pca.plot[sgm_pca.plot$Dom.Genus=="Saccharomyces",]

summary(aov(sgm_sc.plot$PC1 ~ sgm_sc.plot$Condition))
summary(aov(sgm_sc.plot$PC2 ~ sgm_sc.plot$Condition))

save.image("E:/R/Metabolomics_final.RData")
