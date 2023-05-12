# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.7. METABOLITE PROFILING OF FERMENTED WINES


library(ggplot2)
library(dplyr)
library(cowplot)

#load("E:/R/pca_cond.RData")


sample_df <- read_csv("C:/Users/Windows/OneDrive/Escritorio/sample_GM_volat.csv")
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Origin <- factor(sample_df$Origin, levels = c("LM", "M", "R1", "R2", "R3A", "R3B", "R3C", "RdG", "VLP"))
sample_df$Farming <- factor(sample_df$Farming, levels = c("CONV", "ECO"))
#sample_df <- subset(sample_df, Variety != 'Garnacha') #Si eliminamos la variedad Garnacha, cuyas muest

# Missing data imputation
mice_pca <- mice::mice(sample_df[, 21:25], maxit = 999, method = "pmm", seed = 1) # set.seed para que siempre reemplace los missing values por los mismos valores. Maxit son las iteraciones, no hace falta hacer tantas, pero cuantas más iteraciones mejor
gm_pca_df <- mice::complete(mice_pca, 1)

originales <- sample_df[, 21:25]
imputados <- gm_pca_df[, 21:25]

if (!is.numeric(originales) || !is.numeric(imputados)) {
  # Convertir los dataframes a valores numéricos
  originales <- sapply(originales, as.numeric)
  imputados <- sapply(imputados, as.numeric)

}

# Create density plot of imputed data distribution
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

varnames <- colnames(sample_df[, 21:25])
dens_plot_all(originales, imputados, varnames)
png("E:/Wineteractions/Results_GM_tf/Density_plots.png", width = 2000, height = 2000, res = 120)
dens_plot_all(originales, imputados, varnames)
dev.off()


# Metabolite analysis

gm_pca_df <- cbind.data.frame(sample_df[, 1:20], gm_pca_df)
row.names(gm_pca_df) <- gm_pca_df[, 1]

gm_pca <- prcomp(gm_pca_df[,c(6:11, 14:24)], scale = TRUE)

gm_pca.plot <- as.data.frame(gm_pca$x)
gm_pca.plot$Sample_ID <- gm_pca_df$Sample_ID
row.names(gm_pca.plot) <- gm_pca_df[, 1]
gm_pca.plot <-gm_pca.plot %>%
  relocate(Sample_ID)

gm_pca.plot <- cbind(gm_pca_df, gm_pca.plot[, -1])

PCA_origin <- ggplot(gm_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Origin, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                         "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
                                           xlab(paste("PC1: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 18, color = "black", margin = margin(t = 10, r = 0, b = 30, l = 0)),
        axis.title.y = element_text(size = 18, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 30)),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"))

PCA_origin

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/Metabolites/GM_origin.png", PCA_origin, bg = "white",
       width = 11, height = 7)

PCA_cond <- ggplot(gm_pca.plot) + 
  geom_point(aes(x = PC1, y = PC2, color = Condition, shape = Farming), size = 4) +
  scale_color_manual(values = c("#2ac219", "#1949c2", "#dba54d", "#e02424", 
                                "#c124e0", "#89209e", "#a6165c", "#750f41", "#5c105e")) +
  xlab(paste("PC1: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[1], "%", sep = "")) +
  ylab(paste("PC2: ", round(((gm_pca$sdev)^2 / sum((gm_pca$sdev)^2))*100, 2)[2], "%", sep = "")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 18, color = "black", margin = margin(t = 10, r = 0, b = 30, l = 0)),
        axis.title.y = element_text(size = 18, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 30)),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"))

PCA_cond


ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/Metabolites/GM_cond.png", PCA_cond, bg = "white",
       width = 11, height = 7)

### ANOVA PC1, PC2, Origin, Farming y Condition

summary(aov(gm_pca.plot$PC1 ~ gm_pca.plot$Origin + gm_pca.plot$Farming + gm_pca.plot$Condition))# + gm_pca.plot$Variety))
summary(aov(gm_pca.plot$PC2 ~ gm_pca.plot$Origin + gm_pca.plot$Farming + gm_pca.plot$Condition))# + gm_pca.plot$Variety))

#save.image("E:/R/pca_cond.RData")

