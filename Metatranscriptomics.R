# Author: Benítez Domínguez, Belén
# Date: 11/05/2023
# Project: TFM: Wine yeast community assembly: effects of regionality, agricultural 
# management, and fermentative conditions on their composition and transcriptional regulation


# Section: 5.8.2.7. Differential expression analysis and biological enrichment.

library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggforce)
library(clusterProfiler)
library(reshape2)
library(goseq)
library(vegan)


#load("E:/R/Metatranscriptomics.RData")

## Load Count data + Emapper Annotations to get KEGG Ortholog counts
read.emapper <- function(path.count, path.orth) {
  
  count.list <- NULL
  
  for (smpl in list.files(path.count)) {
    sample <- gsub(".txt", "", smpl)
    
    count_df <- read.table(paste(path.count, smpl, sep = "/"), header = TRUE)[,c(1,7)]
    colnames(count_df) <- c("seed_ortholog", sample)
    
    orth_file <- list.files(paste(path.orth, sample, sep = "/"), full.names = TRUE)
    orth_df <- read.delim(orth_file, sep = "\t", skip = 4)
    orth_df <- unique(orth_df[1:(nrow(orth_df)-3), c(2,12)])
    
    count.deff_df <- merge(count_df, orth_df, by = "seed_ortholog")
    count.deff_df <- subset(count.deff_df, KEGG_ko != "-")
    
    count.deff_df <- aggregate(count.deff_df[,2], list(count.deff_df$KEGG_ko), sum)
    colnames(count.deff_df) <- c("KEGG_ko", sample)
    
    count.list[[sample]] <- count.deff_df
    
  }
  
  ko_df <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), count.list)
  ko_df[is.na(ko_df)] <- 0
  
  return(ko_df)
  
}

## Load Gene Ontology Data
read.GOdata <- function(file.ext, path.annot, cols) {
  
  kegg.go_df <- NULL
  
  for (smpl in list.files(path.annot)) {
    
    #Open Annotation files and select KO and GO columns
    annot_df <- read.delim(list.files(paste(path.annot, smpl, sep = "/"), full.names = TRUE), sep = "\t", skip = 4)[,cols]
    annot_df <- unique(annot_df[annot_df$KEGG_ko != "" & annot_df$KEGG_ko != "-", ])
    
    kegg.go_df <- rbind.data.frame(kegg.go_df, annot_df)
    kegg.go_df <- unique(kegg.go_df)
    
  }
  
  #Format the table into a useful form
  go_df <- NULL
  for (n in 1:nrow(annot_df)) {
    
    go_df <- rbind.data.frame(go_df, cbind(gene_id = annot_df$KEGG_ko[n],
                                           category = unlist(strsplit(as.character(annot_df$GOs[n]),",", fixed = TRUE))))
    
    go_df <- unique(go_df)
    
  }
  
  return(go_df)
  
}

# Load data

ko_df <- read.emapper(path.count = "E:/Wineteractions/Count", 
                      path.orth = "E:/Wineteractions/Annotation")

ko_df.f <- ko_df[!grepl(",", ko_df$KEGG_ko),] # delete rows with more than one assigned KEGG_ko


# Annotation
kegg_df <- read.table("E:/Wineteractions/KEGG_names.txt", 
                      sep = "\t", header = TRUE, quote = "") 

go_df <- read.GOdata(".emapper.annotations", "E:/Wineteractions/Annotation", c("KEGG_ko", "GOs"))

saveRDS(go_df, "E:/Wineteractions/go_df.rds")

# Sample data 
sample_df <- read.table("C:/Users/Windows/OneDrive/Escritorio/sample_SGM.txt", sep = "\t", header = TRUE)
sample_df$Condition <- factor(sample_df$Condition, levels = c("Control", "18C", "NH4", "SO2"))
sample_df$Genus <- factor(sample_df$Genus, levels = c("Citeromyces", "Hanseniaspora", "Kluyveromyces", "Lachancea", 
                                                      "Saccharomyces", "Other"))


sample_df$Specie <- gsub("Hanseniaspora ", "H.", sample_df$Specie)

## COLOURS
col.cond <- c("#bf2c45", "#1e74eb", "#ebb249", "#93bf2c")
col.gen <- c("#caf55d", "#cc3939", "#cab2d6", "#8da0cb", "#1b9e77", "gray70")



#### DIFFERENTIAL ANALYSIS - GLOBAL ####
dds_tot <- DESeqDataSetFromMatrix(countData = ko_df.f, 
                                  colData = sample_df, 
                                  design = ~ Condition + Genus + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_tot)) >= 10
dds_tot <- dds_tot[keep,]

dds_tot <- DESeq(dds_tot, fitType = "local") 

## PCA
vst_tot <- vst(dds_tot, blind = TRUE) 
pcaData_tot <- plotPCA(vst_tot, intgroup = "Condition", returnData = TRUE)
rownames(pcaData_tot) <- pcaData_tot[, 1]

percentVar_tot <- round(100 * attr(pcaData_tot, "percentVar"), 2)
pcaData_tot <- merge(pcaData_tot, sample_df[,-4], by.x = "name", by.y = "Sample_ID")

gg.pca_tot_cond <- ggplot(pcaData_tot, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_tot[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_tot[2],"% variance")) + 
  theme_bw() + 
  theme(
    aspect.ratio = 1,
    axis.title.y = element_text(size = 17, color = "black"),
    axis.title.x = element_text(size = 17, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 17, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
    legend.key.width = unit(0.5, "cm")
  ) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")

gg.pca_tot_cond


# Samples are separated by dominant species
gg.pca_tot_gen <- ggplot(pcaData_tot, aes(PC1, PC2, color = Genus)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_tot[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_tot[2],"% variance")) + 
  theme_bw() + theme(plot.title = element_text(size=22)) +
  theme(aspect.ratio = 1,
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0,5), "cm")) +
  scale_color_manual(values = col.gen) +
  labs(color = "Genus")

gg.pca_tot_gen

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/RNA/total_PCA.png", 
       plot_grid(gg.pca_tot_gen, gg.pca_tot_cond, align = "v", labels = c("a", "b"), label_size = 20), 
       width = 15.2, height = 5.4, dpi = 300, bg = "white")

summary(aov(pcaData_tot$PC1~pcaData_tot$Origin * pcaData_tot$Genus))
summary(aov(pcaData_tot$PC1~pcaData_tot$Genus * pcaData_tot$Condition))

summary(aov(pcaData_tot$PC2~pcaData_tot$Origin * pcaData_tot$Genus))
summary(aov(pcaData_tot$PC2~pcaData_tot$Genus * pcaData_tot$Condition))



#### DIFFERENTIAL ANALYSIS - GLOBAL.Genus ####

# KEEP SAMPLES DOMINATED BY ONE OF SAID SPECIES
sample_gen <- sample_df[sample_df$Genus %in% c("Saccharomyces", "Lachancea", "Hanseniaspora"),] # Filter sample_df by dominant genus
ko_df.gen <- ko_df.f[,c("KEGG_ko", sample_gen$Sample_ID)] # Filters samples (columns) corresponding with sample_gen

dds_gen <- DESeqDataSetFromMatrix(countData = ko_df.gen, 
                                  colData = sample_gen, 
                                  design = ~ Condition + Origin + Genus, tidy = TRUE)

dds_gen$Genus <- relevel(dds_gen$Genus, "Saccharomyces") 

keep <- rowSums(counts(dds_gen) >= 10) >= 5
dds_gen <- dds_gen[keep,]

dds_gen <- DESeq(dds_gen, fitType = "local", betaPrior = FALSE)
resultsNames(dds_gen)


## DIFFERENTIAL EXPRESSION

# Saccharomyces vs Lachancea

res_lt.sc <- results(dds_gen, contrast = c("Genus", "Lachancea", "Saccharomyces"), alpha = 0.05) 
res_lt.sc <- lfcShrink(dds_gen, res = res_lt.sc, type = "ashr") 

ress_lt.sc <- as.data.frame(res_lt.sc)
ress_lt.sc$KEGG_ko <- row.names(ress_lt.sc)

ress_lt.sc <- merge(ress_lt.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_lt.sc <- ress_lt.sc[order(ress_lt.sc$padj),]
ress_lt.sc$DEO <- ifelse(abs(ress_lt.sc$log2FoldChange) > 1 & ress_lt.sc$padj <= 0.05, 1, 0)

# Saccharomyces vs Hanseniaspora

res_hs.sc <- results(dds_gen, contrast = c("Genus", "Hanseniaspora", "Saccharomyces"), alpha = 0.05)
res_hs.sc <- lfcShrink(dds_gen, res = res_hs.sc, type = "ashr")
summary(res_hs.sc)

ress_hs.sc <- as.data.frame(res_hs.sc)
ress_hs.sc$KEGG_ko <- row.names(ress_hs.sc)

ress_hs.sc <- merge(ress_hs.sc, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_hs.sc <- ress_hs.sc[order(ress_hs.sc$padj),]
ress_hs.sc$DEO <- ifelse(abs(ress_hs.sc$log2FoldChange) > 1 & ress_hs.sc$padj <= 0.05, 1, 0)

# Lachancea and Hanseniaspora comparative of LFC

lfc_hs <- ress_hs.sc$log2FoldChange
lfc_lt <- ress_lt.sc$log2FoldChange
mean_hs <- mean(lfc_hs)
mean_lt <- mean(lfc_lt)
t.test(lfc_hs, lfc_lt)


## SUMMARY
# DE Orthologs - Venn diagram
venn.df_gen <- merge(ress_lt.sc[,c(1, 9)], ress_hs.sc[,c(1, 9)], by = "KEGG_ko")
colnames(venn.df_gen) <- c("KEGG_ko", "Lt.Sc", "Hs.Sc")
venn.df_gen <- as.data.frame(venn.df_gen)

venn.plot_gen <- rbind.data.frame(cbind(Lt.Sc = 0, Hs.Sc = 0, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 0)),
                                  cbind(Lt.Sc = 0, Hs.Sc = 1, 
                                        Counts = sum(venn.df_gen[,2] == 0 & venn.df_gen[,3] == 1)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 0, 
                                        Counts = sum(venn.df_gen[,2] == 1 & venn.df_gen[,3] == 0)),
                                  cbind(Lt.Sc = 1, Hs.Sc = 1, 
                                        Counts = sum(rowSums(venn.df_gen[,-1]) == 2)))



venn.plot_gen <- cbind.data.frame(venn.plot_gen, x = c(2, 1.4, -1.4, 0), y = c(-1.5, 0, 0, 0))

venn.out2 <- data.frame(x = c(-0.75, 0.75), y = c(0, 0), labels = c("Lachancea", "Hanseniaspora"))
venn.out2$labels <- factor(venn.out2$labels, levels = c("Lachancea", "Hanseniaspora"))

gg.venn_gen <- ggplot() +
  geom_circle(data = venn.out2, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_gen, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 13, color  = "black"),
        legend.title  = element_text(size = 15, color  = "black"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  labs(fill = 'Genus')

gg.venn_gen 


ggsave("E:/Wineteractions/Results_Metatranscriptomics/Venn_Hs_Lt.png", gg.venn_gen, 
       width = 7.6, height = 5.4, dpi = 300, bg = "white")


# Accumulated LFC - Histogram
hist_gen <- rbind(cbind(subset(ress_lt.sc, DEO == 1), Genus = "Lachancea"),
                  cbind(subset(ress_hs.sc, DEO == 1), Genus = "Hanseniaspora"))

hist_gen$Genus <- factor(hist_gen$Genus, levels = c("Lachancea", "Hanseniaspora"))

gg.hist_gen <- ggplot(hist_gen, aes(x = abs(log2FoldChange), fill = Genus)) + 
  geom_histogram(binwidth = 2, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE Orthologs") +
  scale_fill_manual(values = c("#8da0cb", "#cc3939"))  +
  theme(aspect.ratio = 1,
        legend.position = "none",
        axis.text.y = element_text(size = 15, color  = "black"),
        axis.title.x = element_text(size = 17, color  = "black"),
        axis.title.y = element_text(size = 17, color  = "black"),
        legend.text = element_text(size = 15, color  = "black"),
        legend.title = element_text(size = 17, color  = "black"),
        axis.text.x = element_text(size = 15, color  = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

gg.hist_gen


ggsave("E:/Wineteractions/Results_Metatranscriptomics/LFC_Hs_Lt.png", gg.hist_gen, 
       width = 12, height = 7.2, dpi = 300, bg = "white")


### EXPORT DEG (differentially expressed genes) LIST 

DEG_lt.sc_df <- subset(ress_lt.sc, padj <= 0.05 & abs(log2FoldChange) >= 1)[,c(1,3)]
DEG_lt.sc_df$OverExpr <- ifelse(DEG_lt.sc_df$log2FoldChange > 0, "Lachancea", "Saccharomyces")

DEG_hs.sc_df <- subset(ress_hs.sc, padj <= 0.05 & abs(log2FoldChange) >= 1)[,c(1,3)]
DEG_hs.sc_df$OverExpr <- ifelse(DEG_hs.sc_df$log2FoldChange > 0, "Hanseniaspora", "Saccharomyces")

ko.deg_df <- merge(DEG_lt.sc_df[,-2], DEG_hs.sc_df[,-2], by = "KEGG_ko", all = TRUE)
ko.deg_df$OverExpr <- ifelse(is.na(ko.deg_df$OverExpr.x) & !is.na(ko.deg_df$OverExpr.y), ko.deg_df$OverExpr.y, 
                             ifelse(is.na(ko.deg_df$OverExpr.y) & !is.na(ko.deg_df$OverExpr.x), ko.deg_df$OverExpr.x, 
                                    ifelse(ko.deg_df$OverExpr.x == "Saccharomyces", ko.deg_df$OverExpr.y, 
                                           ifelse(ko.deg_df$OverExpr.y == "Saccharomyces", ko.deg_df$OverExpr.x, "Lt&Hs"))))

saveRDS(ko.deg_df, "E:/Wineteractions/ko.deg_df.rds")


## BIOLOGICAL ENRICHMENT - Gene Ontology BP
count.bias <- rowSums(ko_df.gen[,-1])
names(count.bias) <- ko_df.gen[,1]

DEG_lt.sc <- as.integer(ress_lt.sc$padj < 0.05 & abs(ress_lt.sc$log2FoldChange) >= 1)
names(DEG_lt.sc) <- ress_lt.sc$KEGG_ko

pwf_lt.sc <- nullp(DEgenes = DEG_lt.sc, bias.data = count.bias[names(DEG_lt.sc)])
go_lt.sc <- goseq(pwf_lt.sc, gene2cat = go_df, test.cats = c("GO:BP")) 
go_lt.sc <- subset(go_lt.sc, ontology == "BP" & numDEInCat > 5)

enrichGO_lt.sc <- go_lt.sc[p.adjust(go_lt.sc$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_hs.sc <- as.integer(ress_hs.sc$padj < 0.05 & abs(ress_hs.sc$log2FoldChange) >= 1)
names(DEG_hs.sc) <- ress_hs.sc$KEGG_ko

pwf_hs.sc <- nullp(DEgenes = DEG_hs.sc, bias.data = jitter(rep(1000, length(DEG_hs.sc))))

go_hs.sc <- goseq(pwf_hs.sc, gene2cat = go_df) 
go_hs.sc <- subset(go_hs.sc, ontology == "BP" & numDEInCat > 5) 
enrichGO_hs.sc <- go_hs.sc[p.adjust(go_hs.sc$over_represented_pvalue, method = "fdr") < 0.05,]

# Plot
go_gen <- rbind.data.frame(cbind(enrichGO_lt.sc, Comparison = "Lachancea"),
                           cbind(enrichGO_hs.sc, Comparison = "Hanseniaspora"))

go_gen$Genus <- factor(go_gen$Comparison, levels = c("Lachancea", "Hanseniaspora"))
go_gen$comp.cat <- ifelse(duplicated(go_gen$term) | duplicated(go_gen$term, fromLast = TRUE), "Both", go_gen$Comparison)

go_gen$term <- factor(go_gen$term, 
                      levels = unique(go_gen$term[order(go_gen$comp.cat, go_gen$numDEInCat, decreasing = TRUE)]))

gg.go_gen <- ggplot(go_gen) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Genus), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#8da0cb", "#cc3939")) +
  coord_flip() +
  theme_bw()+
  theme(legend.position = "bottom",
        aspect.ratio = 0.62,
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black")) + 
  ylab("DE Orthologs")  + xlab("")

gg.go_gen


ggsave("E:/Wineteractions/Results_Metatranscriptomics/genus_enrichGO.png", gg.go_gen, 
       width = 10, height = 9, dpi = 300, bg = "white")

gg.summary_gen <- plot_grid(gg.venn_gen, gg.hist_gen, nrow = 1, rel_widths = c(2, 3), labels = c("a", "b"), label_size = 20)

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/RNA/gen_summary.png", 
       plot_grid(gg.summary_gen, gg.go_gen, nrow = 2, labels = c("", "c"), label_size = 20, 
                 rel_widths = c(1,2)), 
       width = 14, height = 14, dpi = 300, bg = "white")



#### DIFFERENTIAL ANALYSIS - SACCHAROMYCES ####
sample_sc <- sample_df[sample_df$Genus == "Saccharomyces",]
sample_sc <- sample_sc[-c(5,6),] 
ko_sc <- ko_df.f[,c("KEGG_ko", sample_sc$Sample_ID)]

dds_sc <- DESeqDataSetFromMatrix(countData = ko_sc, 
                                 colData = sample_sc, 
                                 design = ~ Condition + Origin, tidy = TRUE)

keep <- rowSums(counts(dds_sc)) >= 10
dds_sc <- dds_sc[keep,]

dds_sc <- DESeq(dds_sc, fitType = "local")
resultsNames(dds_sc)

## PCA
vst_sc <- vst(dds_sc, blind = TRUE)

pcaData_sc <- plotPCA(vst_sc, intgroup = "Condition", returnData = TRUE)
percentVar_sc <- round(100 * attr(pcaData_sc, "percentVar"), 2)
pcaData_sc <- merge(pcaData_sc, sample_sc[,-4], by.x = "name", by.y = "Sample_ID")

gg.pca_sc <- ggplot(pcaData_sc, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_sc[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sc[2],"% variance")) + 
  theme_bw() +

  theme(aspect.ratio = 1, 
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Condition")

gg.pca_sc

gg.pca_sc.or <- ggplot(pcaData_sc, aes(PC1, PC2, color = Origin)) +
  geom_point(size = 4) + 
  xlab(paste0("PC1: ", percentVar_sc[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar_sc[2],"% variance")) + 
  theme_bw() +
  
  theme(aspect.ratio = 1, 
        axis.title.y = element_text(size = 17, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black")) +
  scale_color_manual(values = col.cond) +
  labs(color = "Origin")

gg.pca_sc.or

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/RNA/sc_PCA.png", 
       plot_grid(gg.pca_sc, gg.pca_sc.or, align = "v", labels = c("a", "b"), label_size = 20), 
       width = 10, height = 5, dpi = 300, bg = "white")

#
## DIFFERENTIAL EXPRESSION

# Low Temperature
res_sc.18C <- results(dds_sc, contrast = c("Condition", "18C", "Control"), alpha = 0.05) 
res_sc.18C <- lfcShrink(dds_sc, res = res_sc.18C, contrast = c("Condition", "18C", "Control"), type = "normal") 

ress_sc.18C <- as.data.frame(res_sc.18C)
ress_sc.18C$KEGG_ko <- row.names(ress_sc.18C)

ress_sc.18C <- merge(ress_sc.18C, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.18C <- ress_sc.18C[order(ress_sc.18C$padj),]
ress_sc.18C$DEO <- ifelse(abs(ress_sc.18C$log2FoldChange) > 1 & ress_sc.18C$padj <= 0.05, 1, 0)

# High Ammonia
res_sc.NH4 <- results(dds_sc, contrast = c("Condition", "NH4", "Control"), alpha = 0.05)
res_sc.NH4 <- lfcShrink(dds_sc, res = res_sc.NH4, contrast = c("Condition", "NH4", "Control"), type = "normal")
summary(res_sc.NH4)

ress_sc.NH4 <- as.data.frame(res_sc.NH4)
ress_sc.NH4$KEGG_ko <- row.names(ress_sc.NH4)

ress_sc.NH4 <- merge(ress_sc.NH4, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.NH4 <- ress_sc.NH4[order(ress_sc.NH4$padj),]
ress_sc.NH4$DEO <- ifelse(abs(ress_sc.NH4$log2FoldChange) > 1 & ress_sc.NH4$padj <= 0.05, 1, 0)

# High Sulfites
res_sc.SO2 <- results(dds_sc, contrast = c("Condition", "SO2", "Control"), alpha = 0.05)
res_sc.SO2 <- lfcShrink(dds_sc, res = res_sc.SO2, contrast = c("Condition", "SO2", "Control"), type = "normal")
summary(res_sc.SO2)

ress_sc.SO2 <- as.data.frame(res_sc.SO2)
ress_sc.SO2$KEGG_ko <- row.names(ress_sc.SO2)

ress_sc.SO2 <- merge(ress_sc.SO2, kegg_df, by = "KEGG_ko", all.x = TRUE)
ress_sc.SO2 <- ress_sc.SO2[order(ress_sc.SO2$padj),]
ress_sc.SO2$DEO <- ifelse(abs(ress_sc.SO2$log2FoldChange) > 1 & ress_sc.SO2$padj <= 0.05, 1, 0)


# Conditions comparative of LFC
lfc_18C <- ress_sc.18C$log2FoldChange
lfc_NH4 <- ress_sc.NH4$log2FoldChange
lfc_SO2 <- ress_sc.SO2$log2FoldChange
mean_18C <- mean(lfc_18C)
mean_NH4 <- mean(lfc_NH4)
mean_SO2 <- mean(lfc_SO2)
anova_res <- aov(c(lfc_18C, lfc_NH4, lfc_SO2) ~ c(rep("18C", length(lfc_18C)), rep("NH4", length(lfc_NH4)), rep("SO2", length(lfc_SO2))))
summary(anova_res)


## SUMMARY
# DE Orthologs - Venn diagram
venn.df_sc <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "KEGG_ko"), 
                     list(ress_sc.18C[,c(1,10)], ress_sc.NH4[,c(1,10)], ress_sc.SO2[,c(1,10)]))
colnames(venn.df_sc) <- c("KEGG_ko", "Sc.18C", "Sc.NH4", "Sc.SO2")
venn.df_sc[is.na(venn.df_sc)] <- 0

venn.plot_sc <- rbind.data.frame(cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 0,
                                       Counts = sum(rowSums(venn.df_sc[,-1]) == 0)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,2] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,3] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 0, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,4] == 1 & rowSums(venn.df_sc[,-1]) == 1)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 0,
                                       Counts = sum(venn.df_sc[,4] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 0, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,3] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 0, Sc.NH4 = 1, Sc.SO2 = 1,
                                       Counts = sum(venn.df_sc[,2] == 0 & rowSums(venn.df_sc[,-1]) == 2)),
                                 cbind(Sc.18C = 1, Sc.NH4 = 1, Sc.SO2 = 1,
                                       Counts = sum(rowSums(venn.df_sc[,-1]) == 3)))

venn.plot_sc <- cbind.data.frame(venn.plot_sc, 
                                 x = c(2.1, 0, -1.5, 1.5, -0.85, 0.85, 0, 0), 
                                 y = c(-2, 1.5, -0.5, -0.5, 0.5, 0.5, -1, 0))

venn.out3 <- data.frame(x = c(0, -0.75, 0.75), y = c(1, -0.5, -0.5), labels = c("18C", "NH4", "SO2"))
venn.out3$labels <- factor(venn.out3$labels, levels = c("18C", "NH4", "SO2"))

gg.venn_sc <- ggplot() +
  geom_circle(data = venn.out3, aes(x0 = x, y0 = y, r = 1.5, fill = labels), 
              alpha = 0.85, linewidth = 1, colour = "gray30") + 
  geom_text(data = venn.plot_sc, aes(x = x, y = y, label = Counts), size = 7) +
  coord_fixed() + theme_void() +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 13, color  = "black"),
        legend.title = element_text(size = 15, color = "black")) +
  scale_fill_manual(values = col.cond[-1]) +
  labs(fill = 'Condition')

gg.venn_sc

# Accumulated LFC - Histogram
hist_cond <- rbind(cbind(subset(ress_sc.18C, DEO == 1), Condition = "18C"),
                   cbind(subset(ress_sc.NH4, DEO == 1), Condition = "NH4"),
                   cbind(subset(ress_sc.SO2, DEO == 1), Condition = "SO2"))

gg.hist_sc <- ggplot(hist_cond, aes(x = abs(log2FoldChange), fill = Condition)) + 
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.75, color = "gray30") +
  theme_bw() + xlab("Absolute log2 Fold Change") + ylab("Number of DE genes") +
  scale_fill_manual(values = col.cond[-1]) +

  theme(axis.text.y = element_text(size = 15, color  = "black"), 
        axis.title.x = element_text(size = 17, color  = "black", margin = margin(r = 30)),
        axis.title.y = element_text(size = 17, color  = "black", margin = margin(t = 30)),
        legend.text = element_text(size = 13, color  = "black"),
        legend.title = element_text(size = 15, color  = "black"),
        axis.text.x = element_text(size = 15, color  = "black")) 

gg.hist_sc

gg.summary_sc <- plot_grid(gg.venn_sc, gg.hist_sc, nrow = 1, rel_widths = c(2, 3), labels = c("a", "b"), label_size = 20)
gg.summary_sc


## BIOLOGICAL ENRICHMENT - Gene Ontology BP
count.bias <- rowSums(ko_sc[,-1])
names(count.bias) <- ko_sc[,1]

DEG_sc.18C <- as.integer(ress_sc.18C$padj < 0.05 & abs(ress_sc.18C$log2FoldChange) >= 1)
names(DEG_sc.18C) <- ress_sc.18C$KEGG_ko

pwf_sc.18C <- nullp(DEgenes = DEG_sc.18C, bias.data = count.bias[names(DEG_sc.18C)])
go_sc.18C <- goseq(pwf_sc.18C, gene2cat = go_df, test.cats = c("GO:BP"))
go_sc.18C <- subset(go_sc.18C, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.18C <- go_sc.18C[p.adjust(go_sc.18C$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_sc.NH4 <- as.integer(ress_sc.NH4$padj < 0.05 & abs(ress_sc.NH4$log2FoldChange) >= 1)
names(DEG_sc.NH4) <- ress_sc.NH4$KEGG_ko

pwf_sc.NH4 <- nullp(DEgenes = DEG_sc.NH4, bias.data = count.bias[names(DEG_sc.NH4)])
go_sc.NH4 <- goseq(pwf_sc.NH4, gene2cat = go_df)
go_sc.NH4 <- subset(go_sc.NH4, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.NH4 <- go_sc.NH4[p.adjust(go_sc.NH4$over_represented_pvalue, method = "fdr") < 0.05,]

DEG_sc.SO2 <- as.integer(ress_sc.SO2$padj < 0.05 & abs(ress_sc.SO2$log2FoldChange) >= 1)
names(DEG_sc.SO2) <- ress_sc.SO2$KEGG_ko
DEG_sc.SO2[is.na(DEG_sc.SO2)] <- 0

pwf_sc.SO2 <- nullp(DEgenes = DEG_sc.SO2, bias.data = jitter(rep(1000, length(DEG_sc.SO2))))
go_sc.SO2 <- goseq(pwf_sc.SO2, gene2cat = go_df)
go_sc.SO2 <- subset(go_sc.SO2, ontology == "BP" & numDEInCat > 0)

enrichGO_sc.SO2 <- go_sc.SO2[p.adjust(go_sc.SO2$over_represented_pvalue, method = "fdr") < 0.05,] 

# Plot
go_sc <- rbind.data.frame(cbind(enrichGO_sc.18C, Comparison = "18C"),
                          cbind(enrichGO_sc.NH4, Comparison = "NH4"))

go_sc$comp.cat <- ifelse(duplicated(go_sc$term) | duplicated(go_sc$term, fromLast = TRUE), "Both", go_sc$Comparison)

go_sc$term <- factor(go_sc$term, levels = unique(go_sc$term[order(go_sc$comp.cat, go_sc$numDEInCat, decreasing = TRUE)]))

gg.go_sc <- ggplot(go_sc) +
  geom_bar(aes(x = term, y = numDEInCat, fill = Comparison), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#1e74eb", "#ebb249")) +
  coord_flip() +
  theme_bw() +
  theme(aspect.ratio = 0.62,
        legend.position = "bottom",
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 13, color = "black")) + 
  ylab("DE genes")  + xlab("")+
  labs(fill = 'Condition')
  

gg.go_sc

ggsave("C:/Users/Windows/OneDrive/Escritorio/Final_scripts/RNA/sc_summary.png", 
       plot_grid(gg.summary_sc, gg.go_sc, nrow = 2, labels = c("", "c"), label_size = 20, 
                 rel_widths = c(1,2)), 
       width = 14, height = 14, dpi = 300, bg = "white")


ko.n_df <- counts(dds_tot, normalized = TRUE) 
saveRDS(ko.n_df, "E:/Wineteractions/ko.n_df.rds")

#save.image("E:/R/Metatranscriptomics.RData")
