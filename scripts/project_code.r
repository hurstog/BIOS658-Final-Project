
#### Final Project Code ####

##### Load required packages 
library(GEOquery)
library(SummarizedExperiment)
library(ggplot2)
library(edgeR)
library(genefilter)
library(dplyr)
library(pheatmap)
library(pander)
library(ggrepel)
library(sva)
library(biomaRt)
library(tidyverse)
library(writexl)
library(enrichR) 
library(annotables)
library(clusterProfiler)
library(DOSE)
library(HGNChelper)
library(msigdbr)
library(nichenetr)
library(enrichplot)
library(ensembldb)
library(readxl)
library(ggVolcano)

##### Downloading and Preparing data
#Load Pheno data from GEO
GSE <- "GSE267064"
gse <- getGEO(GSE, GSEMatrix = TRUE, getGPL = FALSE)
expres_set <- gse[[1]]
pheno_data <- pData(expres_set)

#Load the raw data of gene expressions
raw_data <- read.csv("GSE267064_processed_counts.csv", row.names = 1)
count_data <- as.matrix(raw_data[, c(4:18)])

#We will remove the samples with 662 transfection since we are not interested in them
pheno_data <- pheno_data[-c(which(grepl("662", pheno_data[, 50]))), ]
count_data <- count_data[, -c(which(grepl("662", colnames(count_data))))]

#Create log2 transformed data for visualization 
log_data <- log2(count_data + 1)

#Create summarized experiment
#Change rownames of pheno data to match the samples in the count data (make sure they match)
rownames(pheno_data) <- pheno_data$title
rownames(pheno_data)
colnames(count_data)
#Must change the colnames of count data to match the row names of pheno data
columns <- colnames(count_data)
columns <- gsub("\\.", "-", columns)
colnames(count_data) <- columns 

#Now create the SE object
se <- SummarizedExperiment(
  assays = list(counts = count_data), colData = pheno_data)

#Create vectors for experimental conditions
treatment <- c("Control", "Treatment", "Control", "Treatment", "Control", 
               "Treatment", "Control", "Treatment", "Control", "Treatment",
               "Control", "Treatment")
genotype <- pheno_data[, 48]
#Combine them to one data frame
sample_annotations <- data.frame(
  Treatment = treatment,
  Genotype = genotype)

#Quality assessment:
#Library sizes:
#Calculate the sizes
lib_size <- matrix(data = NA, nrow = 12, ncol = 2)
lib_size[, 1] <- columns
colnames(lib_size) <- c("Sample", "Sum")
for (i in 1:12) {
  lib_size[i, 2] <- sum(count_data[, i])
}
#Create a dataframe and then visualize the library size
lib_size <- as.data.frame(lib_size)
lib_size$Sum <- as.numeric(lib_size$Sum)
ggplot(lib_size, aes(x = Sample, y = Sum)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample", y = "Sum of Counts", title = "Sum of Counts by Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
#Min and max library size
min(lib_size$Sum)
max(lib_size$Sum)
#Count Distribution: 
#Use the log2 transformed count data and create a boxplot of the count distribution of sample 
c_dist_data <- log_data
colnames(c_dist_data) <- c("Control1", "Treat1", "Control2", "Treat2",
                           "Control3", "Treat3", "Control4", "Treat4",
                           "Control5", "Treat5", "Control6", "Treat6")
boxplot(c_dist_data,
        las = 2,
        outline = FALSE,
        main = "Log2 Gene Counts per Sample",
        ylab = "log2(count + 1)")

#Sample correlation:
#Filter and normalize the gene counts
dge = DGEList(counts = count_data, genes = rownames(count_data))
keep <- filterByExpr(dge)  #Keep 12899 genes
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge = calcNormFactors((dge), method = "TMM")
dge = estimateDisp(dge)
#Get normalized CPMs 
mtx <- cpm(dge, log = TRUE, normalized.lib.sizes = TRUE)
rownames(mtx) <- dge@.Data[[3]][, "genes"]

#Plot the normalized gene counts to compare with raw data
cpm_data <- mtx
colnames(cpm_data) <- c("Control1", "Treat1", "Control2", "Treat2",
                           "Control3", "Treat3", "Control4", "Treat4",
                           "Control5", "Treat5", "Control6", "Treat6")
boxplot(cpm_data,
        las = 2,
        outline = FALSE,
        main = "Counts per million per Sample",
        ylab = "CPM")


#Filter out genes with low variation
mtx_to_plot <- varFilter(mtx, var.cutoff = 0.75)
#Now calculate correlation matrix and present as a heatmap
IAC <- mtx_to_plot %>% cor(. , use = "pairwise.complete.obs",  method = "pearson") 
plot_annotations <- data.frame(Sample = sample_annotations$Treatment)
rownames(plot_annotations) <- rownames(IAC)
pheatmap(IAC,annotation_col = plot_annotations, clustering_method = "ward.D")

#PCA:
#Perform PCA
pca <- prcomp(scale(t(mtx)))
summary(pca)$importance[,1:5] %>% pander
#Present a plot of the results
colorby <- "cells"
ggplot(
  data = data.frame(pca$x, cells = sample_annotations$Treatment, samples = sample_annotations$Treatment),
  aes(x = PC1, y = PC2, label = samples)
) +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle("PCA colored by treatment group") +
  geom_point(aes(color = cells), size = 3) +
  geom_text_repel(size = 3) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  scale_x_continuous(
    name = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 2), "%)")
  ) +
  scale_y_continuous(
    name = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 2), "%)")
  )
ggplot(
  data = data.frame(pca$x, cells = sample_annotations$Genotype, samples = sample_annotations$Genotype),
  aes(x = PC1, y = PC2, label = samples)
) +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle("PCA colored by sample genotype")) +
  geom_point(aes(color = cells), size = 3) +
  geom_text_repel(size = 3) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  scale_x_continuous(
    name = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 2), "%)")
  ) +
  scale_y_continuous(
    name = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 2), "%)")
  )

#Batch Effect Correction
#Batches are unknown, so we will use SVA
#Create model matrices
mod <- model.matrix(~ Treatment + Genotype, data = sample_annotations) #with biological covariates
mod0 <- model.matrix(~ 1, data = sample_annotations) #null model
# Estimate number of surrogate variables automatically (Leek method)
n.sv <- num.sv(mtx, mod, method = "leek")
mtx_norm <- cpm(edgeR.dgelist, log = FALSE, normalized.lib.sizes = TRUE) #Provide function with normalized but not log transformed data
svseq <- svaseq(dat = mtx_norm, mod = mod, mod0 = mod0, n.sv = n.sv)

#Perform hierarchical clustering on top 10% of genes (by variance)
gene_var <- apply(mtx, 1, var)
cutoff <- quantile(gene_var, 0.90)
top10_counts <- mtx[gene_var >= cutoff, ]
hclust.average <- hclust(dist(t(top10_counts)), method = "average")
hclust.average$height
plot(hclust.average, main = "Average Linkage Dendrogram",
     xlab = "Observation", ylab = "Distance")
#Note that we did not need to redo the PCA or hierarchical clustering since we found no significant batches

#Perform differential expression analysis using edgeR
#We use the dge from above and create a design matrix for the conditions of interest
design <- model.matrix(~ Genotype + Treatment, data = sample_annotations)
dge2 <- estimateDisp(dge, design)
fit <- glmFit(dge2, design)
qlf <- glmLRT(fit, coef= ncol(design))   
edgeR_full <- topTags(qlf, n = Inf)
dge_res <- edgeR_full$table

#Calculate the total number of DEGs
sum(dge_res$FDR < 0.05) #696

#Add gene symbols to the table
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
bm_annotations <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"), 
  filters = "ensembl_gene_id", 
  values = dge_res$genes, 
  mart = mart)
colnames(bm_annotations) <- c("Gene_Symbol", "genes")
dge_res <- left_join(dge_res, bm_annotations, by = "genes")
#Create box plots of the top 10 DEGs
#We first need to record which are the 10 most significant and create a dataframe to create the plots
dge_res <- dge_res[order(dge_res$FDR), ]
top10_genes <- dge_res$genes[1:10]
top10_expr <- mtx[top10_genes, ]
plot_df <- as.data.frame(top10_expr) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression")
group <- rep(c("Control", "Treatment"), times = 60)
plot_df <- cbind(plot_df, group)
ggplot(plot_df, aes(x = gene, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(title = "Top 10 Differentially Expressed Genes",
       y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
#Export results as CSV file
write.csv(dge_res, file = "DGE_results.csv", row.names = FALSE)
#We will also create a volcano plot of the results 
dge_res$FDR <- as.numeric(dge_res$FDR)
data <- add_regulate(data = dge_res, log2FC_name = "logFC",
                     fdr_name = "FDR", log2FC = 1, fdr = 0.3)

# plot
volplot <- gradual_volcano(data, x = "log2FoldChange", y = "padj",
                label = "Gene_Symbol", label_number = 10,
                fills = c("#53a4cf", "#aec5d1", "white", "#e2ae96", "#cc917c", "#a96d5d", "#91403e"),
                colors = "#404950",
                log2FC_cut = 0, FDR_cut = 1, add_line = FALSE, pointSizeRange = c(2, 2))
volplot + theme(legend.position = "right")

#Create bar plots for the functional analysis results
#functional_analysis.R contains the code used to generate the results used here
#Bar charts

file_path <- paste0(getwd(), "/GSEA2_proj_deg.xlsx")

# Get sheet names
sheet_names <- excel_sheets(file_path)

# Create an output folder for plots
output_dir <- getwd()

# Loop over each sheet
for (sheet in sheet_names) {
  
  # Read data
  df <- read_excel(file_path, sheet = sheet)
  
  # Expect columns: "Pathway" and "pvalue" (adjust if needed)
  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)
  
  # Clean up pathway names for readability
  df <- df %>%
    mutate(description = str_trunc(description, width = 60)) %>%
    mutate(description = factor(description, levels = rev(description)))
  
  # Plot
  p <- ggplot(df, aes(x = description, y = -1*log10(as.numeric(p.adjust)))) +
    geom_bar(stat = "identity", fill = "#2C7BB6", width = 0.7) +
    coord_flip() +
    labs(
      title = paste("Top 10 Pathways by P-value —", sheet),
      x = "Pathway",
      y = expression(-log[10]("adjusted p-value"))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(
    filename = paste0(output_dir, "/", sheet, "_top10_pathways.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}
