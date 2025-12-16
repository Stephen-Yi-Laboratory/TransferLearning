# Integration with Quantile Normalization - Final Version
# Scale Kinker dataset to match Original dataset distribution

library(Seurat)
library(ggplot2)
library(Matrix)
library(symphony)
library(caret)
library(pbapply)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(ggforce)
library(patchwork)
library(reshape2)
library(data.table)
library(preprocessCore)  # For quantile normalization
library(pvclust) # For LOOCV

source('https://raw.githubusercontent.com/dosorio/utilities/master/data.frame2matrix.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/ggColors.R')

# ===== CONFIGURATION =====
N_TOP_GENES <- 3000

# Loading dictionary of gene symbols
ENSEMBL <- read.table('./data/hsaENSEMBL-GENES.txt', sep = '\t', header = TRUE)
ENSEMBL <- ENSEMBL[ENSEMBL$Gene_name != '',]
GENEID <- ENSEMBL$Gene_name
names(GENEID) <- ENSEMBL$Gene.stable.ID
ENSEMBL <- ENSEMBL[!duplicated(ENSEMBL$Gene.stable.ID), ]

# Function to extract cell line IDs
getCellLines <- function(X){
  cellLines <- colnames(X)
  cellLines <- unlist(lapply(strsplit(cellLines, '_'), function(X){X[1]}))
  return(cellLines)
}

# Quantile Normalization (Scale Kinker to match Original distribution)
quantile_normalize_to_reference <- function(reference_data, query_data) {
  cat("Applying quantile normalization to scale query to reference...\n")
  
  # Convert sparse to dense if needed
  if(inherits(reference_data, "sparseMatrix")) {
    reference_data <- as.matrix(reference_data)
  }
  if(inherits(query_data, "sparseMatrix")) {
    query_data <- as.matrix(query_data)
  }
  
  # For each gene, match Kinker's distribution to Original's distribution
  normalized_query <- query_data
  
  for(i in 1:nrow(query_data)) {
    if(i %% 500 == 0) cat(paste("Processing gene", i, "of", nrow(query_data), "\n"))
    
    # Get reference distribution for this gene
    ref_values <- reference_data[i, ]
    query_values <- query_data[i, ]
    
    # Rank-based transformation
    # Map Kinker's ranks to Original's quantiles
    query_ranks <- rank(query_values, ties.method = "average")
    query_quantiles <- query_ranks / length(query_values)
    
    # Use Original's quantile function to transform
    normalized_query[i, ] <- quantile(ref_values, probs = query_quantiles, type = 7)
  }
  
  rownames(normalized_query) <- rownames(query_data)
  colnames(normalized_query) <- colnames(query_data)
  
  return(normalized_query)
}

# Load original 32 cell line dataset

X_orig <- readRDS('./data/RAW.UMI.counts.BC.cell.lines.rds')
X_orig <- X_orig[rownames(X_orig) %in% ENSEMBL$Gene.stable.ID, , drop = FALSE]

gene_symbols <- as.vector(GENEID[rownames(X_orig)])
gene_symbols <- make.unique(gsub("_", "-", gene_symbols))
rownames(X_orig) <- gene_symbols

cat(paste("Original dataset:", nrow(X_orig), "genes x", ncol(X_orig), "cells\n"))

orig_cell_lines <- unique(getCellLines(X_orig))
cat(paste("\nOriginal dataset contains", length(orig_cell_lines), "cell lines:\n"))
print(sort(orig_cell_lines))

# Load Kinker et al. 2020 dataset

kinker_meta <- fread('./data/SCP542/metadata/Metadata.txt',
                     header = TRUE,
                     sep = '\t')

all_breast_lines <- unique(kinker_meta$Cell_line[grepl("_BREAST", kinker_meta$Cell_line)])
cat(paste("Found", length(all_breast_lines), "breast cancer cell lines in Kinker dataset:\n"))
print(all_breast_lines)

kinker_breast_meta <- kinker_meta[kinker_meta$Cell_line %in% all_breast_lines, ]
kinker_breast_cells <- kinker_breast_meta$NAME

header <- fread('./data/SCP542/other/UMIcount_data.txt', 
                header = TRUE,
                sep = '\t',
                nrows = 0)
all_cell_names <- colnames(header)[-1]

breast_cell_indices <- which(all_cell_names %in% kinker_breast_cells)
columns_to_load <- c(1, breast_cell_indices + 1)

cat("Loading UMI count data...\n")
kinker_breast_dt <- fread('./data/SCP542/other/UMIcount_data.txt',
                          header = TRUE,
                          sep = '\t',
                          select = columns_to_load,
                          check.names = FALSE)

gene_names <- kinker_breast_dt[[1]]
kinker_breast_dt <- kinker_breast_dt[, -1]
kinker_breast <- as.data.frame(kinker_breast_dt)
rownames(kinker_breast) <- gene_names

kinker_breast <- apply(kinker_breast, 2, as.numeric)
rownames(kinker_breast) <- gene_names

# Clean Kinker data
cat("\nCleaning Kinker data...\n")
na_genes <- rowSums(is.na(kinker_breast)) > 0
if(sum(na_genes) > 0) {
  cat(paste("  Removing", sum(na_genes), "genes with NA values\n"))
  kinker_breast <- kinker_breast[!na_genes, ]
}

gene_vars <- apply(kinker_breast, 1, var)
zero_var_genes <- gene_vars == 0 | is.na(gene_vars)
if(sum(zero_var_genes) > 0) {
  cat(paste("  Removing", sum(zero_var_genes), "genes with zero variance\n"))
  kinker_breast <- kinker_breast[!zero_var_genes, ]
}

all_zero_genes <- rowSums(kinker_breast) == 0
if(sum(all_zero_genes) > 0) {
  cat(paste("  Removing", sum(all_zero_genes), "all-zero genes\n"))
  kinker_breast <- kinker_breast[!all_zero_genes, ]
}

cat(paste("  After filtering:", nrow(kinker_breast), "genes x", ncol(kinker_breast), "cells\n"))

# Create proper cell names
cell_line_for_each_cell <- kinker_breast_meta$Cell_line[match(colnames(kinker_breast), kinker_breast_meta$NAME)]
cell_line_for_each_cell <- gsub("_BREAST", "", cell_line_for_each_cell)

new_names <- ave(cell_line_for_each_cell, cell_line_for_each_cell, 
                 FUN = function(x) paste0(x[1], "_", seq_along(x)))
colnames(kinker_breast) <- new_names

kinker_cell_lines <- unique(gsub("_.*", "", colnames(kinker_breast)))
cat(paste("\nKinker dataset contains", length(kinker_cell_lines), "cell lines:\n"))
print(sort(kinker_cell_lines))

# Identify overlapping cell lines

orig_cl_upper <- toupper(orig_cell_lines)
kinker_cl_upper <- toupper(kinker_cell_lines)

overlapping_lines <- intersect(orig_cl_upper, kinker_cl_upper)
cat(paste("Found", length(overlapping_lines), "overlapping cell lines:\n"))
print(overlapping_lines)

unique_to_original <- orig_cell_lines[!toupper(orig_cell_lines) %in% kinker_cl_upper]
unique_to_kinker <- kinker_cell_lines[!toupper(kinker_cell_lines) %in% orig_cl_upper]

cat(paste("\nCell lines unique to Original dataset:", length(unique_to_original), "\n"))
cat(paste("Cell lines unique to Kinker dataset:", length(unique_to_kinker), "\n"))
cat(paste("Overlapping cell lines (will merge):", length(overlapping_lines), "\n"))

kinker_breast_all <- kinker_breast

# Find common genes for intersection

common_genes <- intersect(rownames(X_orig), rownames(kinker_breast_all))
cat(paste("Common genes between datasets:", length(common_genes), "\n"))

X_orig_common <- X_orig[common_genes, ]
kinker_common <- kinker_breast_all[common_genes, ]

cat(paste("Original dataset subset:", nrow(X_orig_common), "genes x", ncol(X_orig_common), "cells\n"))
cat(paste("Kinker dataset subset:", nrow(kinker_common), "genes x", ncol(kinker_common), "cells\n"))

# Filter for highly expressed genes

if(!is.null(N_TOP_GENES)) {
  if(inherits(X_orig_common, "sparseMatrix")) {
    orig_mean <- Matrix::rowMeans(X_orig_common)
  } else {
    orig_mean <- rowMeans(X_orig_common)
  }
  
  if(inherits(kinker_common, "sparseMatrix")) {
    kinker_mean <- Matrix::rowMeans(kinker_common)
  } else {
    kinker_mean <- rowMeans(kinker_common)
  }
  
  combined_mean <- (orig_mean + kinker_mean) / 2
  n_genes_to_keep <- min(N_TOP_GENES, length(combined_mean))
  top_genes <- names(sort(combined_mean, decreasing = TRUE)[1:n_genes_to_keep])
  
  cat(paste("Filtering to top", length(top_genes), "most highly expressed genes\n"))
  
  X_orig_common <- X_orig_common[top_genes, ]
  kinker_common <- kinker_common[top_genes, ]
  
  cat(paste("\nAfter expression filtering:\n"))
  cat(paste("  Original dataset:", nrow(X_orig_common), "genes x", ncol(X_orig_common), "cells\n"))
  cat(paste("  Kinker dataset:", nrow(kinker_common), "genes x", ncol(kinker_common), "cells\n"))
}

# Quantile Normalization

kinker_scaled <- quantile_normalize_to_reference(X_orig_common, kinker_common)

# Check distribution match
cat("\nChecking distribution alignment:\n")
cat("Original dataset statistics:\n")
cat(paste("Mean of gene means:", round(mean(rowMeans(as.matrix(X_orig_common))), 2), "\n"))
cat(paste("SD of gene means:", round(sd(rowMeans(as.matrix(X_orig_common))), 2), "\n"))

cat("Kinker dataset statistics (after scaling):\n")
cat(paste("Mean of gene means:", round(mean(rowMeans(as.matrix(kinker_scaled))), 2), "\n"))
cat(paste("SD of gene means:", round(sd(rowMeans(as.matrix(kinker_scaled))), 2), "\n"))

# Dataset integration
X_combined <- cbind(X_orig_common, kinker_scaled)
cat(paste("Combined dataset:", nrow(X_combined), "genes x", ncol(X_combined), "cells\n"))

# Add metadata
dataset_source <- c(rep("Original", ncol(X_orig_common)), 
                    rep("Kinker", ncol(kinker_scaled)))
names(dataset_source) <- colnames(X_combined)

cell_lines_combined <- getCellLines(X_combined)

cat("\nCombined dataset summary:\n")
cat(paste("  Total unique cell lines:", length(unique(cell_lines_combined)), "\n"))
cat(paste("  From original only:", length(unique_to_original), "\n"))
cat(paste("  From Kinker only:", length(unique_to_kinker), "\n"))
cat(paste("  Overlapping (merged from both):", length(overlapping_lines), "\n"))

if(length(overlapping_lines) > 0) {
  cat("\nDetails for overlapping cell lines:\n")
  for(line in overlapping_lines) {
    orig_match <- orig_cell_lines[toupper(orig_cell_lines) == line]
    kinker_match <- kinker_cell_lines[toupper(kinker_cell_lines) == line]
    
    orig_cells <- sum(grepl(paste0("^", orig_match, "_"), colnames(X_orig_common)))
    kinker_cells <- sum(grepl(paste0("^", kinker_match, "_"), colnames(kinker_scaled)))
    
    cat(paste("  ", line, ": ", orig_cells, " cells from Original + ", 
              kinker_cells, " cells from Kinker (scaled) = ", 
              orig_cells + kinker_cells, " total\n", sep=""))
  }
}

# Create seurat object

X_combined_sparse <- as(X_combined, "sparseMatrix")

X <- CreateSeuratObject(X_combined_sparse, min.cells = 3, min.features = 200)

X$dataset_source <- dataset_source[colnames(X)]
X$cell_line <- getCellLines(X)

cat(paste("After Seurat filtering:", ncol(X), "cells,", nrow(X), "genes\n"))
cat(paste("Cell lines remaining:", length(unique(X$cell_line)), "\n"))

cat("\nApplying Seurat normalization (log)...\n")
X <- NormalizeData(X)

cat(paste("Final dataset:", ncol(X), "cells from", length(unique(X$cell_line)), "cell lines\n"))

# Summary
cat(paste("Dimensions:", nrow(X), "genes x", ncol(X), "cells\n"))
cat("\nCells by dataset source:\n")
print(table(X$dataset_source))
cat("\nCells by cell line:\n")
print(sort(table(X$cell_line), decreasing = TRUE))

# Save
cat("\nSaving integrated Seurat object...\n")
saveRDS(X, paste0('./results/integrated_breast_cancer_atlas', '.rds'))

metadata_df <- data.frame(
  cell_id = colnames(X),
  cell_line = X$cell_line,
  dataset_source = X$dataset_source,
  stringsAsFactors = FALSE
)
write.csv(metadata_df, './results/integrated_metadata.csv', row.names = FALSE)

# Clean up
rm(X_orig, kinker_breast, kinker_breast_dt, gene_names, header, all_cell_names,
   breast_cell_indices, columns_to_load, X_orig_common, kinker_common,
   kinker_scaled, X_combined, X_combined_sparse, kinker_breast_all)
gc()

# Show list of breast cancer types

clMD <- data.frame(CL = unique(getCellLines(X)), Type = NA)
clMD$Type <- c('Her2+', 'TNBC-A', 'TNBC-A', 'TNBC-A', 'Lum',
               'TNBC-A', 'TNBC-A', 'Her2+', 'Lum', 'TNBC-A',
               'TNBC-B', 'Her2+', 'TNBC-B', 'Lum', 'TNBC-B',
               'Her2+', 'TNBC-B', 'Lum', 'Lum', 'TNBC-A',
               'Lum', 'Her2+', 'TNBC-A', 'Her2+', 'TNBC-A',
               'TNBC-A', 'NI', 'TNBC-A', 'Lum', 'Lum',
               'TNBC-B', 'Her2+', 'TNBA-B', 'Lum', 'Her2+',
               'Her2+')
clMD <- clMD[order(clMD$CL),]
cellSubtype <- clMD$Type
names(cellSubtype) <- clMD$CL
clColor <- ggColors(nrow(clMD))
names(clColor) <- clMD$CL

# Cross Validation

steps <- seq(from = (ncol(X)-1000), to = 1000, by = -1000)
CM <- pbsapply(steps, function(nCells){
  set.seed(0)
  testData <- sample(seq_len(ncol(X)), nCells)
  trainData <- seq_len(ncol(X))[!seq_len(ncol(X)) %in% testData]
  
  # Use GetAssayData() to access the counts
  counts_matrix <- GetAssayData(X, assay = "RNA", slot = "counts")
  
  print(paste("Dimension of counts_matrix:", paste(dim(counts_matrix), collapse="x")))
  print(paste("Number of unique cell lines:", length(unique(getCellLines(counts_matrix[,trainData])))))
  
  metadata_ref <- data.frame(cellLine = getCellLines(counts_matrix[,trainData]))
  metadata_ref$cellLine <- factor(metadata_ref$cellLine)
  
  print(str(metadata_ref))
  
  clBRCA <- buildReference(
    exp_ref = counts_matrix[,trainData],
    metadata_ref = metadata_ref,
    do_umap = TRUE,
    verbose = TRUE,
    d = 50)
  
  testMap <- mapQuery(counts_matrix[,testData],
                      metadata_query = data.frame(cl = getCellLines(counts_matrix[,testData])),
                      ref_obj = clBRCA)
  testMap <- knnPredict(testMap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  CM <- confusionMatrix(data = as.factor(testMap$meta_data[,1]), reference = as.factor(testMap$meta_data[,2]), mode = 'everything')
  CM$overall
})

ACC <- data.frame(t(CM), steps)
ACC$steps <- (ncol(X)-ACC$steps)
ACC$Accuracy <- ACC$Accuracy * 100
ACC$AccuracyLower <- ACC$AccuracyLower * 100
ACC$AccuracyUpper <- ACC$AccuracyUpper * 100
write.csv(ACC, './results/F2B.csv')

P2 <- ggplot(ACC, aes(steps, Accuracy)) +
  stat_smooth(color = 'red') +
  geom_point() +
  geom_errorbar(aes(ymin = AccuracyLower, ymax = AccuracyUpper)) +
  xlab('Cells used for Training') +
  ylab('Symphony Accuracy (%)') +
  theme_bw()
print(P2)

# Create breast cancer atlas

clBRCA <- buildReference(
  exp_ref = GetAssayData(X, assay = "RNA", slot = "counts"),
  metadata_ref = data.frame(cellLine = getCellLines(GetAssayData(X, assay = "RNA", slot = "counts"))),
  do_umap = TRUE,
  verbose = TRUE,
  d = 50, 
  save_uwot_path = 'umapBRCA'
)

umapBRCA <- data.frame(clBRCA$umap$embedding, cl=clBRCA$meta_data$cellLine)
save(clBRCA, umapBRCA, file = './results/refBRCA.RData')

labelPos <- split(umapBRCA, umapBRCA$cl)
labelPos <- lapply(labelPos, function(S){
  mDistance <- mahalanobis(S[,1:2], colMeans(S[,1:2]), cov(S[,1:2]))
  S <- S[!mDistance %in% boxplot.stats(mDistance)$out,]
  as.data.frame(c(apply(S[,1:2],2,median), CL = S[1,3]))
})
labelPos <- t(do.call(cbind.data.frame, labelPos))
labelPos <- as.data.frame(labelPos)
rownames(labelPos) <- NULL
labelPos$UMAP1 <- as.numeric(labelPos$UMAP1)
labelPos$UMAP2 <- as.numeric(labelPos$UMAP2)
labelPos$cl <- factor(labelPos$CL, levels = clMD$CL)
umapBRCA <- umapBRCA[order(umapBRCA$cl, decreasing = TRUE),]
umapBRCA$cl <- factor(umapBRCA$cl, levels = clMD$CL)
write.csv('umapBRCA', './results/F2A.csv')

P1 <- ggplot(umapBRCA, aes(UMAP1, UMAP2, color = cl)) +
  geom_point(cex = 0.01) +
  theme_bw() +
  theme(legend.position = 'None') +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  #scale_color_viridis_d() +
  geom_text_repel(aes(UMAP1, UMAP2, label = CL),
                  labelPos,
                  min.segment.length = 0,
                  nudge_y = 1.5,
                  bg.color = 'white')
print(P1)

# Benchmarking

# ===== TESTING MCF7 =====
cat("\n===== Testing MCF7 =====\n")
MCF7 <- read.csv('./data/MCF7.csv')
MCF7 <- as.matrix(MCF7)

mcf7Map <- mapQuery(MCF7,
                    metadata_query = data.frame(rep('MCF7', ncol(MCF7))),
                    ref_obj = clBRCA, do_umap = TRUE)

mcf7Map <- knnPredict(mcf7Map, clBRCA, clBRCA$meta_data$cellLine, k = 5)
mcf7CM <- confusionMatrix(as.factor(mcf7Map$meta_data[,1]), as.factor(mcf7Map$meta_data[,2]))

mDistance <- mahalanobis(mcf7Map$umap, colMeans(mcf7Map$umap), cov(mcf7Map$umap))
plotData <- rbind(data.frame(clBRCA$umap$embedding, cl='Ref'), data.frame(mcf7Map$umap, cl='MCF7'))
plotData$ct <- c(clBRCA$meta_data$cellLine, ifelse(!mDistance %in% boxplot.stats(mDistance)$out, 'Q', 'O'))
write.csv(plotData, './results/F2C_kinker.csv')

P3 <- ggplot(plotData, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01,
             color = ifelse(plotData$ct %in% c('Q', 'O'), rgb(1,0,0,0.01), 'gray75'),
             alpha = 1) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_mark_ellipse(aes(filter = ct == 'Q', color = 'red'),expand = unit(2,'mm')) +
  annotate(x = median(mcf7Map$umap[,1]) - 3, 
           y = max(mcf7Map$umap[,2]) - 6,
           geom = 'text', 
           label = paste0(round(mcf7CM$overall[1]*100,1), '% MCF7'),
           color = 'red', 
           fontface = 'bold', 
           size = 3) +
  xlab('UMAP 1') +
  ylab('UMAP 2')
P3 <- P3 + labs(title = 'Wild-Type MCF7', 
                subtitle = parse(text = 'italic(n)==14372~Cells')) +
  theme(plot.title = element_text(face = 2))
print(P3)

# MCF7 Accuracy Metrics
cat("\nMCF7 Classification Accuracy:\n")
cat(paste0("Accuracy: ", round(mcf7CM$overall[1]*100, 2), "%\n"))
cat(paste0("95% CI: [", round(mcf7CM$overall[3]*100, 2), "%, ", 
           round(mcf7CM$overall[4]*100, 2), "%]\n"))
cat(paste0("Kappa: ", round(mcf7CM$overall[2], 3), "\n"))

cat("\nPrediction breakdown:\n")
pred_table <- table(mcf7Map$meta_data$cell_type_pred_knn)
print(sort(pred_table[pred_table > 0], decreasing = TRUE))

# ===== TESTING T47D =====
cat("\n===== Testing T47D =====\n")
td47dData <- read.csv('./data/GSM4285803_scRNA_RawCounts.csv.gz', row.names = 1)
td47dData <- t(td47dData)
td47dMetaData <- read.csv('./data/GSM4285803_scRNA_metaInfo.csv.gz')
td47dData <- td47dData[,td47dMetaData$X[grepl('T47D KO', td47dMetaData$CellType)]]
td47dData <- as.matrix(td47dData)

td47Map <- mapQuery(exp_query = td47dData,
                    metadata_query = data.frame(rep('T47D', ncol(td47dData))),
                    ref_obj = clBRCA, do_umap = TRUE)
td47Map <- knnPredict(td47Map, clBRCA, clBRCA$meta_data$cellLine, k = 5)
td47CM <- confusionMatrix(as.factor(td47Map$meta_data[,1]), as.factor(td47Map$meta_data[,2]))

td47_umap <- td47Map$umap
mDistance <- mahalanobis(td47_umap, colMeans(td47_umap), cov(td47_umap))
outlier_threshold <- quantile(mDistance, 0.95)
td47_inliers <- mDistance <= outlier_threshold

plotData <- rbind(data.frame(clBRCA$umap$embedding, cl='Ref'), 
                  data.frame(td47Map$umap, cl='T47D'))
plotData$ct <- c(clBRCA$meta_data$cellLine, 
                 ifelse(td47_inliers, 'Q', 'O'))
write.csv(plotData, './results/F2D_kinker.csv')

P4 <- ggplot(plotData, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01, color = ifelse(plotData$ct %in% c('Q', 'O'), 'red', 'gray75'), alpha = 1) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_mark_ellipse(aes(filter = ct == 'Q', color = 'red'), expand = unit(0.5,'mm')) +
  annotate(x = median(td47Map$umap[td47_inliers,1]), 
           y = max(td47Map$umap[td47_inliers,2]) + 2,  # Just above the main cluster
           geom = 'text', 
           label = paste0(round(td47CM$overall[1]*100, 1), '% T47D'),
           color = 'red', 
           fontface = 'bold',
           size = 3.5) +
  xlab('UMAP 1') +
  ylab('UMAP 2')

P4 <- P4 + labs(title = 'T47D with CDH1 knockout',
                subtitle = parse(text = 'italic(n)==491~Cells')) +
  theme(plot.title = element_text(face = 2))
print(P4)

# T47D Accuracy Metrics
cat("\nT47D Classification Accuracy:\n")
cat(paste0("Accuracy: ", round(td47CM$overall[1]*100, 2), "%\n"))
cat(paste0("95% CI: [", round(td47CM$overall[3]*100, 2), "%, ", 
           round(td47CM$overall[4]*100, 2), "%]\n"))
cat(paste0("Kappa: ", round(td47CM$overall[2], 3), "\n"))

cat("\nPrediction breakdown:\n")
pred_table_t47d <- table(td47Map$meta_data$cell_type_pred_knn)
print(sort(pred_table_t47d[pred_table_t47d > 0], decreasing = TRUE))

# ===== TESTING BT474 =====
cat("\n===== Testing BT474 =====\n")
bt474Data <- read.csv('./data/GSE150949_pooled_watermelon.count.matrix.csv.gz', row.names = 1)
bt474MetaData <- read.csv('./data/GSE150949_pooled_watermelon.metadata.matrix.csv.gz')
bt474MetaData <- bt474MetaData[grepl('BT474', bt474MetaData$cell_line),]
bt474Data <- bt474Data[,gsub('-','.',bt474MetaData$cell)]
bt474Data <- as.matrix(bt474Data)

bt474Map <- mapQuery(exp_query = bt474Data,
                     metadata_query = data.frame(rep('BT474', ncol(bt474Data))),
                     ref_obj = clBRCA, do_umap = TRUE)
bt474Map <- knnPredict(bt474Map, clBRCA, clBRCA$meta_data$cellLine, k = 5)
bt474CM <- confusionMatrix(as.factor(bt474Map$meta_data[,1]), as.factor(bt474Map$meta_data[,2]))

# More aggressive outlier filtering - remove cells far from the main cluster
bt474_umap <- bt474Map$umap
mDistance <- mahalanobis(bt474_umap, colMeans(bt474_umap), cov(bt474_umap))

# Use a stricter threshold - only keep cells within 2 standard deviations
outlier_threshold <- quantile(mDistance, 0.76)  # Keep only the densest 90%
bt474_inliers <- mDistance <= outlier_threshold

# Create plotData
plotData <- rbind(data.frame(clBRCA$umap$embedding, cl='Ref'), 
                  data.frame(bt474Map$umap, cl='BT474'))
plotData$ct <- c(clBRCA$meta_data$cellLine, 
                 ifelse(bt474_inliers, 'Q', 'O'))

write.csv(plotData, './results/F2E_integration.csv')

P5 <- ggplot(plotData, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01,
             color = ifelse(plotData$ct %in% c('Q', 'O'), 'red', 'gray75'),
             alpha = 1) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_mark_ellipse(aes(filter = ct == 'Q', color = 'red'), 
                    expand = unit(0.5,'mm')) +
  annotate(x = median(bt474Map$umap[bt474_inliers,1] - 5), 
           y = max(bt474Map$umap[bt474_inliers,2]) - 3,
           geom = 'text', 
           label = paste0(round(bt474CM$overall[1]*100,1), '% BT474'),
           color = 'red', 
           fontface = 'bold',
           size = 3) +
  xlab('UMAP 1') +
  ylab('UMAP 2')
P5 <- P5 + labs(title = 'BT474 + Lapatinib', 
                subtitle = parse(text = 'italic(n)==131~Cells')) +
  theme(plot.title = element_text(face = 2))
print(P5)

# BT474 Accuracy Metrics
cat("\nBT474 Classification Accuracy:\n")
cat(paste0("Accuracy: ", round(bt474CM$overall[1]*100, 2), "%\n"))
cat(paste0("95% CI: [", round(bt474CM$overall[3]*100, 2), "%, ", 
           round(bt474CM$overall[4]*100, 2), "%]\n"))
cat(paste0("Kappa: ", round(bt474CM$overall[2], 3), "\n"))

cat("\nPrediction breakdown:\n")
pred_table_bt474 <- table(bt474Map$meta_data$cell_type_pred_knn)
print(sort(pred_table_bt474[pred_table_bt474 > 0], decreasing = TRUE))

# Create Combined Figure
png('./figures/F2_kinker_integrated.png', width = 3500, height = 1750, res = 300)
pLayout <- '
AABC
AADE'
P1 + P2 + P3 + P4 + P5 + plot_layout(design = pLayout) + 
  plot_annotation(tag_levels = 'A', 
                  title = 'Reference Validation using Integrated Datasets',
                  theme = theme(plot.tag = element_text(face = 2),
                                plot.title = element_text(face = 2, size = 16)))
dev.off()

# Patient data
query <- readMM('./data/count_matrix_sparse.mtx')
rownames(query) <- readLines('./data/count_matrix_genes.tsv')
colnames(query) <- readLines('./data/count_matrix_barcodes.tsv')
queryMetadata <- read.csv('./data/metadata.csv', row.names = 1)
queryMetadata <- queryMetadata[grepl('Cancer',queryMetadata$celltype_minor),]
query <- query[,rownames(queryMetadata)]
donorSubType <- queryMetadata$subtype
names(donorSubType) <- queryMetadata$orig.ident

clProportion <- lapply(unique(queryMetadata$orig.ident), function(donor){
  donorData <- query[,(queryMetadata$orig.ident %in% donor)]
  donorMD <- data.frame(donor = rep(donor, ncol(donorData)))
  qmap <- mapQuery(donorData, metadata_query = donorMD, ref_obj = clBRCA, do_umap = TRUE, do_normalize = TRUE)
  qmap <- knnPredict(qmap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  qc <- qmap$meta_data
  qc <- table(qc$cell_type_pred_knn)/nrow(qc)
  qc <- data.frame(donor = donor, qc)
  colnames(qc) <- c('donor', 'cellLine', 'proportion')
  qc
})

clProportion <- do.call(rbind.data.frame, clProportion)
clProportion$subtype <- donorSubType[clProportion$donor]
clProportion$cellsubtype <- cellSubtype[as.vector(clProportion$cellLine)]
clProportion$cellsubtype[grepl('TNBC', clProportion$cellsubtype)] <- 'TNBC'
clProportion$cellsubtype[grepl('Lum', clProportion$cellsubtype)] <- 'ER+'
clProportion$cellLine <- factor(clProportion$cellLine, levels = clMD$CL)
write.csv(clProportion, './results/F3C_integrated.csv')

clProportion <- read.csv('./results/F3C_integrated.csv', row.names = 1)
O <- round(acast(data = clProportion, formula = donor~cellLine, value.var = 'proportion') * 100,2)
O <- O[,clMD$CL]
col_fun = colorRamp2(c(0, 50), c("gray99", "red"))

png('./figures/S1.png', width = 4000, height = 2500, res = 300)
Heatmap(O, col = col_fun,
        column_split = clMD$Type,name = '%',
        row_split = donorSubType[rownames(O)],
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.1f", O[i, j]), x, y, gp = gpar(fontsize = 10))})
dev.off()

chm <- acast(clProportion, subtype ~ cellLine, value.var = 'proportion', fun.aggregate = mean)
chm <- chm/rowSums(chm)
chm <- as.matrix(chm)
write.csv(chm, './results/F4.csv')

chm <- read.csv('./results/F4.csv', row.names = 1)
chm <- chm * 100
col_fun = colorRamp2(c(0, 36.5), c("gray99", "red"))
#chm <- scale(chm)
chm <- round(chm,1)
png('./figures/F4.png', width = 3500, height = 750, res = 300)
Heatmap(chm,
        column_split = cellSubtype[colnames(chm)],
        col = col_fun, name = '%',
        show_row_dend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.1f", chm[i, j]), x, y, gp = gpar(fontsize = 10))})
dev.off()

pInfo <- data.frame(donor = clProportion$donor, subtype = clProportion$subtype)
pInfo <- unique(pInfo)
pInfo <- pInfo[order(pInfo$subtype),]
clProportion$donor <- factor(clProportion$donor,levels = pInfo$donor)
clProportion$proportion <- clProportion$proportion * 100
clProportion <- clProportion[order(clProportion$cellLine),]
clProportion$cellLine <- factor(clProportion$cellLine, levels = clMD$CL)

P7 <- ggplot(clProportion, aes(x = donor, y = proportion , fill = cellLine)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(face = 2)) +
  xlab('Donor') +
  ylab('%') +
  labs(fill = 'Reference\nCell Line')
print(P7)

# Donor specific profile
donor <- 'CID44971'
donorData <- query[,(queryMetadata$orig.ident %in% donor)]
donorMD <- data.frame(donor = rep(donor, ncol(donorData)))
qmap <- mapQuery(donorData, metadata_query = donorMD, ref_obj = clBRCA, do_umap = TRUE, do_normalize = TRUE)
qmap <- knnPredict(qmap, clBRCA, clBRCA$meta_data$cellLine, k = 5)

qc <- qmap$meta_data

A <- data.frame(clBRCA$umap$embedding, cl = NA)
B <- data.frame(qmap$umap, cl = qmap$meta_data$cell_type_pred_knn)
qmap <- rbind(A,B)
qmap$cl <- factor(qmap$cl, levels = clMD$CL)
lInfo <- round((table(qmap$cl)/sum(table(qmap$cl)))*100,2)
lInfo <- paste0(lInfo, '% ', names(lInfo))
lInfo <- gsub('0% ', '', lInfo)
names(lInfo) <- names(table(qmap$cl))
labelPos$pct <- lInfo[as.vector(labelPos$cl)]
labelPos$pct[!grepl('%',labelPos$pct)] <- NA
write.csv(qmap, './results/F3A.csv')

P6 <- ggplot(qmap, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01, color = ifelse(is.na(qmap$cl), 'gray75', rgb(1,0,0,1))) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_text_repel(aes(UMAP1, UMAP2, label = pct),
                  labelPos,
                  min.segment.length = 0,
                  nudge_y = 1.7,
                  bg.color = 'white',
                  col = 'black', size = 4) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  labs(title = donor,
       subtitle = parse(text = paste0('italic(n)==', ncol(donorData), '~Cancer~Cells'))) +
  theme(plot.title = element_text(face = 2))
P6

source('./S2-LOOCV-Integration.R')
LOOCV <- read.csv('./results/ccLOOCV_integrated.csv', row.names = 1)
cvMean <- round((apply(LOOCV,1,mean)/ncol(donorData)) * 100,2)
cvSD <- round((apply(LOOCV,1,sd)/ncol(donorData)) * 100,2)
RMSE <- sqrt(mean((round(sort(table(qc$cell_type_pred_knn)/nrow(qc)*100, decreasing = TRUE),2)-sort(cvMean, decreasing = TRUE))^2))
CTest <- cor.test(table(qc$cell_type_pred_knn)/nrow(qc)*100, cvMean, method = 'sp', continuity = TRUE)
CTest$p.value
DF <- data.frame(CL = rownames(LOOCV), M = cvMean, LB = cvMean-cvSD, UB = cvMean+cvSD)
DF <- DF[order(DF$M),]
DF$CL <- factor(DF$CL, levels = DF$CL)
DF$CL2 <- factor(DF$CL, levels = clMD$CL)
DF <- DF[DF$M > 0,]
write.csv(DF, './results/F3B.csv')

P8 <- ggplot(DF, aes(M,CL)) +
  geom_bar(stat = 'identity', fill = clColor[DF$CL2]) +
  geom_errorbarh(mapping = aes(xmin = LB, xmax = UB), height = .25, size = 0.5) +
  theme_bw() +
  theme(legend.position = 'None') +
  xlab('%') +
  ylab('Cell Line')
print(P8)


png('./figures/F3.png', width = 3500, height = 2000, res = 300)
pLayout <- '
AAACCC
AAACCC
BBBCCC
'
P6 + P8 + P7 + plot_layout(design = pLayout) + plot_annotation(tag_levels = 'A', theme = theme(plot.tag = element_text(face = 2)))
dev.off()