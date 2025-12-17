# Leave-one-out Cross Validation
# Tests the stability of cell line predictions by iteratively removing one cell and predicting on the remaining cells

# Load patient data
query <- readMM('./data/count_matrix_sparse.mtx')
rownames(query) <- readLines('./data/count_matrix_genes.tsv')
colnames(query) <- readLines('./data/count_matrix_barcodes.tsv')
queryMetadata <- read.csv('./data/metadata.csv', row.names = 1)

# Create pseudobulk profiles for all cells and cancer cells only
DC <- pbsapply(unique(queryMetadata$orig.ident), function(Donor){
  rowSums(query[,grepl(Donor, colnames(query))])
})
write.csv(DC, './results/allCellsPatientsPseudobulkProfiles.csv')

queryMetadata <- queryMetadata[grepl('Cancer',queryMetadata$celltype_minor),]
query <- query[,rownames(queryMetadata)]

DCC <- pbsapply(unique(queryMetadata$orig.ident), function(Donor){
  rowSums(query[,grepl(Donor, colnames(query))])
})
write.csv(DCC, './results/cancerCellsPatientPseudobulkProfiles.csv')

donorSubType <- queryMetadata$subtype
names(donorSubType) <- queryMetadata$orig.ident

# Select donor for LOOCV analysis
donor <- 'CID44971'
donorData <- query[,(queryMetadata$orig.ident %in% donor)]
cNames <- colnames(donorData)

cat(paste("Running LOOCV for donor", donor, "with", length(cNames), "cells\n"))

# Perform LOOCV
LOOCV <- pbsapply(colnames(donorData), function(C){
  # Remove one cell at a time
  donorData_loo <- donorData[, !cNames %in% C]
  donorMD <- data.frame(donor = rep(donor, ncol(donorData_loo)))
  
  # Map query to the integrated reference (clBRCA from first file)
  qmap <- mapQuery(donorData_loo, 
                   metadata_query = donorMD, 
                   ref_obj = clBRCA, 
                   do_umap = TRUE, 
                   do_normalize = TRUE)
  
  # Predict cell line identity using k-NN
  qmap <- knnPredict(qmap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  
  # Return prediction counts
  table(qmap$meta_data$cell_type_pred_knn)
})

write.csv(LOOCV, './results/ccLOOCV.csv')

cat("LOOCV analysis complete!\n")
cat(paste("Results saved to: ./results/ccLOOCV.csv\n\n"))

# Calculate LOOCV statistics
cvMean <- round((apply(LOOCV, 1, mean)/ncol(donorData)) * 100, 2)
cvSD <- round((apply(LOOCV, 1, sd)/ncol(donorData)) * 100, 2)

cat("LOOCV Summary Statistics:\n")
cat("Mean proportions per cell line (%):\n")
print(sort(cvMean, decreasing = TRUE))
cat("\nStandard deviations (%):\n")
print(sort(cvSD, decreasing = TRUE))

# Hierarchical Clustering Analysis

# Create pseudobulk profiles from integrated reference
cat("Creating pseudobulk profiles from integrated reference\n")
CL <- pbsapply(unique(X$cell_line), function(CT) {
  rowSums(GetAssayData(X, assay = "RNA", slot = "counts")[, X$cell_line %in% CT])
})
colnames(CL) <- unique(X$cell_line)
write.csv(CL, './results/cellLinePseudobulkProfiles.csv')

cat(paste("Created pseudobulk profiles for", ncol(CL), "cell lines\n"))

# Clustering with all patient cells
cat("\nClustering: Cell Lines + All Patient Cells\n")
geneList <- intersect(rownames(CL), rownames(DC))
COMBN <- data.frame(CL[geneList,], DC[geneList,])
COMBN <- as.matrix(COMBN)
COMBN <- log1p((t(t(COMBN)/colSums(COMBN)))*1e4)

spCor <- function(x){as.dist(cor(x, method = 'sp'))}
O_all <- pvclust(COMBN, method.dist = spCor, parallel = TRUE, nboot = 1000, method.hclust = 'complete')

png('./figures/S3.png', width = 4000, height = 1000, res = 300)
par(mar=c(1,4,1,1))
plot(O_all, print.pv = 'bp', print.num = FALSE, main = '', sub = '', xlab = '')
pvrect(O_all, alpha = 0.95)
dev.off()

cat("Saved: ./figures/S3.png\n")

# Clustering with cancer cells only
cat("\nClustering: Cell Lines + Cancer Cells Only\n")
geneList <- intersect(rownames(CL), rownames(DCC))
COMBN <- data.frame(CL[geneList,], DCC[geneList,])
COMBN <- as.matrix(COMBN)
COMBN <- log1p((t(t(COMBN)/colSums(COMBN)))*1e4)

O_cancer <- pvclust(COMBN, method.dist = spCor, parallel = TRUE, nboot = 1000, method.hclust = 'complete')

png('./figures/S4.png', width = 4000, height = 1000, res = 300)
par(mar=c(1,4,1,1))
plot(O_cancer, print.pv = 'bp', print.num = FALSE, main = '', sub = '', xlab = '')
pvrect(O_cancer, alpha = 0.95)
dev.off()

cat("Saved: ./figures/S4.png\n")