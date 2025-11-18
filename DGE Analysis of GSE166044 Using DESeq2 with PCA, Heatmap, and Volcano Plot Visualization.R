# -------------------------------------------------------------------
# Project: DGE Analysis of GSE166044 Using DESeq2 with PCA, 
#          Heatmap, and Volcano Plot Visualization
#
# Dataset: GSE166044
# Title:   Whole transcriptome analysis of breast tissues prior to 
#          breast cancer diagnosis [RNA-seq]
# -------------------------------------------------------------------

# Define the GEO Accession ID
geo_id <- "GSE166044"

# -------------------------------------------------------------------
#  Step 1: Package Installation
# -------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2", "pheatmap", "EnhancedVolcano", "ggplot2"), update = FALSE)

# -------------------------------------------------------------------
#  Step 2: Load Libraries
# -------------------------------------------------------------------
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(utils) # For untar()

# -------------------------------------------------------------------
# Step 2b: Define Output Path and Create Directory
# -------------------------------------------------------------------
# This is your requested default path.
output_dir <- "D:/DOWNLOADS"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------------------------------------------------
# Step 3: Fetch Metadata & Download Count Matrix
# -------------------------------------------------------------------

# 3a. Fetch Metadata (pdata) from GEO
gse <- getGEO(geo_id, GSEMatrix = TRUE)
gse <- gse[[1]]
pdata <- pData(gse) # pdata = Phenotype (metadata)

cat(paste("Metadata (pdata) dimensions for", geo_id, ":\n"))
print(dim(pdata)) # Should be 30 samples


# 3b. Download the Supplementary Count Matrix file
# This downloads "GSE166044_RAW.tar" *into* a new folder named "GSE166044"
getGEOSuppFiles(geo_id)


# -------------------------------------------------------------------
# Step 3c: (Corrected) Process TAR, Filter, and Convert to Numeric
# -------------------------------------------------------------------
cat("\nProcessing .tar file...\n")

count_dir <- geo_id # "GSE166044"
tar_file <- file.path(count_dir, paste0(geo_id, "_RAW.tar")) # "GSE166044/GSE166044_RAW.tar"

if (!file.exists(tar_file)) {
  stop(paste("ERROR: Tar file not found at", tar_file))
}

# Extract the .tar file *inside* its own directory
untar(tar_file, exdir = count_dir)

# Get the list of all the individual count files
count_files <- list.files(path = count_dir, 
                          pattern = "*.txt.gz", 
                          full.names = TRUE) 

cat(paste("Found", length(count_files), "count files.\n"))

# --- Start of v4 Fix ---

# Read first file to create the filter mask
first_file_data <- read.delim(count_files[1], header = FALSE, col.names = c("GeneID", "Count"))

# Create a 'mask' to find all rows that do NOT start with "__"
valid_rows_mask <- !grepl("^__", first_file_data$GeneID)

# Get the clean list of gene names
gene_names <- first_file_data$GeneID[valid_rows_mask]
cat(paste("Found", length(gene_names), "valid genes, filtering out summary lines.\n"))

# Loop through all files
all_counts_list <- lapply(count_files, function(file) {
  # Read the data. 'Count' column will be forced to CHARACTER type 
  data <- read.delim(file, header = FALSE, col.names = c("GeneID", "Count"))
  
  # 1. Filter to get the "good" rows (which are still text)
  valid_counts_character <- data$Count[valid_rows_mask]
  
  # 2. (THE FIX) Explicitly convert the "good" text back to numbers
  valid_counts_numeric <- as.numeric(as.character(valid_counts_character))
  
  return(valid_counts_numeric)
})
# --- End of Fix ---

# Combine the list of (now purely numeric) counts into one matrix
exprSet <- do.call(cbind, all_counts_list)

# Assign Gene Names as rownames
rownames(exprSet) <- gene_names

# Clean up sample names 
basenames <- basename(count_files)
sample_names <- gsub("_.*", "", basenames) 
colnames(exprSet) <- sample_names

# Ensure counts are integers (this will now work)
exprSet <- round(exprSet)

cat(paste("\nExpression matrix built successfully for", geo_id, ":\n"))
print(dim(exprSet))


# -------------------------------------------------------------------
# Step 4: (Corrected) Define Sample Condition
# -------------------------------------------------------------------
# This uses the correct 'characteristics_ch1' column
pdata$condition <- ifelse(grepl("Susceptible", pdata$characteristics_ch1, ignore.case = TRUE), 
                          "Susceptible", 
                          "Normal")

pdata$condition <- as.factor(pdata$condition)
cat("\nSample Group Distribution (Must show 2 groups):\n")
print(table(pdata$condition))

# -------------------------------------------------------------------
# Step 5: (Corrected) Prepare DESeq2 Dataset (Remove NAs)
# -------------------------------------------------------------------

# 5a. Remove rows (genes) with NA values
cat(paste("Original gene count (before NA removal):", nrow(exprSet), "\n"))

# na.omit() removes any row containing at least one NA
exprSet_clean <- na.omit(exprSet)

cat(paste("New gene count (after NA removal):", nrow(exprSet_clean), "\n"))

# 5b. Align metadata (pdata) to the *clean* count matrix
pdata_subset <- pdata[colnames(exprSet_clean), ]

# 5c. CRITICAL CHECK: Ensure order matches
stopifnot(all(colnames(exprSet_clean) == rownames(pdata_subset)))

# 5d. Create the DESeqDataSet object using the *clean* matrix
dds <- DESeqDataSetFromMatrix(countData = exprSet_clean,  # <-- Using exprSet_clean
                              colData = pdata_subset,
                              design = ~ condition)

# 5e. Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]
cat(paste("\nDESeqDataSet object created successfully for", geo_id, ".\n"))


# -------------------------------------------------------------------
# Step 6: Run DESeq2 Pipeline
# -------------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
summary(res)

# -------------------------------------------------------------------
# Step 7: Variance Stabilizing Transformation (VST)
# -------------------------------------------------------------------
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
cat("\nVST data dimensions:\n")
print(dim(assay(vsd)))

# -------------------------------------------------------------------
# Step 8: PCA Plot (Dimensional Reduction) & Save
# -------------------------------------------------------------------
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggtitle(paste("PCA:", geo_id, "- Susceptible vs Normal"))

pca_filename <- file.path(output_dir, paste0(geo_id, "_PCA_plot.png"))
ggsave(filename = pca_filename, plot = pca_plot)
cat(paste("\nPCA plot saved to:", pca_filename, "\n"))

print(pca_plot)

# -------------------------------------------------------------------
# Step 9: Heatmap of Top 50 Variable Genes & Save
# -------------------------------------------------------------------
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

heatmap_filename <- file.path(output_dir, paste0(geo_id, "_Heatmap_Top50.png"))

pheatmap(assay(vsd)[topVarGenes, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "condition", drop = FALSE]),
         main = paste(geo_id, ": Top 50 Variable Genes"),
         filename = heatmap_filename)

cat(paste("Heatmap saved to:", heatmap_filename, "\n"))

# -------------------------------------------------------------------
# Step 10: Volcano Plot & Save
# -------------------------------------------------------------------
volcano_plot <- EnhancedVolcano(res,
                                lab = rownames(res),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                pCutoff = 0.05,
                                FCcutoff = 1.5,
                                title = paste(geo_id, ': Volcano Plot - Susceptible vs Normal'),
                                subtitle = 'Differential Gene Expression',
                                legendLabels = c('NS', 'Log2FC', 'p-value', 'p-value & Log2FC'))

volcano_filename <- file.path(output_dir, paste0(geo_id, "_Volcano_plot.png"))
ggsave(filename = volcano_filename, plot = volcano_plot)
cat(paste("Volcano plot saved to:", volcano_filename, "\n"))

print(volcano_plot)

# Craft a prompt summarizing your findings and asking for interpretation
results_summary_prompt <- paste0(
  "The DESeq2 analysis on GSE166044 found 1951 up-regulated and 2028 down-regulated genes (p-adjusted < 0.1) when comparing Susceptible to Normal tissue. Given this context, what are the next steps for pathway analysis, and what kind of biological processes might these genes be involved in?"
)

# Run the query
gemini(results_summary_prompt)
