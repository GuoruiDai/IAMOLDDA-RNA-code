library(DESeq2)
library(biomaRt)
library(pheatmap)


# Specify input path
count_matrix_path = '/Users/gdai/Desktop/BulkRNAseqAnalysis/8patient_subset/raw_gene_cts_8patient_subset.csv'
condition_mapping_path = '/Users/gdai/Desktop/BulkRNAseqAnalysis/8patient_subset/patient_metadata_8patient_subset.csv'

# Read data
counts_data <- read.csv(count_matrix_path, sep=",", row.names = 1)
meta_data <- read.csv(condition_mapping_path, sep=",", row.names = 1)

# Create a DESeq2 obj
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = meta_data,
  design = ~ condition #+ batch
)
# Specify comparing from abnormal group to the normal group
dds$condition <- relevel(dds$condition, ref = 'Normal DLco')


# Run the DESeq function
dds <- DESeq(dds)
# Retrieve the results
res <- results(dds)
# Re-order based on p value
res <- res[order(res$padj),] 
res_copy = res


# -------Export DEG analysis results---------
# convert ensemble id to gene symbol in the result
ensembl_ids <- rownames(res_copy)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = ensembl_ids,
  mart = ensembl
)
gene_symbol_map <- setNames(gene_symbols$hgnc_symbol, gene_symbols$ensembl_gene_id)
# keep the original Ensembl ID if not found
res_genes_mapped <- ifelse(ensembl_ids %in% names(gene_symbol_map), gene_symbol_map[ensembl_ids], ensembl_ids)
res_genes_mapped[res_genes_mapped == ""] <- ensembl_ids[res_genes_mapped == ""]
rownames(res_copy) <- res_genes_mapped

result_export_path = '/Users/gdai/Desktop/BulkRNAseqAnalysis/8patient_subset/DEG_results_8patient_subset.csv'
write.csv(as.data.frame(res_copy), file=result_export_path) # save as csv
# -------------------------------------------------




# Filter out NA values in padj, then filter significant genes and order by log fold change
significant_genes <- res_copy[!is.na(res_copy$padj) & res_copy$padj < 0.05, ]
significant_genes <- significant_genes[order(significant_genes$log2FoldChange, decreasing = TRUE), ]

# Convert gene IDs for all significant genes first
ensembl_ids <- rownames(significant_genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = ensembl_ids,
  mart = ensembl
)
gene_symbol_map <- setNames(gene_symbols$hgnc_symbol, gene_symbols$ensembl_gene_id)

# Keep the original Ensembl ID if not found
res_genes_mapped <- ifelse(ensembl_ids %in% names(gene_symbol_map), 
                           gene_symbol_map[ensembl_ids], 
                           ensembl_ids)
res_genes_mapped[res_genes_mapped == ""] <- ensembl_ids[res_genes_mapped == ""]

# Filter out unmapped genes (those that are still Ensembl IDs)
mapped_genes <- res_genes_mapped != ensembl_ids  # Create logical vector of successfully mapped genes
significant_genes_filtered <- significant_genes[mapped_genes, ]  # Filter the significant genes

# Now select top k genes from the filtered set
top_gene_indices <- head(rownames(significant_genes_filtered), n=20)
res_top_genes <- assay(vst(dds))[top_gene_indices, ]
rownames(res_top_genes) <- res_genes_mapped[mapped_genes][1:20]  # Update with mapped names

# Re-order patient columns based on condition
ordered_cols <- colnames(res_top_genes)[order(meta_data$condition)]

# Create annotation dataframe with the same ordering
annotation_col <- data.frame(
  Condition = meta_data$condition[order(meta_data$condition)],
  row.names = ordered_cols
)

# Create the heatmap
pheatmap(res_top_genes[, ordered_cols],
         scale = "row",
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         main = 'Up-regulated DE genes in the DLco<LLN patient group',
         annotation_col = annotation_col,
         annotation_colors = list(
           Condition = c(`Normal DLco` = "#90EE90",
                         `Abnormal DLco` = "#F08080")
         ),
         gaps_col = sum(meta_data$condition == "Abnormal DLco"))

