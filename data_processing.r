library(DESeq2)


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
  design = ~ batch + condition  # adjust this based on which data to process
)


# Normalize raw counts matrix
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# save normalized matrix as csv
save_path = '/Users/gdai/Desktop/BulkRNAseqAnalysis/8patient_subset/normalized_gene_cts_8patient_subset.csv'
write.table(normalized_counts, file=save_path, sep=",", col.names = NA, quote = FALSE)
