library(ggplot2)
library(dplyr)


# Read the data
gsea_report_tsv_path = '/Users/gdai/Desktop/BulkRNAseqAnalysis/8patient_subset/results/gsea_C5all_abnormal_upreg.tsv'
gsea_data <- read.delim(gsea_report_tsv_path, check.names = FALSE)  # disable header name conversion


# Clean and filter data
# Remove the empty column at the end, since GSEA tsv file has extra empty column
gsea_data <- gsea_data[, colSums(is.na(gsea_data)) != nrow(gsea_data)]

filtered_data <- gsea_data %>%
  filter(`FDR q-val` < 0.05) %>%
  head(20)


# Create the dotplot
ggplot(filtered_data, aes(x = NES, y = reorder(NAME, NES))) +
  geom_point(aes(color = `FDR q-val`), size = 4) +
  scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.05)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    title = "",
    x = "Normalized Enrichment Score",
    y = "",
    color = "FDR q-value"
  )
