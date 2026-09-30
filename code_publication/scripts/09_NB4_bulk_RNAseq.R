# ============================================================================
# Script: 09_NB4_bulk_RNAseq.R
# Purpose: Reproducible differential-expression analysis of the three NB4
#          control and three CD69-overexpression bulk RNA-seq libraries in
#          GSE313879. This script does not infer sample groups from filenames;
#          an explicit sample metadata file is required.
# Inputs:  data/bulk/GSE313879/gene_count_matrix.tsv.gz
#          data/bulk/GSE313879/sample_metadata.tsv
#          Required metadata columns: sample_id, condition (NC or OE)
# Outputs: output/NB4_CD69OE_DESeq2_all_genes.csv
#          output/NB4_CD69OE_DESeq2_significant.csv
#          figures/NB4_CD69OE_PCA.pdf
# Run order: independent of the single-cell/spatial scripts
# NOTE: Exact original alignment and quantification parameters must be taken
#       from the GEO record/methods; this script begins from the deposited
#       processed count matrix and makes the DE model explicit.
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(DESeq2))

counts_file <- file.path(bulk_dir, "GSE313879", "gene_count_matrix.tsv.gz")
metadata_file <- file.path(bulk_dir, "GSE313879", "sample_metadata.tsv")
require_files(c(counts_file, metadata_file))

counts_df <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts_df) < 7L) stop("Expected a gene column plus six GSE313879 samples.")
gene_col <- colnames(counts_df)[1]
if (anyDuplicated(counts_df[[gene_col]])) {
  stop("Duplicated gene identifiers in the processed GSE313879 count matrix.")
}
rownames(counts_df) <- counts_df[[gene_col]]
counts_df[[gene_col]] <- NULL
counts <- as.matrix(counts_df)
if (anyNA(counts) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
  stop("GSE313879 count matrix must contain non-negative integer counts without missing values.")
}
storage.mode(counts) <- "integer"

metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
if (!all(c("sample_id", "condition") %in% colnames(metadata))) {
  stop("sample_metadata.tsv must contain sample_id and condition columns.")
}
if (!setequal(metadata$sample_id, colnames(counts))) {
  stop("sample_metadata.tsv sample IDs do not exactly match count-matrix columns.")
}
metadata <- metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
group_counts <- table(metadata$condition)
if (!all(c("NC", "OE") %in% names(group_counts)) ||
    any(group_counts[c("NC", "OE")] != 3L) || sum(group_counts) != 6L) {
  stop("Expected exactly three NC and three OE samples.")
}
rownames(metadata) <- metadata$sample_id
metadata$condition <- relevel(factor(metadata$condition), ref = "NC")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)
dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "OE", "NC"), alpha = 0.05)
res_df <- data.frame(gene = rownames(res), as.data.frame(res), row.names = NULL)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.csv(res_df, file.path(output_dir, "NB4_CD69OE_DESeq2_all_genes.csv"), row.names = FALSE)
write.csv(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ],
          file.path(output_dir, "NB4_CD69OE_DESeq2_significant.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_classic()
ggsave(file.path(fig_dir, "NB4_CD69OE_PCA.pdf"), p_pca, width = 6, height = 5)

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "NB4_CD69OE_sessionInfo.txt"))
