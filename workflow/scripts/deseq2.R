# Gene-level differential expression with DESeq2.
# Contrast is read from config: deseq2.contrast = [condition_of_interest, reference_condition]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

counts <- read.delim(snakemake@input[["counts"]], check.names = FALSE)
gene_info <- counts[, c("gene_id", "gene_name")]
mat <- as.matrix(counts[, setdiff(colnames(counts), c("gene_id", "gene_name"))])
rownames(mat) <- counts$gene_id

samples <- read.delim(snakemake@input[["samples"]], comment.char = "#", colClasses = "character")
rownames(samples) <- samples$sample
samples <- samples[colnames(mat), , drop = FALSE]

contrast <- unlist(snakemake@params[["contrast"]])
covariates <- unlist(snakemake@params[["covariates"]])
alpha <- as.numeric(snakemake@params[["alpha"]])
min_count <- as.numeric(snakemake@params[["min_count"]])

stopifnot(length(contrast) == 2, all(contrast %in% samples$condition))
samples$condition <- relevel(factor(samples$condition), ref = contrast[2])
for (cv in covariates) samples[[cv]] <- factor(samples[[cv]])

# Technical covariates (e.g. sequencing kit) go first, the condition of interest last:
# design ~ kit + condition tests condition while accounting for the kit effect.
design_formula <- as.formula(paste("~", paste(c(covariates, "condition"), collapse = " + ")))
message("Design: ", deparse(design_formula))

dds <- DESeqDataSetFromMatrix(countData = mat, colData = samples, design = design_formula)
# Keep genes with at least min_count reads in at least as many samples as the smallest group
smallest_group <- min(table(samples$condition))
dds <- dds[rowSums(counts(dds) >= min_count) >= smallest_group, ]
message("Genes kept after filtering: ", nrow(dds))

# The parametric dispersion trend needs many genes; tiny test datasets fall back to a constant trend
fit_type <- if (nrow(dds) < 200) "mean" else "parametric"
message("Dispersion fit type: ", fit_type)
dds <- DESeq(dds, fitType = fit_type)
res <- results(dds, contrast = c("condition", contrast[1], contrast[2]), alpha = alpha)
summary(res)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- merge(gene_info, res_df, by = "gene_id")
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.table(res_df, snakemake@output[["results"]], sep = "\t", quote = FALSE, row.names = FALSE)

norm <- as.data.frame(counts(dds, normalized = TRUE))
norm <- cbind(gene_id = rownames(norm), norm)
write.table(norm, snakemake@output[["normalized"]], sep = "\t", quote = FALSE, row.names = FALSE)

# PCA on variance-stabilised counts (vst needs enough genes; fall back to rlog for small tests).
# Page 1: data as measured. Page 2 (only with covariates): after removing the covariate effect
# with limma::removeBatchEffect, for visualisation only (DESeq2 models the covariate directly).
vsd <- tryCatch(vst(dds, blind = FALSE), error = function(e) rlog(dds, blind = FALSE))
shape_by <- if (length(covariates) > 0) covariates[1] else NULL

pca_plot <- function(vsd_obj, title) {
  groups <- c("condition", shape_by)
  pca <- plotPCA(vsd_obj, intgroup = groups, returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  aes_map <- if (is.null(shape_by)) aes(PC1, PC2, color = condition, label = name)
             else aes(PC1, PC2, color = condition, shape = .data[[shape_by]], label = name)
  ggplot(pca, aes_map) +
    geom_point(size = 3) +
    geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"), title = title) +
    theme_bw()
}

pdf(snakemake@output[["pca"]], width = 6.5, height = 5)
print(pca_plot(vsd, "Samples by gene expression (VST)"))
if (length(covariates) > 0) {
  vsd_corrected <- vsd
  covariate_df <- as.data.frame(colData(vsd)[, covariates, drop = FALSE])
  batch_design <- model.matrix(~condition, data = as.data.frame(colData(vsd)))
  assay(vsd_corrected) <- limma::removeBatchEffect(
    assay(vsd),
    batch = covariate_df[[1]],
    batch2 = if (ncol(covariate_df) > 1) covariate_df[[2]] else NULL,
    design = batch_design
  )
  print(pca_plot(vsd_corrected, paste("After removing", paste(covariates, collapse = " + "), "effect")))
}
dev.off()

pdf(snakemake@output[["ma"]], width = 6, height = 5)
plotMA(res, alpha = alpha, main = paste(contrast[1], "vs", contrast[2]))
dev.off()
