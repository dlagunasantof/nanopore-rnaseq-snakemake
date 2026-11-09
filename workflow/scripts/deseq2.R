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
alpha <- as.numeric(snakemake@params[["alpha"]])
min_count <- as.numeric(snakemake@params[["min_count"]])

stopifnot(length(contrast) == 2, all(contrast %in% samples$condition))
samples$condition <- relevel(factor(samples$condition), ref = contrast[2])

dds <- DESeqDataSetFromMatrix(countData = mat, colData = samples, design = ~condition)
# Keep genes with at least min_count reads in at least as many samples as the smallest group
smallest_group <- min(table(samples$condition))
dds <- dds[rowSums(counts(dds) >= min_count) >= smallest_group, ]
message("Genes kept after filtering: ", nrow(dds))

dds <- DESeq(dds)
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

# PCA on variance-stabilised counts (vst needs enough genes; fall back to rlog for small tests)
vsd <- tryCatch(vst(dds, blind = TRUE), error = function(e) rlog(dds, blind = TRUE))
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"),
       title = "Samples by gene expression (VST)") +
  theme_bw()
ggsave(snakemake@output[["pca"]], p, width = 6, height = 5)

pdf(snakemake@output[["ma"]], width = 6, height = 5)
plotMA(res, alpha = alpha, main = paste(contrast[1], "vs", contrast[2]))
dev.off()
