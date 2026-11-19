# Tissue-identity check for LI2024_VIBRIO_GILL: is every library gill?
#
# The 0 h baseline was 3-10x more variable than the infected groups in the
# replicate pre-screen and at full depth. This script records the evidence
# used to exclude 0h_replicate3 (SRR19046239):
#   pca_and_correlation.tsv  PC1-PC3 of VST expression (2,000 most variable
#                            genes) and each library's mean Spearman r with
#                            the libraries of each time point.
#   marker_genes_tpm.tsv     TPM of mantle/shell-formation markers and of
#                            genes the gill libraries share.
# None of these measure the infection contrast itself.
#
# Run from the study worktree root, after the RNA-seq workflow, with the aree
# conda env (DESeq2, tximport):
#   Rscript data/studies/LI2024_VIBRIO_GILL/tissue_identity_check/tissue_identity_check.R
suppressMessages({library(tximport); library(DESeq2)})
out <- "data/studies/LI2024_VIBRIO_GILL/tissue_identity_check"
R <- "results/rnaseq/LI2024_VIBRIO_GILL"
dirs <- file.path(R, c("valginolyticus_48h_vs_0h", "valginolyticus_6h_vs_0h"), "rnaseq/salmon")
ss <- read.delim("data/studies/LI2024_VIBRIO_GILL/samplesheet.tsv")
files <- sapply(ss$run_accession, function(r) {
  f <- file.path(dirs, r, "quant.sf"); f[file.exists(f)][1] })
names(files) <- ss$sample_id
t2g <- read.delim("data/reference/GCF_963853765.1/tx2gene.tsv", header = FALSE)
t2g[, 1] <- sub("[.].*", "", t2g[, 1])
txi <- suppressMessages(tximport(files, type = "salmon", tx2gene = t2g, ignoreTxVersion = TRUE))
ss$timepoint <- factor(ss$condition, levels = c("0h", "6h", "48h")); rownames(ss) <- ss$sample_id
dds <- suppressMessages(DESeqDataSetFromTximport(txi, ss, ~timepoint))
dds <- dds[rowSums(counts(dds)) >= 10, ]
m <- assay(vst(dds, blind = TRUE))
top <- head(order(apply(m, 1, var), decreasing = TRUE), 2000)
p <- prcomp(t(m[top, ])); ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
cc <- cor(m[top, ], method = "spearman")
g <- function(l) ss$sample_id[ss$condition == l]
mean_r <- function(s, set) round(mean(cc[s, setdiff(set, s)]), 3)
tab <- data.frame(sample_id = ss$sample_id, run_accession = ss$run_accession, timepoint = ss$condition,
                  PC1 = round(p$x[ss$sample_id, 1], 1), PC2 = round(p$x[ss$sample_id, 2], 1),
                  PC3 = round(p$x[ss$sample_id, 3], 1),
                  mean_r_0h = sapply(ss$sample_id, mean_r, g("0h")),
                  mean_r_6h = sapply(ss$sample_id, mean_r, g("6h")),
                  mean_r_48h = sapply(ss$sample_id, mean_r, g("48h")))
writeLines(sprintf("# PC1 %s%%, PC2 %s%%, PC3 %s%% of variance (VST, 2,000 most variable genes)",
                   ve[1], ve[2], ve[3]), file.path(out, "pca_and_correlation.tsv"))
suppressWarnings(write.table(tab, file.path(out, "pca_and_correlation.tsv"), sep = "\t",
                             quote = FALSE, row.names = FALSE, append = TRUE))
# Marker genes (NCBI GeneIDs, RS_2024_06 annotation).
# The four mantle/shell genes are among the 15 genes most raised in
# 0h_replicate3 relative to the other eight libraries; ectin and prestin are
# among the 10 most depleted (log2 TPM ratio, genes with max TPM > 50).
markers <- c(
  "105319088" = "putative tyrosinase-like protein tyr-3 (shell formation, mantle edge)",
  "105346498" = "putative per-hexamer repeat protein 5 (shell matrix)",
  "105338611" = "calcium-binding protein LPS1-beta (shell matrix)",
  "105331575" = "tryptase-2",
  "117690964" = "ectin (expressed in the eight other libraries)",
  "105336253" = "prestin (expressed in the eight other libraries)")
tpm <- txi$abundance[intersect(names(markers), rownames(txi$abundance)), , drop = FALSE]
mk <- data.frame(gene_id = rownames(tpm), description = markers[rownames(tpm)], round(tpm, 1), check.names = FALSE)
write.table(mk, file.path(out, "marker_genes_tpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("PC1 %s%%; wrote %d samples and %d marker genes to %s\n", ve[1], nrow(tab), nrow(mk), out))
