library(tximport); library(DESeq2); library(tidyverse)
library(pheatmap); library(ggrepel)

meta <- tribble(
  ~sample,        ~group,
  "SRR33308535", "NTG_HA", "SRR33308536", "NTG_HA", "SRR33308537", "NTG_HA",
  "SRR33308529", "OXR1_act", "SRR33308530", "OXR1_act", "SRR33308531", "OXR1_act",
  "SRR33308532", "NTG_KO", "SRR33308533", "NTG_KO", "SRR33308534", "NTG_KO",
  "SRR33308526", "EMC4_KO", "SRR33308527", "EMC4_KO", "SRR33308528", "EMC4_KO"
) %>% as.data.frame()
meta$group <- factor(meta$group, levels = c("NTG_HA", "OXR1_act", "NTG_KO", "EMC4_KO"))
rownames(meta) <- meta$sample

qdir <- path.expand("~/RNAseq/sra-toolkit/sra/quants")

files <- sapply(meta$sample, function(s) {
  f <- list.files(file.path(qdir, s), pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
  if (length(f) != 1) stop("Expected exactly one quant.sf for ", s, ", found ", length(f))
  f
})
names(files) <- meta$sample
stopifnot(all(file.exists(files)))

# tx2gene from the Ensembl 110 FASTA you indexed
fa <- path.expand("~/sra/cDNA-reference/Homo_sapiens.GRCh38.cdna.all.fa")
stopifnot(file.exists(fa))
hdr <- system(paste("grep '^>'", shQuote(fa)), intern = TRUE)

tx2gene <- tibble(
  tx     = sub("^>(\\S+).*", "\\1", hdr),
  gene   = sub(".*gene:(\\S+).*", "\\1", hdr),
  symbol = ifelse(grepl("gene_symbol:", hdr), sub(".*gene_symbol:(\\S+).*", "\\1", hdr), NA)
)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx", "gene")])

dds <- DESeqDataSetFromTximport(txi, meta, ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# 1. PCA: do samples cluster by group? (expect the two backgrounds to separate)
plotPCA(vsd, intgroup = "group")

# 2. perturbation check
sym <- distinct(tx2gene, gene, symbol) %>% filter(!is.na(symbol))
get_res <- function(num, den) {
  lfcShrink(dds, contrast = c("group", num, den), type = "ashr") %>%
    as.data.frame() %>% rownames_to_column("gene") %>% left_join(sym, by = "gene")
}
res_oxr1 <- get_res("OXR1_act", "NTG_HA")
res_emc4 <- get_res("EMC4_KO", "NTG_KO")
filter(res_oxr1, symbol == "OXR1")   # expect positive log2FC
filter(res_emc4, symbol == "EMC4")   # expect negative log2FC

# 3. volcano
volcano <- function(res, title) {
  ggplot(res, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(colour = padj < 0.05 & abs(log2FoldChange) > 1), alpha = .5) +
    geom_text_repel(data = slice_min(res, padj, n = 10, with_ties = FALSE), aes(label = symbol)) +
    labs(title = title) + theme_minimal() + theme(legend.position = "none")
}
volcano(res_oxr1, "OXR1 activation vs NTG (HA)")
volcano(res_emc4, "EMC4 KO vs NTG (KO)")
