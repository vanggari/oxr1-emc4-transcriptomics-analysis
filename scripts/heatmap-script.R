# 1. Tidy a results table: drop unannotated and very low-count genes,
#    and keep one row per symbol (the most significant)
tidy_res <- function(res, min_base = 10) {
  res %>%
    filter(!is.na(padj), !is.na(symbol), baseMean >= min_base) %>%
    group_by(symbol) %>%
    slice_min(padj, n = 1, with_ties = FALSE) %>%
    ungroup()
}
t_oxr1 <- tidy_res(res_oxr1)
t_emc4 <- tidy_res(res_emc4)

# 2. Volcano that always labels the perturbed gene plus the top hits
volcano <- function(res, target, title, n_label = 8) {
  res <- res %>% mutate(
    sig = padj < 0.05 & abs(log2FoldChange) > 1,
    lab = symbol %in% c(target, slice_min(res, padj, n = n_label, with_ties = FALSE)$symbol)
  )
  ggplot(res, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(colour = sig), alpha = 0.5) +
    geom_point(data = filter(res, symbol == target), colour = "black", size = 3) +
    geom_text_repel(data = filter(res, lab), aes(label = symbol),
                    max.overlaps = Inf, box.padding = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c("grey60", "firebrick")) +
    labs(title = title, x = "log2 fold change", y = "-log10(adjusted p)") +
    theme_minimal(base_size = 14) + theme(legend.position = "none")
}
p1 <- volcano(t_oxr1, "OXR1", "OXR1 activation vs NTG (HA)")
p2 <- volcano(t_emc4, "EMC4", "EMC4 KO vs NTG (KO)")
ggsave("volcano_OXR1.png", p1, width = 7, height = 5.5, dpi = 300)
ggsave("volcano_EMC4.png", p2, width = 7, height = 5.5, dpi = 300)

# 3. Heatmaps: top 30 per contrast, using only that contrast's samples
heat <- function(res, groups, title, file) {
  smp <- meta$sample[meta$group %in% groups]
  top <- slice_min(res, padj, n = 30, with_ties = FALSE)
  mat <- assay(vsd)[top$gene, smp]
  rownames(mat) <- top$symbol
  pheatmap(mat, scale = "row", main = title,
           annotation_col = meta[smp, "group", drop = FALSE],
           filename = file, width = 7, height = 8)
}
heat(t_oxr1, c("NTG_HA", "OXR1_act"), "OXR1 activation vs NTG (HA)", "heatmap_OXR1.png")
heat(t_emc4, c("NTG_KO", "EMC4_KO"), "EMC4 KO vs NTG (KO)", "heatmap_EMC4.png")