#!/usr/bin/env Rscript

# ============================================================
# Title: Protein DEA Pipeline (limma)
# Author: Keerthanaa Balasubramanian Shanthi
# Date: 5.8.2025
# Description: Differential expression analysis of proteomics data from DIA proteomics
# ============================================================

# --- 1. Load Packages ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(imputeLCMD)
  library(pheatmap)
  library(RColorBrewer)
})

# --- 2. Configurable Paths ---
input_dir <- "data"
output_dir <- "results"

protein_file <- file.path(input_dir, "Protein_Quant.tsv")
sampleinfo_file <- file.path(input_dir, "SampleInfo.tsv")

dir.create(output_dir, showWarnings = FALSE)

# --- 3. Load Data ---
protein_data <- read_tsv(protein_file)
sample_info <- read_tsv(sampleinfo_file)

sample_cols <- paste0("directLFQ.", sample_info$SampleID)
intensity_matrix <- protein_data %>%
  select(Protein.Ids, all_of(sample_cols)) %>%
  column_to_rownames("Protein.Ids") %>%
  as.matrix()

colnames(intensity_matrix) <- sample_info$SampleID

# --- 4. Log2 + MinDet Imputation ---
log2_matrix <- log2(intensity_matrix + 1)
log2_matrix[is.infinite(log2_matrix)] <- NA
log2_imputed <- impute.MinDet(log2_matrix)

# --- 5. limma Design ---
group <- factor(sample_info$Group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log2_imputed, design)

# --- 6. Contrasts ---
group_levels <- levels(group)
contrast_names <- combn(group_levels, 2, FUN = function(x) paste0(x[2], "_vs_", x[1]))
contrast_formulas <- combn(group_levels, 2, FUN = function(x) paste0(x[2], "-", x[1]))
contrast_matrix <- makeContrasts(contrasts = contrast_formulas, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# --- 7. Save Results & Plots ---
all_up_down_combined <- list()
all_contrast_results <- list()

for (i in seq_along(contrast_formulas)) {
  contrast <- contrast_formulas[i]
  res <- topTable(fit2, coef = contrast, number = Inf, sort.by = "P") %>%
    rownames_to_column("ProteinID") %>%
    mutate(Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"),
           contrast = contrast)

  all_contrast_results[[contrast]] <- res

  write_csv(res, file.path(output_dir, paste0("DE_results_", contrast, ".csv")))

  ggplot(res, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = paste("Volcano Plot:", contrast), x = "log2 FC", y = "-log10 P-value") +
    theme_minimal()
  ggsave(file.path(output_dir, paste0("volcano_", contrast, ".png")), width = 6, height = 5)

  up <- res %>% filter(logFC > 1, adj.P.Val < 0.05) %>% slice_max(logFC, n = 20)
  down <- res %>% filter(logFC < -1, adj.P.Val < 0.05) %>% slice_min(logFC, n = 20)
  all_up_down_combined[[contrast]] <- bind_rows(up, down)

  contrast_data <- all_up_down_combined[[contrast]] %>%
    mutate(Regulation = ifelse(logFC > 0, "Upregulated", "Downregulated"))

  ggplot(contrast_data, aes(x = reorder(ProteinID, logFC), y = logFC, fill = Regulation)) +
    geom_col() +
    scale_fill_manual(values = c("Upregulated" = "darkgreen", "Downregulated" = "darkred")) +
    coord_flip() +
    labs(title = paste("Top 20 Up & Down Proteins:", contrast),
         x = "Protein", y = "log2 Fold Change") +
    theme_minimal()
  ggsave(file.path(output_dir, paste0("barplot_up_down_", contrast, ".png")), width = 7, height = 6)
}

# --- 8. Heatmaps ---
top_all <- bind_rows(all_up_down_combined) %>%
  distinct(ProteinID, .keep_all = TRUE) %>%
  slice_max(abs(logFC), n = 200)

heatmap_data <- log2_imputed[rownames(log2_imputed) %in% top_all$ProteinID, ]
annotation_col <- data.frame(Condition = sample_info$Group)
rownames(annotation_col) <- sample_info$SampleID

pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  filename = file.path(output_dir, "heatmap_top200_combined.png"),
  width = 9,
  height = 10
)

combined_res <- bind_rows(all_contrast_results) %>%
  group_by(ProteinID) %>%
  arrange(adj.P.Val) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Significance = case_when(
    logFC > 1 & adj.P.Val < 0.05 ~ "Upregulated",
    logFC < -1 & adj.P.Val < 0.05 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

ordered_proteins <- combined_res %>%
  arrange(factor(Significance, levels = c("Upregulated", "Downregulated", "Not Significant")),
          desc(abs(logFC))) %>%
  pull(ProteinID)

ordered_proteins <- intersect(ordered_proteins, rownames(log2_imputed))
heatmap_matrix <- log2_imputed[ordered_proteins, , drop = FALSE]

hm_height_all <- pmin(30, 0.015 * nrow(heatmap_matrix) + 4)

pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",
  show_rownames = FALSE,
  cluster_rows = FALSE,
  clustering_distance_cols = "euclidean",
  filename = file.path(output_dir, "heatmap_up_down_nonsig_ordered.png"),
  width = 9,
  height = hm_height_all
)

# --- 9. Save session info ---
writeLines(capture.output(sessionInfo()), file.path(output_dir, "../session_info.txt"))
