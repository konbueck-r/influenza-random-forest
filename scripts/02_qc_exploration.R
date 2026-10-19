# Filter lowly expressed genes, rlog-transform the counts and check sample quality and clustering.
# Input:  results/intermediate/counts_raw.rds, results/sample_table.csv
# Output: results/intermediate/rlog_blind.rds, results/qc_sample_summary.csv,
#         results/qc_pca_variance.csv, figures/qc/*.png

######### Libraries ##########
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(ggdendro)
  library(patchwork)
  library(pheatmap)
})
source("functions.R")

######### Configuration ##########
N_TOP_VAR_GENES <- 500  # genes used for the PCA, as in DESeq2's plotPCA()
FIG_DIR <- file.path("figures", "qc")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

######### Data loading ##########
counts_raw <- readRDS(file.path("results", "intermediate", "counts_raw.rds"))
sample_table <- fread(file.path("results", "sample_table.csv"))
sample_table[, group := factor(group, levels = names(GROUP_COLOURS))]
stopifnot(identical(colnames(counts_raw), sample_table$sample))

######### Preprocessing ##########
keep <- filter_expressed_genes(counts_raw)
counts_filtered <- counts_raw[keep, ]

# Design ~ 1 because this object is only used for size factors and the blind transformation.
dds_qc <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = data.frame(sample_table, row.names = sample_table$sample, check.names = FALSE),
  design = ~ 1
)
dds_qc <- estimateSizeFactors(dds_qc)
rlog_blind <- assay(rlog(dds_qc, blind = TRUE))

######### Analysis ##########
# PCA on the most variable genes, centred but not scaled (same as plotPCA()).
top_var <- order(apply(rlog_blind, 1, var), decreasing = TRUE)[seq_len(N_TOP_VAR_GENES)]
pca <- prcomp(t(rlog_blind[top_var, ]), center = TRUE, scale. = FALSE)
pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_scores <- data.table(sample = rownames(pca$x), pca$x[, 1:3])[sample_table, on = "sample"]

# Euclidean distances between samples on all filtered genes, and average-linkage clustering.
sample_dist <- dist(t(rlog_blind))
sample_clust <- hclust(sample_dist, method = "average")

# Outlier check: a sample whose nearest neighbour is not a replicate of the same group would be
# suspicious. Also compare its mean distance to its own replicates with the closest other group.
dist_mat <- as.matrix(sample_dist)
outlier_check <- rbindlist(lapply(sample_table$sample, function(s) {
  own_group <- sample_table[sample == s, group]
  others <- setdiff(colnames(dist_mat), s)
  nearest <- others[which.min(dist_mat[s, others])]
  replicates <- setdiff(sample_table[group == own_group, sample], s)
  other_groups <- setdiff(levels(sample_table$group), as.character(own_group))
  closest_other <- min(vapply(other_groups, function(g) {
    mean(dist_mat[s, sample_table[group == g, sample]])
  }, numeric(1)))
  data.table(
    sample = s,
    nearest_neighbour = nearest,
    nearest_is_replicate = nearest %in% replicates,
    mean_dist_replicates = mean(dist_mat[s, replicates]),
    mean_dist_closest_other_group = closest_other
  )
}))

# Does the grouping hold on the most variable genes? Shallow libraries add noise mainly to weakly
# expressed genes, which dominate distances computed on all genes.
cluster_top_var <- cutree(hclust(dist(t(rlog_blind[top_var, ])), method = "average"), k = 4)
cluster_all <- cutree(sample_clust, k = 4)

qc_summary <- data.table(
  sample = sample_table$sample,
  group = sample_table$group,
  library_size = colSums(counts_raw),
  genes_detected = colSums(counts_raw > 0),
  share_filtered_genes_below_5_counts = colMeans(counts_filtered < 5),
  size_factor = sizeFactors(dds_qc),
  cluster_top_var_genes = cluster_top_var[sample_table$sample],
  cluster_all_genes = cluster_all[sample_table$sample]
)[pca_scores[, .(sample, PC1, PC2, PC3)], on = "sample"][outlier_check, on = "sample"]

######### Visualization ##########
group_fill <- scale_fill_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL)
group_colour <- scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL)

p_libsize <- ggplot(qc_summary, aes(factor(sample, levels = unique(sample)), library_size / 1e6,
                                    fill = group)) +
  geom_col(width = 0.75) +
  group_fill +
  labs(x = NULL, y = "Assigned reads (million)", title = "Library sizes") +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rlog_long <- data.table(
  sample = rep(colnames(rlog_blind), each = nrow(rlog_blind)),
  value = as.vector(rlog_blind)
)[sample_table[, .(sample, group)], on = "sample"]
p_box <- ggplot(rlog_long, aes(factor(sample, levels = unique(sample)), value, fill = group)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  group_fill +
  labs(x = NULL, y = "rlog expression", title = "Expression distribution after rlog") +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_pca <- function(x_pc, y_pc) {
  ggplot(pca_scores, aes(.data[[x_pc]], .data[[y_pc]], colour = group)) +
    geom_point(size = 2.5) +
    group_colour +
    labs(
      x = sprintf("%s (%.1f %% variance)", x_pc, pct_var[as.integer(sub("PC", "", x_pc))]),
      y = sprintf("%s (%.1f %% variance)", y_pc, pct_var[as.integer(sub("PC", "", y_pc))])
    ) +
    theme_project()
}
p_pca <- (plot_pca("PC1", "PC2") + plot_pca("PC2", "PC3")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("PCA of the %d most variable genes (rlog)", N_TOP_VAR_GENES),
    theme = theme_project()
  )

scree <- data.table(pc = seq_along(pct_var), pct = pct_var)[pc <= 11]
p_scree <- ggplot(scree, aes(factor(pc), pct)) +
  geom_col(fill = "grey45", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -0.4, size = 2.6, family = "Arial") +
  labs(x = "Principal component", y = "Variance explained (%)", title = "Scree plot") +
  theme_project()

# dendro_data() queries par(), which would open a default device and write Rplots.pdf; a
# file-less pdf device absorbs that.
pdf(NULL)
dendro <- dendro_data(as.dendrogram(sample_clust))
invisible(dev.off())
dendro_labels <- as.data.table(label(dendro))[sample_table[, .(label = sample, group)],
                                               on = "label"]
p_dendro <- ggplot() +
  geom_segment(data = segment(dendro), aes(x, y, xend = xend, yend = yend), linewidth = 0.4) +
  geom_text(data = dendro_labels, aes(x, y - 1, label = label, colour = group),
            angle = 90, hjust = 1, size = 2.8, family = "Arial", show.legend = FALSE) +
  group_colour +
  scale_y_continuous(breaks = function(limits) pretty(c(0, limits[2])),
                     expand = expansion(mult = c(0.45, 0.05))) +
  labs(x = NULL, y = "Euclidean distance (rlog)",
       title = "Hierarchical clustering of samples (average linkage)") +
  theme_project() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank())

######### Output and summary ##########
saveRDS(rlog_blind, file.path("results", "intermediate", "rlog_blind.rds"))
fwrite(qc_summary, file.path("results", "qc_sample_summary.csv"))
fwrite(data.table(pc = seq_along(pct_var), percent_variance = pct_var),
       file.path("results", "qc_pca_variance.csv"))

save_plot(p_libsize, file.path(FIG_DIR, "library_sizes.png"), 14, 8)
save_plot(p_box, file.path(FIG_DIR, "rlog_boxplot.png"), 14, 8)
save_plot(p_pca, file.path(FIG_DIR, "pca.png"), 17, 7.5)
save_plot(p_scree, file.path(FIG_DIR, "pca_scree.png"), 10, 7)
save_plot(p_dendro, file.path(FIG_DIR, "sample_dendrogram.png"), 14, 9)

annotation <- data.frame(Group = sample_table$group, row.names = sample_table$sample)
levels(annotation$Group) <- GROUP_LABELS[levels(annotation$Group)]
# Open the ragg device first so pheatmap lays out the heatmap on it; build silently and draw
# explicitly, because pheatmap's own drawing left the PNG empty and wrote to Rplots.pdf.
ragg::agg_png(file.path(FIG_DIR, "sample_distance_heatmap.png"),
              width = 14, height = 11.5, units = "cm", res = 300)
distance_heatmap <- pheatmap(
  dist_mat,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  clustering_method = "average",
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = list(Group = setNames(GROUP_COLOURS, GROUP_LABELS)),
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames = FALSE,  # the matrix is symmetric; the column labels name the samples
  cellwidth = 11,
  cellheight = 11,  # square cells for a square distance matrix
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  color = colorRampPalette(c("#08306B", "#6BAED6", "white"))(100),
  fontsize = 8,  # ragg draws the default "sans" family as Arial on Windows
  main = "Sample-to-sample distances (rlog, all filtered genes)",
  silent = TRUE
)
grid::grid.newpage()
grid::grid.draw(distance_heatmap$gtable)
invisible(dev.off())

# Plot objects for the composite report figures (scripts/07_report_figures.R).
dir.create(file.path("results", "intermediate", "plots"), showWarnings = FALSE)
saveRDS(list(pca = p_pca, scree = p_scree, library_sizes = p_libsize, rlog_boxplot = p_box,
             dendrogram = p_dendro, distance_heatmap = distance_heatmap$gtable),
        file.path("results", "intermediate", "plots", "plots_02_qc.rds"))

cat("Genes before filtering:", nrow(counts_raw), "| after (>= 10 counts in >= 3 samples):",
    nrow(counts_filtered), "| removed:", sum(!keep), "\n")
cat("Variance explained by PC1-PC3 (%):", pct_var[1:3], "\n")
cat("\nOutlier check (nearest neighbour and distances):\n")
print(qc_summary[, .(sample, library_size, size_factor = round(size_factor, 2),
                     low_count_share = round(share_filtered_genes_below_5_counts, 3),
                     nearest_neighbour, nearest_is_replicate,
                     dist_repl = round(mean_dist_replicates, 1),
                     dist_other = round(mean_dist_closest_other_group, 1),
                     cluster_top_var_genes, cluster_all_genes)])
matches_groups <- function(clusters) {
  all(tapply(clusters, sample_table$group, function(x) length(unique(x))) == 1) &&
    length(unique(clusters)) == 4
}
cat("4 clusters reproduce the 4 groups - top", N_TOP_VAR_GENES, "genes:",
    matches_groups(cluster_top_var), "| all filtered genes:", matches_groups(cluster_all), "\n")
cat("\nWritten: results/intermediate/rlog_blind.rds, results/qc_sample_summary.csv,",
    length(list.files(FIG_DIR)), "figures in", FIG_DIR, "\n")
