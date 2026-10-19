# Differential expression with DESeq2 for infection, genotype and their interaction.
# Input:  results/intermediate/counts_raw.rds, results/sample_table.csv,
#         results/intermediate/rlog_blind.rds (heatmap only)
# Output: results/deg/*.csv, results/intermediate/dds.rds, figures/deg/*.png

######### Libraries ##########
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ashr)  # used by lfcShrink(type = "ashr")
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})
source("functions.R")

######### Configuration ##########
ALPHA <- 0.05          # BH-adjusted p-value cutoff
N_VOLCANO_LABELS <- 10
N_HEATMAP_GENES <- 50
DEG_DIR <- file.path("results", "deg")
FIG_DIR <- file.path("figures", "deg")
dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# References are wt and mock, so the main-effect coefficients describe wt cells and uninfected
# cells, and the interaction describes how 6KR changes the infection response.
CONTRASTS <- list(
  infection_wt = list(name = "infection_infected_vs_mock"),
  infection_6KR = list(contrast = list(c("infection_infected_vs_mock",
                                         "genotype6KR.infectioninfected"))),
  genotype_mock = list(name = "genotype_6KR_vs_wt"),
  genotype_infected = list(contrast = list(c("genotype_6KR_vs_wt",
                                             "genotype6KR.infectioninfected"))),
  interaction = list(name = "genotype6KR.infectioninfected")
)
CONTRAST_TITLES <- c(
  infection_wt = "Infected vs mock in wt cells",
  infection_6KR = "Infected vs mock in 6KR cells",
  genotype_mock = "6KR vs wt in mock cells",
  genotype_infected = "6KR vs wt in infected cells",
  interaction = "Interaction: change of the infection effect in 6KR"
)

######### Data loading ##########
counts_raw <- readRDS(file.path("results", "intermediate", "counts_raw.rds"))
rlog_blind <- readRDS(file.path("results", "intermediate", "rlog_blind.rds"))
sample_table <- fread(file.path("results", "sample_table.csv"))
sample_table[, `:=`(
  genotype = factor(genotype, levels = c("wt", "6KR")),
  infection = factor(infection, levels = c("mock", "infected")),
  group = factor(group, levels = names(GROUP_COLOURS))
)]
stopifnot(identical(colnames(counts_raw), sample_table$sample))

######### Preprocessing ##########
counts_filtered <- counts_raw[filter_expressed_genes(counts_raw), ]
stopifnot(identical(rownames(counts_filtered), rownames(rlog_blind)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = data.frame(sample_table, row.names = sample_table$sample, check.names = FALSE),
  design = ~ genotype + infection + genotype:infection
)
dds <- DESeq(dds, quiet = TRUE)
stopifnot(setequal(
  resultsNames(dds),
  c("Intercept", "genotype_6KR_vs_wt", "infection_infected_vs_mock",
    "genotype6KR.infectioninfected")
))

# Gene symbols from org.Hs.eg.db; the Ensembl ID stays the key. Several symbols per ID: keep the
# first (the count is printed below). IDs without a symbol stay in the tables with an empty symbol.
symbol_list <- suppressMessages(mapIds(org.Hs.eg.db, keys = rownames(dds), keytype = "ENSEMBL",
                                       column = "SYMBOL", multiVals = "list"))
symbol_map <- data.table(
  gene_id = names(symbol_list),
  symbol = vapply(symbol_list, function(x) if (all(is.na(x))) "" else x[1], character(1)),
  n_symbols = vapply(symbol_list, function(x) sum(!is.na(x)), integer(1))
)

######### Analysis ##########
run_contrast <- function(spec) {
  res <- do.call(results, c(list(object = dds, alpha = ALPHA), spec))
  shrunk <- suppressMessages(lfcShrink(dds, res = res, type = "ashr", quiet = TRUE))
  tab <- data.table(
    gene_id = rownames(res),
    baseMean = res$baseMean,
    log2FC = res$log2FoldChange,
    log2FC_shrunk = shrunk$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj
  )
  tab <- symbol_map[, .(gene_id, symbol)][tab, on = "gene_id"]
  tab[, direction := fifelse(!is.na(padj) & padj < ALPHA,
                             fifelse(log2FC > 0, "up", "down"), "ns")]
  setorder(tab, padj, na.last = TRUE)
  tab[]
}
deg_tables <- lapply(CONTRASTS, run_contrast)

deg_counts <- rbindlist(lapply(names(deg_tables), function(contrast) {
  tab <- deg_tables[[contrast]]
  data.table(
    contrast = contrast,
    tested = sum(!is.na(tab$padj)),
    up = sum(tab$direction == "up"),
    down = sum(tab$direction == "down"),
    significant = sum(tab$direction != "ns")
  )
}))

# Consistency check: the infection effect in 6KR must equal the wt infection effect plus the
# interaction (unshrunken log2 fold changes of the same fit).
check_genes <- deg_tables$interaction$gene_id[1:5]
lfc_of <- function(contrast) deg_tables[[contrast]][match(check_genes, gene_id), log2FC]
lfc_sum_error <- max(abs(lfc_of("infection_6KR") - (lfc_of("infection_wt") +
                                                     lfc_of("interaction"))))
stopifnot(lfc_sum_error < 1e-6)

# Direction check: classic interferon-stimulated genes must be induced by infection in wt cells.
isg_check <- deg_tables$infection_wt[symbol %in% c("IFIT1", "ISG15", "MX1", "OASL", "IFNB1"),
                                     .(symbol, log2FC_shrunk, padj)]

######### Visualization ##########
plot_volcano <- function(contrast) {
  tab <- deg_tables[[contrast]][!is.na(padj)]
  tab[, neg_log10_padj := -log10(pmax(padj, .Machine$double.xmin))]
  labelled <- tab[direction != "ns"][seq_len(min(.N, N_VOLCANO_LABELS))]
  labelled[, label := fifelse(symbol == "", gene_id, symbol)]
  ggplot(tab, aes(log2FC_shrunk, neg_log10_padj, colour = direction)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
    geom_text_repel(data = labelled, aes(label = label), size = 2.5, family = "Arial",
                    colour = "black", max.overlaps = Inf, min.segment.length = 0,
                    segment.size = 0.2, seed = 42) +
    scale_colour_manual(
      values = c(up = "#B2182B", down = "#2166AC", ns = "grey70"),
      labels = c(up = sprintf("up (%d)", sum(tab$direction == "up")),
                 down = sprintf("down (%d)", sum(tab$direction == "down")),
                 ns = "not significant"),
      breaks = c("up", "down", "ns"), name = NULL
    ) +
    labs(x = "log2 fold change (ashr-shrunken)", y = "-log10 adjusted p-value",
         title = CONTRAST_TITLES[[contrast]]) +
    theme_project() +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
}
volcanoes <- lapply(setNames(names(CONTRASTS), names(CONTRASTS)), plot_volcano)

deg_counts_long <- melt(deg_counts[, .(contrast, up, down)], id.vars = "contrast",
                        variable.name = "direction", value.name = "n_genes")
deg_counts_long[, contrast := factor(contrast, levels = rev(names(CONTRASTS)),
                                     labels = rev(CONTRAST_TITLES[names(CONTRASTS)]))]
p_counts <- ggplot(deg_counts_long, aes(n_genes, contrast, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n_genes), position = position_dodge(width = 0.8), hjust = -0.15,
            size = 2.6, family = "Arial") +
  scale_fill_manual(values = c(up = "#B2182B", down = "#2166AC"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = sprintf("Differentially expressed genes (adjusted p < %.2f)", ALPHA), y = NULL,
       title = "Number of DEGs per contrast") +
  theme_project()

######### Output and summary ##########
for (contrast in names(deg_tables)) {
  fwrite(deg_tables[[contrast]], file.path(DEG_DIR, paste0(contrast, "_all.csv")))
  fwrite(deg_tables[[contrast]][direction != "ns"],
         file.path(DEG_DIR, paste0(contrast, "_significant.csv")))
  save_plot(volcanoes[[contrast]], file.path(FIG_DIR, paste0("volcano_", contrast, ".png")),
            13, 10)
}
fwrite(deg_counts, file.path(DEG_DIR, "deg_counts.csv"))
save_plot(p_counts, file.path(FIG_DIR, "deg_counts.png"), 15, 7)
saveRDS(dds, file.path("results", "intermediate", "dds.rds"))

# Heatmap of the strongest interaction genes on the blind rlog values, scaled per gene.
top_interaction <- deg_tables$interaction[direction != "ns"][seq_len(min(.N, N_HEATMAP_GENES))]
heat_values <- rlog_blind[top_interaction$gene_id, ]
rownames(heat_values) <- fifelse(top_interaction$symbol == "", top_interaction$gene_id,
                                 top_interaction$symbol)
annotation <- data.frame(Group = sample_table$group, row.names = sample_table$sample)
levels(annotation$Group) <- GROUP_LABELS[levels(annotation$Group)]
# Open the ragg device before building, so pheatmap lays out on it (see 02_qc_exploration.R).
ragg::agg_png(file.path(FIG_DIR, "heatmap_interaction_top_genes.png"),
              width = 15, height = 20, units = "cm", res = 300)
interaction_heatmap <- pheatmap(
  heat_values,
  scale = "row",
  clustering_method = "average",
  annotation_col = annotation,
  annotation_colors = list(Group = setNames(GROUP_COLOURS, GROUP_LABELS)),
  annotation_names_col = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  fontsize = 8,
  fontsize_row = 6,
  main = sprintf("Top %d interaction genes (rlog, row z-score)", nrow(heat_values)),
  silent = TRUE
)
grid::grid.newpage()
grid::grid.draw(interaction_heatmap$gtable)
invisible(dev.off())

# Plot objects for the composite report figures (scripts/07_report_figures.R).
dir.create(file.path("results", "intermediate", "plots"), showWarnings = FALSE)
saveRDS(list(volcanoes = volcanoes, deg_counts = p_counts,
             interaction_heatmap = interaction_heatmap$gtable),
        file.path("results", "intermediate", "plots", "plots_03_deg.rds"))

cat("Genes tested:", nrow(dds), "\n")
cat("Symbol mapping: without symbol", sum(symbol_map$n_symbols == 0),
    "| with several symbols (first kept)", sum(symbol_map$n_symbols > 1), "\n\n")
cat("DEGs per contrast (BH-adjusted p <", ALPHA, "):\n")
print(deg_counts)
cat("\nInfection effect in 6KR = wt effect + interaction, max deviation:",
    signif(lfc_sum_error, 3), "\n")
cat("\nInterferon-stimulated genes in infection_wt (should be up):\n")
print(isg_check)
cat("\nWritten:", length(list.files(DEG_DIR)), "tables in", DEG_DIR, "|",
    length(list.files(FIG_DIR)), "figures in", FIG_DIR, "\n")
