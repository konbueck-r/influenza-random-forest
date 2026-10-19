# Candidate biomarker table and expression plot from the Random Forest results.
# Candidates: genes of the final Model 2 (factor_balanced) that rank in its top 20 by importance
# and were selected in at least 10 of the 12 LOOCV folds.
# Input:  results/random_forest/*.csv, results/deg/*_all.csv, results/enrichment/*.csv,
#         results/intermediate/rlog_blind.rds, results/sample_table.csv
# Output: results/biomarkers/candidate_genes.csv, figures/biomarkers/candidate_expression.png

######### Libraries ##########
suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
  library(ggplot2)
})
source("functions.R")

######### Configuration ##########
CANDIDATE_MODEL <- "factor_balanced"
COMPARISON_MODEL <- "group_lrt"
MAX_IMPORTANCE_RANK <- 20
MIN_FOLDS_SELECTED <- 10
MAX_TERMS_LISTED <- 5
DE_CONTRASTS <- c("infection_wt", "genotype_mock", "interaction")
# Ribosomal protein and OXPHOS genes: their apparent induction after infection is most likely
# relative, because they resist the virus-induced host shutoff (see the report, section 4.3).
SHUTOFF_PATTERN <- "^(RP[LS][0-9]+[AX]?$|RPLP[0-2]$|RPSA$|NDUF[ABCSV][0-9]+|COX[4-8]|ATP5(F1|M|P)|UQCR|SDH[ABCD]$)"
OUT_DIR <- file.path("results", "biomarkers")
FIG_DIR <- file.path("figures", "biomarkers")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

######### Data loading ##########
importance <- fread(file.path("results", "random_forest", "final_model_importance.csv"))
stability <- fread(file.path("results", "random_forest", "gene_selection_stability.csv"))
enrichment <- fread(file.path("results", "enrichment", "enrichment_results.csv"))
rlog_blind <- readRDS(file.path("results", "intermediate", "rlog_blind.rds"))
sample_table <- fread(file.path("results", "sample_table.csv"))
sample_table[, group := factor(group, levels = names(GROUP_COLOURS))]
deg <- lapply(setNames(DE_CONTRASTS, DE_CONTRASTS), function(contrast) {
  fread(file.path("results", "deg", paste0(contrast, "_all.csv")))
})

######### Candidate selection ##########
model_importance <- importance[model == CANDIDATE_MODEL]
model_stability <- stability[model == CANDIDATE_MODEL, .(gene_id, folds_selected)]
candidates <- model_stability[model_importance, on = "gene_id"]
candidates[is.na(folds_selected), folds_selected := 0L]
candidates <- candidates[importance_rank <= MAX_IMPORTANCE_RANK &
                           folds_selected >= MIN_FOLDS_SELECTED]
comparison_top <- importance[model == COMPARISON_MODEL & importance_rank <= MAX_IMPORTANCE_RANK,
                             gene_id]

######### Annotation ##########
candidates[, gene_name := suppressMessages(mapIds(org.Hs.eg.db, keys = gene_id,
                                                  keytype = "ENSEMBL", column = "GENENAME",
                                                  multiVals = "first"))]
candidates[, in_model1_top20 := gene_id %in% comparison_top]
candidates[, host_shutoff_caveat := !is.na(symbol) & grepl(SHUTOFF_PATTERN, symbol)]

for (contrast in DE_CONTRASTS) {
  de <- deg[[contrast]][, .(gene_id, log2FC_shrunk, padj)]
  setnames(de, c("log2FC_shrunk", "padj"),
           paste0(contrast, c("_log2FC_shrunk", "_padj")))
  candidates <- de[candidates, on = "gene_id"]
}

# Enriched terms (all ORA lists) that contain the gene, most significant first.
term_hits <- enrichment[, .(gene_id = unlist(strsplit(intersection, ",", fixed = TRUE))),
                        by = .(query, source, term_name, p_adj)]
term_hits <- term_hits[gene_id %in% candidates$gene_id]
setorder(term_hits, gene_id, p_adj)
term_summary <- term_hits[, .(
  n_enriched_terms = .N,
  top_enriched_terms = paste(head(sprintf("%s [%s, %s]", term_name, source, query),
                                  MAX_TERMS_LISTED), collapse = "; ")
), by = gene_id]
candidates <- term_summary[candidates, on = "gene_id"]
candidates[is.na(n_enriched_terms), n_enriched_terms := 0L]

setorder(candidates, importance_rank)
candidate_table <- candidates[, c(
  "importance_rank", "gene_id", "symbol", "gene_name", "dominant_effect", "folds_selected",
  "mean_decrease_accuracy", "best_separated_class", "in_model1_top20", "host_shutoff_caveat",
  paste0(rep(DE_CONTRASTS, each = 2), c("_log2FC_shrunk", "_padj")),
  "n_enriched_terms", "top_enriched_terms"
), with = FALSE]

######### Visualization ##########
expression_long <- data.table(
  gene_id = rep(candidate_table$gene_id, times = ncol(rlog_blind)),
  sample = rep(colnames(rlog_blind), each = nrow(candidate_table)),
  rlog = as.vector(rlog_blind[candidate_table$gene_id, ])
)[sample_table[, .(sample, group)], on = "sample"]
expression_long <- candidate_table[, .(gene_id, symbol, dominant_effect)][
  expression_long, on = "gene_id"]
expression_long[, panel := factor(sprintf("%s\n(%s)", symbol, dominant_effect),
                                  levels = sprintf("%s\n(%s)", candidate_table$symbol,
                                                   candidate_table$dominant_effect))]

p_expression <- ggplot(expression_long, aes(group, rlog, colour = group)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.6, linewidth = 0.3,
               colour = "grey40") +
  geom_point(size = 1.6, position = position_jitter(width = 0.12, seed = 1)) +
  facet_wrap(~ panel, scales = "free_y", ncol = 5) +
  scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = "rlog expression",
       title = "Candidate biomarker genes (Model 2, stable and important)") +
  theme_project() +
  theme(axis.ticks.x = element_blank(), legend.position = "bottom",
        strip.text = element_text(size = 7))

######### Output and summary ##########
fwrite(candidate_table, file.path(OUT_DIR, "candidate_genes.csv"))
save_plot(p_expression, file.path(FIG_DIR, "candidate_expression.png"), 18, 12)
# Plot object for the composite report figures (scripts/07_report_figures.R).
dir.create(file.path("results", "intermediate", "plots"), showWarnings = FALSE)
saveRDS(list(candidate_expression = p_expression),
        file.path("results", "intermediate", "plots", "plots_06_biomarkers.rds"))

cat("Candidates (Model 2, importance rank <=", MAX_IMPORTANCE_RANK, "and selected in >=",
    MIN_FOLDS_SELECTED, "of 12 folds):", nrow(candidate_table), "\n")
print(candidate_table[, .(importance_rank, symbol, dominant_effect, folds_selected,
                          in_model1_top20, host_shutoff_caveat,
                          infection_wt = round(infection_wt_log2FC_shrunk, 2),
                          genotype_mock = round(genotype_mock_log2FC_shrunk, 2),
                          interaction = round(interaction_log2FC_shrunk, 2),
                          n_enriched_terms)])
cat("\nWithout symbol:", candidate_table[is.na(symbol) | symbol == "", .N],
    "| without gene name:", candidate_table[is.na(gene_name), .N], "\n")
