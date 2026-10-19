# Random Forest classification of the four groups with nested leave-one-out cross-validation.
# In every fold, genes are selected with DESeq2 on the training samples only, so the left-out
# sample never influences which genes the model sees. Two gene-selection rules are compared:
#   group_lrt       - planned model: top genes of a likelihood-ratio test ~ group vs ~ 1
#   factor_balanced - added after the planned model's errors turned out to be genotype errors:
#                     half the genes from the genotype effect, half from the infection effect
#                     of an additive model ~ genotype + infection (post-hoc, so optimistic)
# Input:  results/intermediate/counts_raw.rds, results/intermediate/rlog_blind.rds,
#         results/sample_table.csv, results/deg/<contrast>_all.csv (symbols, DE status)
# Output: results/random_forest/*.csv, figures/random_forest/*.png

######### Libraries ##########
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(randomForest)
  library(ggplot2)
  library(patchwork)
})
source("functions.R")

######### Configuration ##########
RANDOM_SEED <- 42
N_FEATURES <- 50                          # genes per model, fixed before any result was seen
N_FEATURES_SENSITIVITY <- c(20, 50, 200)  # descriptive check only, no tuning
N_TREES <- 1000
N_IMPORTANCE_PLOT <- 20
MODELS <- c("group_lrt", "factor_balanced")
MODEL_TITLES <- c(group_lrt = "Model 1: LRT across the four groups",
                  factor_balanced = "Model 2: genotype + infection genes")
OUT_DIR <- file.path("results", "random_forest")
FIG_DIR <- file.path("figures", "random_forest")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

######### Data loading ##########
counts_raw <- readRDS(file.path("results", "intermediate", "counts_raw.rds"))
rlog_blind <- readRDS(file.path("results", "intermediate", "rlog_blind.rds"))
sample_table <- fread(file.path("results", "sample_table.csv"))
sample_table[, `:=`(
  genotype = factor(genotype, levels = c("wt", "6KR")),
  infection = factor(infection, levels = c("mock", "infected")),
  group = factor(group, levels = names(GROUP_COLOURS))
)]
read_deg <- function(contrast) fread(file.path("results", "deg", paste0(contrast, "_all.csv")))
symbols <- read_deg("infection_wt")[, .(gene_id, symbol)]

# Dominant effect of each gene in the full-data DE analysis, used only to annotate plots: the
# factor with the larger significant shrunken fold change (infection in either genotype vs
# genotype in either infection state). Almost every gene is an infection DEG somewhere, so a
# plain "significant or not" label would not tell the two apart.
max_significant_lfc <- function(contrasts) {
  tab <- rbindlist(lapply(contrasts, function(x) read_deg(x)[, .(gene_id, log2FC_shrunk,
                                                                 direction)]))
  tab[, .(max_abs_lfc = max(fifelse(direction != "ns", abs(log2FC_shrunk), 0))), by = gene_id]
}
infection_lfc <- max_significant_lfc(c("infection_wt", "infection_6KR"))
genotype_lfc <- max_significant_lfc(c("genotype_mock", "genotype_infected"))
dominant_effect <- merge(infection_lfc[, .(gene_id, infection = max_abs_lfc)],
                         genotype_lfc[, .(gene_id, genotype = max_abs_lfc)], by = "gene_id")
dominant_effect[, dominant_effect := fcase(
  infection == 0 & genotype == 0, "not significant",
  infection >= genotype, "infection",
  default = "genotype"
)]

# The expression filter uses no labels, so applying it once to all samples does not leak class
# information; the same holds for the blind rlog transformation.
counts_filtered <- counts_raw[filter_expressed_genes(counts_raw), ]
stopifnot(identical(rownames(counts_filtered), rownames(rlog_blind)),
          identical(colnames(counts_filtered), sample_table$sample))
samples <- sample_table$sample
labels <- setNames(sample_table$group, samples)

######### Helper functions ##########
fit_deseq <- function(train_samples, design, ...) {
  col_data <- data.frame(sample_table[match(train_samples, sample)], row.names = train_samples,
                         check.names = FALSE)
  dds <- DESeqDataSetFromMatrix(counts_filtered[, train_samples], col_data, design = design)
  DESeq(dds, quiet = TRUE, ...)
}

rank_by_pvalue <- function(res) rownames(res)[order(res$pvalue, -abs(res$stat), na.last = NA)]

# Gene rankings of one fold for both selection rules, from the training samples only.
rank_genes <- function(train_samples) {
  dds_lrt <- fit_deseq(train_samples, ~ group, test = "LRT", reduced = ~ 1)
  dds_add <- fit_deseq(train_samples, ~ genotype + infection)
  list(
    group_lrt = rank_by_pvalue(results(dds_lrt)),
    genotype = rank_by_pvalue(results(dds_add, name = "genotype_6KR_vs_wt")),
    infection = rank_by_pvalue(results(dds_add, name = "infection_infected_vs_mock")),
    fitted_samples = union(colnames(dds_lrt), colnames(dds_add))
  )
}

select_genes <- function(ranking, model, n) {
  if (model == "group_lrt") return(ranking$group_lrt[seq_len(n)])
  # Half from the genotype ranking, then fill up with infection genes not already chosen.
  genotype_part <- ranking$genotype[seq_len(n %/% 2)]
  infection_part <- setdiff(ranking$infection, genotype_part)[seq_len(n - length(genotype_part))]
  c(genotype_part, infection_part)
}

fit_forest <- function(genes, train_samples, seed) {
  set.seed(seed)
  randomForest(x = t(rlog_blind[genes, train_samples]), y = labels[train_samples],
               ntree = N_TREES, importance = TRUE)
}

# Sensitivity, specificity and balanced accuracy per class from a truth/prediction pair.
per_class_metrics <- function(truth, predicted) {
  rbindlist(lapply(levels(truth), function(cls) {
    tp <- sum(predicted == cls & truth == cls)
    fn <- sum(predicted != cls & truth == cls)
    tn <- sum(predicted != cls & truth != cls)
    fp <- sum(predicted == cls & truth != cls)
    data.table(class = cls, sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
               balanced_accuracy = (tp / (tp + fn) + tn / (tn + fp)) / 2)
  }))
}

infection_of <- function(group) sub(".*_", "", as.character(group))
genotype_of <- function(group) sub("_.*", "", as.character(group))

summarise_predictions <- function(pred) {
  n_correct <- sum(pred$predicted == pred$truth)
  n_samples <- nrow(pred)
  ci <- binom.test(n_correct, n_samples)$conf.int
  # One-sided test against guessing among four balanced classes (no-information rate 0.25).
  p_chance <- binom.test(n_correct, n_samples, p = 0.25, alternative = "greater")$p.value
  infected_truth <- infection_of(pred$truth) == "infected"
  infected_pred <- infection_of(pred$predicted) == "infected"
  infected_rows <- pred[infection_of(truth) == "infected"]
  mock_rows <- pred[infection_of(truth) == "mock"]
  data.table(
    correct = n_correct,
    samples = n_samples,
    accuracy = n_correct / n_samples,
    accuracy_ci_lower = ci[1],
    accuracy_ci_upper = ci[2],
    p_vs_chance = p_chance,
    infection_accuracy = mean(infected_truth == infected_pred),
    infection_sensitivity = sum(infected_truth & infected_pred) / sum(infected_truth),
    infection_specificity = sum(!infected_truth & !infected_pred) / sum(!infected_truth),
    genotype_accuracy_mock = mean(genotype_of(mock_rows$truth) ==
                                    genotype_of(mock_rows$predicted)),
    genotype_accuracy_infected = mean(genotype_of(infected_rows$truth) ==
                                        genotype_of(infected_rows$predicted))
  )
}

######### Analysis: nested LOOCV ##########
fold_rankings <- setNames(vector("list", length(samples)), samples)
for (i in seq_along(samples)) {
  ranking <- rank_genes(samples[-i])
  # Leakage guard: the left-out sample must not be part of the fold's DESeq2 fits.
  stopifnot(!samples[i] %in% ranking$fitted_samples, length(ranking$fitted_samples) == 11)
  fold_rankings[[samples[i]]] <- ranking
  cat("Fold", i, "of", length(samples), "- gene rankings done for left-out", samples[i], "\n")
}

run_loocv <- function(model, n_features) {
  rbindlist(lapply(seq_along(samples), function(i) {
    test_sample <- samples[i]
    genes <- select_genes(fold_rankings[[test_sample]], model, n_features)
    forest <- fit_forest(genes, samples[-i], seed = RANDOM_SEED + i)
    votes <- predict(forest, t(rlog_blind[genes, test_sample, drop = FALSE]), type = "prob")
    data.table(
      model = model,
      n_features = n_features,
      sample = test_sample,
      truth = labels[[test_sample]],
      predicted = factor(colnames(votes)[which.max(votes)], levels = levels(labels)),
      votes_true_class = votes[1, as.character(labels[[test_sample]])],
      as.data.table(votes)
    )
  }))
}
runs <- CJ(model = MODELS, n_features = N_FEATURES_SENSITIVITY)
loocv_all <- rbindlist(Map(run_loocv, runs$model, runs$n_features))
loocv <- loocv_all[n_features == N_FEATURES]

######### Analysis: performance ##########
summary_all <- loocv_all[, summarise_predictions(.SD), by = .(model, n_features)]
# BH adjustment over the two main-model tests (one per selection rule).
summary_main <- summary_all[n_features == N_FEATURES]
summary_main[, p_vs_chance_bh := p.adjust(p_vs_chance, method = "BH")]
sensitivity_table <- summary_all[, .(model, n_features, correct, samples, accuracy)]

class_metrics <- loocv[, per_class_metrics(truth, predicted), by = model]
confusion <- loocv[, as.data.table(table(predicted = predicted, truth = truth)), by = model]

######### Analysis: final models and gene stability ##########
full_ranking <- rank_genes(samples)
final_fits <- lapply(setNames(MODELS, MODELS), function(model) {
  genes <- select_genes(full_ranking, model, N_FEATURES)
  list(genes = genes, forest = fit_forest(genes, samples, seed = RANDOM_SEED))
})

importance_of <- function(model) {
  forest <- final_fits[[model]]$forest
  # randomForest's default: the mean decrease in accuracy is divided by its standard error.
  imp <- importance(forest)
  tab <- data.table(model = model, gene_id = rownames(imp),
                    mean_decrease_accuracy = imp[, "MeanDecreaseAccuracy"],
                    mean_decrease_gini = imp[, "MeanDecreaseGini"])
  class_imp <- imp[, levels(labels), drop = FALSE]
  tab[, best_separated_class := colnames(class_imp)[max.col(class_imp, "first")]]
  tab <- dominant_effect[, .(gene_id, dominant_effect)][tab, on = "gene_id"]
  tab <- symbols[tab, on = "gene_id"]
  setorder(tab, -mean_decrease_accuracy)
  tab[, importance_rank := .I]
  tab[]
}
importance_table <- rbindlist(lapply(MODELS, importance_of))

stability <- rbindlist(lapply(MODELS, function(model) {
  selected <- unlist(lapply(fold_rankings, select_genes, model = model, n = N_FEATURES))
  tab <- data.table(gene_id = selected)[, .(folds_selected = .N), by = gene_id]
  tab[, `:=`(model = model, in_final_model = gene_id %in% final_fits[[model]]$genes)]
  tab
}))
stability <- symbols[stability, on = "gene_id"]
setcolorder(stability, "model")
setorder(stability, model, -folds_selected, gene_id)

oob_confusion <- rbindlist(lapply(MODELS, function(model) {
  conf <- final_fits[[model]]$forest$confusion
  data.table(model = model, truth = rownames(conf), conf)
}))

######### Visualization ##########
plot_confusion <- function(model_name) {
  dt <- confusion[model == model_name]
  dt[, `:=`(predicted = factor(predicted, levels = rev(levels(labels))),
            truth = factor(truth, levels = levels(labels)))]
  n_correct <- dt[as.character(predicted) == as.character(truth), sum(N)]
  ggplot(dt, aes(truth, predicted, fill = N)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = N), size = 3.2, family = "Arial") +
    scale_fill_gradient(low = "white", high = "#1F5A96", limits = c(0, 3), guide = "none") +
    scale_x_discrete(labels = GROUP_LABELS) +
    scale_y_discrete(labels = GROUP_LABELS) +
    labs(x = "True group", y = "Predicted group",
         title = sprintf("%s\n%d of %d correct", MODEL_TITLES[[model_name]], n_correct,
                         sum(dt$N))) +
    theme_project() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(size = 9))
}
p_confusion <- plot_confusion("group_lrt") + plot_confusion("factor_balanced")

# Model order is stored in the data, so the saved plot object does not depend on MODELS.
sensitivity_plot_data <- copy(sensitivity_table)[, model := factor(model, levels = MODELS)]
p_sensitivity <- ggplot(sensitivity_plot_data, aes(factor(n_features), correct, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%d/12", correct)), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 2.6, family = "Arial") +
  scale_fill_manual(values = c(group_lrt = "grey55", factor_balanced = "#1F5A96"),
                    labels = MODEL_TITLES, name = NULL) +
  scale_y_continuous(limits = c(0, 12.8), breaks = seq(0, 12, 3)) +
  labs(x = "Genes per fold", y = "Correct LOOCV predictions",
       title = "Number of selected genes (descriptive, not tuned)") +
  theme_project() +
  theme(legend.position = "bottom", legend.direction = "vertical")

plot_importance <- function(model_name) {
  top <- importance_table[model == model_name][seq_len(N_IMPORTANCE_PLOT)]
  top[, label := fifelse(is.na(symbol) | symbol == "", gene_id, symbol)]
  top[, label := factor(label, levels = rev(label))]
  ggplot(top, aes(mean_decrease_accuracy, label, colour = dominant_effect)) +
    geom_segment(aes(x = 0, xend = mean_decrease_accuracy, yend = label), linewidth = 0.4) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(infection = "#B2182B", genotype = "#2166AC",
                                   `not significant` = "grey60"),
                        name = "Dominant DE effect", drop = FALSE) +
    labs(x = "Mean decrease in accuracy (scaled)", y = NULL,
         title = sprintf("%s: top %d genes", MODEL_TITLES[[model_name]], N_IMPORTANCE_PLOT)) +
    theme_project()
}

######### Output and summary ##########
fwrite(loocv_all, file.path(OUT_DIR, "loocv_predictions.csv"))
fwrite(confusion, file.path(OUT_DIR, "loocv_confusion_matrices.csv"))
fwrite(class_metrics, file.path(OUT_DIR, "loocv_class_metrics.csv"))
fwrite(summary_main, file.path(OUT_DIR, "loocv_summary_metrics.csv"))
fwrite(sensitivity_table, file.path(OUT_DIR, "loocv_n_features_sensitivity.csv"))
fwrite(importance_table, file.path(OUT_DIR, "final_model_importance.csv"))
fwrite(stability, file.path(OUT_DIR, "gene_selection_stability.csv"))
fwrite(oob_confusion, file.path(OUT_DIR, "final_model_oob_confusion.csv"))

save_plot(p_confusion, file.path(FIG_DIR, "loocv_confusion_matrices.png"), 17, 8.5)
save_plot(p_sensitivity, file.path(FIG_DIR, "loocv_n_features_sensitivity.png"), 10, 9)
# Plot objects for the composite report figures (scripts/07_report_figures.R).
dir.create(file.path("results", "intermediate", "plots"), showWarnings = FALSE)
saveRDS(list(confusion = list(group_lrt = plot_confusion("group_lrt"),
                              factor_balanced = plot_confusion("factor_balanced")),
             sensitivity = p_sensitivity,
             importance = lapply(setNames(MODELS, MODELS), plot_importance)),
        file.path("results", "intermediate", "plots", "plots_05_rf.rds"))
for (model in MODELS) {
  save_plot(plot_importance(model),
            file.path(FIG_DIR, paste0("final_model_importance_", model, ".png")), 15, 11)
}

cat("\nLOOCV summary with", N_FEATURES, "genes per fold:\n")
print(summary_main[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
cat("\nSensitivity to the number of genes (descriptive, not tuned):\n")
print(dcast(sensitivity_table, n_features ~ model, value.var = "correct"))
cat("\nMisclassified samples (", N_FEATURES, "genes):\n")
print(loocv[predicted != truth, .(model, sample, predicted,
                                  votes_true_class = round(votes_true_class, 2))])
cat("\nFinal-model OOB confusion:\n")
print(oob_confusion)
cat("\nGene selection stability (genes selected in all 12 / >= 10 folds / at least once):\n")
print(stability[, .(all_12 = sum(folds_selected == 12), at_least_10 = sum(folds_selected >= 10),
                    distinct = .N), by = model])
cat("\nTop 10 genes per final model:\n")
print(importance_table[importance_rank <= 10, .(model, symbol, mean_decrease_accuracy = round(
  mean_decrease_accuracy, 2), dominant_effect)])
