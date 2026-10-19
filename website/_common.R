# Shared setup for the website pages: result tables, number formatting and plot helpers.
# Every number on the site is read from ../results at render time, as in the report.

suppressPackageStartupMessages({
  library(data.table)
  library(plotly)
})

######### Paths ##########
results_path <- function(...) file.path("..", "results", ...)

######### Formatting ##########
n <- function(x) formatC(x, format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)

# Same group colours and labels as the report figures (functions.R of the pipeline).
GROUP_COLOURS <- c(wt_mock = "#8DB6DD", wt_infected = "#1F5A96",
                   `6KR_mock` = "#F2AB6D", `6KR_infected` = "#B4520C")
GROUP_LABELS <- c(wt_mock = "wt mock", wt_infected = "wt infected",
                  `6KR_mock` = "6KR mock", `6KR_infected` = "6KR infected")
DIRECTION_COLOURS <- c(up = "#B2182B", down = "#2166AC", ns = "#BDBDBD")

######### Result tables ##########
counts_raw <- readRDS(results_path("intermediate", "counts_raw.rds"))
sample_table <- fread(results_path("sample_table.csv"))
qc <- fread(results_path("qc_sample_summary.csv"))
pca_var <- fread(results_path("qc_pca_variance.csv"))$percent_variance

contrasts <- c("infection_wt", "infection_6KR", "genotype_mock", "genotype_infected",
               "interaction")
deg <- lapply(setNames(contrasts, contrasts), function(x) {
  fread(results_path("deg", paste0(x, "_all.csv")))
})
deg_counts <- fread(results_path("deg", "deg_counts.csv"))
dc <- function(contrast, column) deg_counts[[column]][match(contrast, deg_counts$contrast)]
n_tested <- nrow(deg$infection_wt)

ora <- fread(results_path("enrichment", "enrichment_results.csv"))
ora_counts <- fread(results_path("enrichment", "enrichment_term_counts.csv"))

rf <- fread(results_path("random_forest", "loocv_summary_metrics.csv"))
m1 <- rf[model == "group_lrt"]
m2 <- rf[model == "factor_balanced"]
rf_sens <- fread(results_path("random_forest", "loocv_n_features_sensitivity.csv"))
sens <- function(model_name, n_genes) rf_sens[model == model_name & n_features == n_genes, correct]
rf_pred <- fread(results_path("random_forest", "loocv_predictions.csv"))
importance <- fread(results_path("random_forest", "final_model_importance.csv"))
stability <- fread(results_path("random_forest", "gene_selection_stability.csv"))
candidates <- fread(results_path("biomarkers", "candidate_genes.csv"))

######### Plot helpers ##########
# Interactive volcano plot: one point per tested gene, gene name and numbers on hover.
# WebGL keeps 13,000 points responsive in the browser.
plot_volcano <- function(contrast, x_title) {
  tab <- deg[[contrast]][!is.na(padj)]
  tab[, label := fifelse(symbol == "", gene_id, symbol)]
  tab[, neg_log10_padj := -log10(pmax(padj, 1e-300))]
  tab[, hover := sprintf("<b>%s</b><br>log2 FC (shrunken): %.2f<br>adjusted p: %.2g",
                         label, log2FC_shrunk, padj)]
  tab[, direction := factor(direction, levels = c("up", "down", "ns"))]
  plot_ly(tab, x = ~log2FC_shrunk, y = ~neg_log10_padj, color = ~direction,
          colors = DIRECTION_COLOURS, text = ~hover, hoverinfo = "text",
          type = "scattergl", mode = "markers", marker = list(size = 5, opacity = 0.7)) |>
    layout(xaxis = list(title = x_title, zeroline = FALSE),
           yaxis = list(title = "-log10 adjusted p-value"),
           legend = list(title = list(text = "Direction")),
           shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                              y0 = -log10(0.05), y1 = -log10(0.05),
                              line = list(dash = "dash", width = 1, color = "#555555")))) |>
    config(displaylogo = FALSE)
}
