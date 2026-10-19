# Over-representation analysis (ORA) of the DEG lists with g:Profiler.
# Input:  results/deg/<contrast>_all.csv
# Output: results/enrichment/*.csv, figures/enrichment/*.png
# Needs internet access (g:Profiler web API).

######### Libraries ##########
suppressPackageStartupMessages({
  library(data.table)
  library(gprofiler2)
  library(ggplot2)
})
source("functions.R")

######### Configuration ##########
ALPHA <- 0.05          # BH-adjusted p-value cutoff for DEGs and terms
MIN_ABS_LFC <- 1       # shrunken log2 fold-change cutoff for the large DEG lists
MIN_TERM_SIZE <- 10    # term size limits as in the course, to drop very broad and tiny terms
MAX_TERM_SIZE <- 1000
N_TERMS_PER_SOURCE <- 5  # plots show the best terms of each source, so one source cannot fill them
SOURCES <- c("GO:BP", "KEGG", "REAC")
DEG_DIR <- file.path("results", "deg")
OUT_DIR <- file.path("results", "enrichment")
FIG_DIR <- file.path("figures", "enrichment")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# The interaction list is small, so it gets no fold-change cutoff.
QUERY_SPECS <- data.table(
  contrast = c("infection_wt", "genotype_mock", "interaction"),
  min_abs_lfc = c(MIN_ABS_LFC, MIN_ABS_LFC, 0)
)
QUERY_TITLES <- c(
  infection_wt_up = "Up after infection (wt cells)",
  infection_wt_down = "Down after infection (wt cells)",
  genotype_mock_up = "Up in 6KR vs wt (mock cells)",
  genotype_mock_down = "Down in 6KR vs wt (mock cells)",
  interaction_up = "Interaction up (higher in 6KR)",
  interaction_down = "Interaction down (lower in 6KR)"
)

######### Data loading ##########
deg_tables <- lapply(setNames(QUERY_SPECS$contrast, QUERY_SPECS$contrast), function(contrast) {
  fread(file.path(DEG_DIR, paste0(contrast, "_all.csv")))
})

# Background: every gene that passed the expression filter and was tested.
background <- deg_tables$infection_wt$gene_id

######### Preprocessing ##########
queries <- list()
for (i in seq_len(nrow(QUERY_SPECS))) {
  contrast <- QUERY_SPECS$contrast[i]
  tab <- deg_tables[[contrast]][direction != "ns" & abs(log2FC_shrunk) > QUERY_SPECS$min_abs_lfc[i]]
  queries[[paste0(contrast, "_up")]] <- tab[direction == "up", gene_id]
  queries[[paste0(contrast, "_down")]] <- tab[direction == "down", gene_id]
}
query_sizes <- data.table(query = names(queries), n_genes = lengths(queries))

######### Analysis ##########
# One call with a named list gives one result block per query. Ordered = FALSE: plain ORA.
gost_result <- gost(
  query = queries,
  organism = "hsapiens",
  sources = SOURCES,
  domain_scope = "custom",
  custom_bg = background,
  correction_method = "fdr",   # Benjamini-Hochberg
  user_threshold = ALPHA,
  significant = TRUE,
  evcodes = TRUE,
  multi_query = FALSE
)
version_info <- get_version_info(organism = "hsapiens")

enrichment <- as.data.table(gost_result$result)
# g:Profiler's p_value is already the BH-adjusted value; rename it to say so.
setnames(enrichment, "p_value", "p_adj")
enrichment <- enrichment[term_size >= MIN_TERM_SIZE & term_size <= MAX_TERM_SIZE]
enrichment[, gene_ratio := intersection_size / query_size]
list_cols <- names(enrichment)[vapply(enrichment, is.list, logical(1))]
for (col in list_cols) {
  enrichment[, (col) := vapply(get(col), paste, character(1), collapse = ",")]
}
setorder(enrichment, query, p_adj)
enrichment_out <- enrichment[, .(query, source, term_id, term_name, term_size, query_size,
                                 intersection_size, gene_ratio, p_adj, intersection)]

term_counts <- enrichment[, .N, by = .(query, source)]
term_counts <- dcast(term_counts, query ~ source, value.var = "N", fill = 0)
term_counts <- merge(query_sizes, term_counts, by = "query", all.x = TRUE, sort = FALSE)

######### Visualization ##########
shorten <- function(x, width = 55) {
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

plot_query <- function(query_name) {
  top <- enrichment[query == query_name][, head(.SD, N_TERMS_PER_SOURCE), by = source]
  setorder(top, p_adj)
  top[, label := sprintf("%s [%s]", shorten(term_name), source)]
  top[, label := factor(label, levels = rev(unique(label)))]
  ggplot(top, aes(gene_ratio, label, size = intersection_size, colour = -log10(p_adj))) +
    geom_point() +
    scale_colour_gradient(low = "#9ECAE1", high = "#08306B", name = "-log10 adj. p") +
    # Gene counts are whole numbers, so the size legend shows integer breaks only.
    scale_size_continuous(range = c(1.5, 5), name = "Genes",
                          breaks = function(limits) unique(floor(pretty(limits)))) +
    labs(x = "Gene ratio (list genes in term / list size)", y = NULL,
         title = sprintf("%s (%d genes)", QUERY_TITLES[[query_name]],
                         length(queries[[query_name]]))) +
    theme_project()
}

######### Output and summary ##########
fwrite(enrichment_out, file.path(OUT_DIR, "enrichment_results.csv"))
fwrite(term_counts, file.path(OUT_DIR, "enrichment_term_counts.csv"))
writeLines(capture.output(str(version_info)), file.path(OUT_DIR, "gprofiler_version.txt"))

plotted <- character(0)
for (query_name in names(queries)) {
  if (enrichment[query == query_name, .N] == 0) next
  save_plot(plot_query(query_name), file.path(FIG_DIR, paste0("ora_", query_name, ".png")),
            17, 9.5)
  plotted <- c(plotted, query_name)
}
# Plot objects for the composite report figures (scripts/07_report_figures.R).
dir.create(file.path("results", "intermediate", "plots"), showWarnings = FALSE)
saveRDS(lapply(setNames(plotted, plotted), plot_query),
        file.path("results", "intermediate", "plots", "plots_04_enrichment.rds"))

cat("g:Profiler:", version_info$gprofiler_version, "| database:",
    version_info$biomart_version, "| background genes:", length(background), "\n\n")
cat("Significant terms per query (BH-adjusted p <", ALPHA, ", term size",
    MIN_TERM_SIZE, "-", MAX_TERM_SIZE, "):\n")
print(term_counts)
cat("\nQueries without a significant term (no plot):",
    paste(setdiff(names(queries), plotted), collapse = ", "), "\n")
for (query_name in plotted) {
  cat("\nTop terms,", query_name, ":\n")
  print(enrichment[query == query_name][1:min(.N, 8), .(source, term_name, intersection_size,
                                                         p_adj = signif(p_adj, 2))])
}
