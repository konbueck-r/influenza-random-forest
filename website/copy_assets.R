# Copy the figures the website shows from the pipeline output into website/assets/.
# Runs as Quarto pre-render script with website/ as working directory.
# Input:  ../figures/ (written by run_all.R)
# Output: assets/figures/*.png

######### Configuration ##########
FIGURES <- c(
  file.path("report", sprintf("figure_%d.png", 1:5)),
  file.path("report", c("supplementary_1.png", "supplementary_2.png", "supplementary_5.png")),
  file.path("deg", "heatmap_interaction_top_genes.png"),
  file.path("random_forest", "final_model_importance_group_lrt.png")
)
SOURCE_DIR <- file.path("..", "figures")
TARGET_DIR <- file.path("assets", "figures")
REPORT_PDF <- file.path("assets", "influenza_rf_biomarkers_report.pdf")

######### Copy ##########
dir.create(TARGET_DIR, recursive = TRUE, showWarnings = FALSE)
sources <- file.path(SOURCE_DIR, FIGURES)
missing <- sources[!file.exists(sources)]
if (length(missing) > 0) {
  stop("Figures missing, run run_all.R first: ", paste(missing, collapse = ", "))
}
file.copy(sources, file.path(TARGET_DIR, basename(FIGURES)), overwrite = TRUE)

if (!file.exists(REPORT_PDF)) {
  warning("Public report PDF missing (", REPORT_PDF, "); run make_report_pdf.ps1 first.")
}
cat("Copied", length(sources), "figures to", TARGET_DIR, "\n")
