# Run the whole pipeline in order from the project root and record the session.
# Usage: Rscript run_all.R

scripts <- c(
  "scripts/01_build_count_matrix.R",
  "scripts/02_qc_exploration.R",
  "scripts/03_differential_expression.R",
  "scripts/04_enrichment.R",
  "scripts/05_random_forest.R",
  "scripts/06_biomarker_summary.R",
  "scripts/07_report_figures.R"
)

for (script in scripts) {
  cat("\n=====", script, "=====\n")
  source(script, echo = FALSE, local = new.env())
}

writeLines(capture.output(sessionInfo()), file.path("results", "session_info.txt"))
cat("\nPipeline finished; session info written to results/session_info.txt\n")
