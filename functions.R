# Helpers shared across the pipeline scripts: gene filter, group colours, plot theme, plot saving.

# One colour per group: hue marks the genotype (blue = wt, orange = 6KR), shade marks infection
# (light = mock, dark = infected), so the 2 x 2 design can be read off every figure.
GROUP_COLOURS <- c(
  wt_mock = "#8DB6DD",
  wt_infected = "#1F5A96",
  `6KR_mock` = "#F2AB6D",
  `6KR_infected` = "#B4520C"
)

GROUP_LABELS <- c(
  wt_mock = "wt mock",
  wt_infected = "wt infected",
  `6KR_mock` = "6KR mock",
  `6KR_infected` = "6KR infected"
)

# Keep genes with at least `min_count` raw counts in at least `min_samples` samples. The default
# of 3 samples is the smallest group, so genes switched on in a single condition are kept.
filter_expressed_genes <- function(counts, min_count = 10, min_samples = 3) {
  rowSums(counts >= min_count) >= min_samples
}

theme_project <- function(base_size = 9) {
  ggplot2::theme_bw(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )
}

# Save a ggplot as PNG with Arial (ragg device) at a fixed size in cm.
save_plot <- function(plot, path, width_cm, height_cm, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, device = ragg::agg_png,
                  width = width_cm, height = height_cm, units = "cm", dpi = dpi)
}
