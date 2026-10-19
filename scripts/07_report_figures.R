# Compose the main and supplementary report figures from the plot objects of scripts 02-06.
# Input:  results/intermediate/plots/plots_0*.rds
# Output: figures/report/figure_1-5.png and supplementary_*.png (A4 text width, 300 dpi)

######### Libraries ##########
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("functions.R")

######### Configuration ##########
TEXT_WIDTH_CM <- 16
PLOT_DIR <- file.path("results", "intermediate", "plots")
FIG_DIR <- file.path("figures", "report")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

######### Data loading ##########
plots_qc <- readRDS(file.path(PLOT_DIR, "plots_02_qc.rds"))
plots_deg <- readRDS(file.path(PLOT_DIR, "plots_03_deg.rds"))
plots_ora <- readRDS(file.path(PLOT_DIR, "plots_04_enrichment.rds"))
plots_rf <- readRDS(file.path(PLOT_DIR, "plots_05_rf.rds"))
plots_bio <- readRDS(file.path(PLOT_DIR, "plots_06_biomarkers.rds"))

######### Helper functions ##########
# Lower-case bold panel letters as in Nature figures.
tag_panels <- function(figure) {
  figure + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial"))
}

# Panel titles in the composites are short; the figure captions carry the details.
short_title <- function(plot, title) {
  plot + labs(title = title) + theme(plot.title = element_text(size = 9, face = "bold"))
}

# Volcano plots without their legend: the colours are explained in the caption and the
# up/down counts are shown in the bar chart.
volcano <- function(contrast, title) {
  short_title(plots_deg$volcanoes[[contrast]], title) +
    labs(x = "log2 fold change (shrunken)") +
    theme(legend.position = "none")
}

confusion <- function(model_name, model_label) {
  plot <- plots_rf$confusion[[model_name]]
  n_correct <- with(plot$data, sum(N[as.character(predicted) == as.character(truth)]))
  short_title(plot, sprintf("%s: %d of %d correct", model_label, n_correct, sum(plot$data$N)))
}

######### Figure 1: design and sample structure ##########
design <- data.table(
  genotype = factor(c("wt", "wt", "6KR", "6KR"), levels = c("6KR", "wt")),
  infection = factor(c("mock", "infected", "mock", "infected"), levels = c("mock", "infected")),
  group = c("wt_mock", "wt_infected", "6KR_mock", "6KR_infected")
)
design[, label := sprintf("%s\n%s\nn = 3", sub("_.*", "", group), sub(".*_", "", group))]
p_design <- ggplot(design, aes(infection, genotype, fill = group)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = infection == "infected"), size = 2.7, family = "Arial",
            lineheight = 0.9) +
  scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_discrete(labels = c(mock = "mock", infected = "IAV WSN\nMOI 10, 6 hpi")) +
  scale_y_discrete(labels = c(wt = "TRIM28\nwt", `6KR` = "TRIM28\n6KR")) +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_project() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank())
p_design <- short_title(p_design, "Design")

p_pca <- short_title(plots_qc$pca[[1]], "PCA, 500 most variable genes") +
  theme(legend.position = "right")

figure_1 <- (p_design + p_pca + plot_layout(widths = c(1, 1.7))) /
  wrap_elements(full = plots_qc$distance_heatmap) +
  plot_layout(heights = c(1, 1.8))
save_plot(tag_panels(figure_1), file.path(FIG_DIR, "figure_1.png"), TEXT_WIDTH_CM, 15)

######### Figure 2: differential expression ##########
figure_2 <- free(short_title(plots_deg$deg_counts, "DEGs per contrast")) /
  (volcano("infection_wt", "Infected vs mock, wt cells") |
     volcano("interaction", "Interaction (6KR vs wt)")) +
  plot_layout(heights = c(0.75, 1.6))
save_plot(tag_panels(figure_2), file.path(FIG_DIR, "figure_2.png"), TEXT_WIDTH_CM, 14.5)

######### Figure 3: enrichment ##########
figure_3 <- free(short_title(plots_ora$infection_wt_up, "Up after infection, wt cells")) /
  free(short_title(plots_ora$interaction_up, "Interaction up (higher in 6KR)")) /
  free(short_title(plots_ora$genotype_mock_up, "Up in 6KR vs wt, mock cells"))
save_plot(tag_panels(figure_3), file.path(FIG_DIR, "figure_3.png"), TEXT_WIDTH_CM, 18)

######### Figure 4: Random Forest ##########
figure_4 <- (confusion("group_lrt", "Model 1") | confusion("factor_balanced", "Model 2")) /
  free(short_title(plots_rf$sensitivity, "Correct predictions by number of genes") +
         labs(y = "Correct (of 12)") +
         theme(legend.position = "right", legend.direction = "vertical")) +
  plot_layout(heights = c(1.3, 1))
save_plot(tag_panels(figure_4), file.path(FIG_DIR, "figure_4.png"), TEXT_WIDTH_CM, 12)

######### Figure 5: candidate biomarkers ##########
figure_5 <- free(short_title(plots_rf$importance$factor_balanced,
                             "Model 2: top 20 genes by importance")) /
  free(short_title(plots_bio$candidate_expression, "Candidate genes")) +
  plot_layout(heights = c(1, 1.4))
save_plot(tag_panels(figure_5), file.path(FIG_DIR, "figure_5.png"), TEXT_WIDTH_CM, 19.5)

######### Supplementary figures ##########
supp_1 <- (plots_qc$library_sizes + theme(legend.position = "none")) / plots_qc$rlog_boxplot
save_plot(tag_panels(supp_1), file.path(FIG_DIR, "supplementary_1.png"), TEXT_WIDTH_CM, 16)

supp_2 <- plots_qc$dendrogram / plots_qc$scree + plot_layout(heights = c(1.4, 1))
save_plot(tag_panels(supp_2), file.path(FIG_DIR, "supplementary_2.png"), TEXT_WIDTH_CM, 17)

supp_3 <- (volcano("infection_6KR", "Infected vs mock, 6KR cells") |
             volcano("genotype_mock", "6KR vs wt, mock cells")) /
  (volcano("genotype_infected", "6KR vs wt, infected cells") | plot_spacer())
save_plot(tag_panels(supp_3), file.path(FIG_DIR, "supplementary_3.png"), TEXT_WIDTH_CM, 18)

supp_5 <- free(short_title(plots_ora$infection_wt_down, "Down after infection, wt cells")) /
  free(short_title(plots_ora$genotype_mock_down, "Down in 6KR vs wt, mock cells")) /
  free(short_title(plots_ora$interaction_down, "Interaction down (lower in 6KR)"))
save_plot(tag_panels(supp_5), file.path(FIG_DIR, "supplementary_5.png"), TEXT_WIDTH_CM, 19)

######### Output and summary ##########
cat("Written:", paste(list.files(FIG_DIR), collapse = ", "), "to", FIG_DIR, "\n")
