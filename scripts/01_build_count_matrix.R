# Build the raw count matrix and the sample table of GEO GSE133329.
# Input:  data/GSE133329_RAW.tar (12 featureCounts tables), data/GSE133329_family.soft,
#         data/reference/course_count_table.csv (supervisor table, optional cross-check)
# Output: results/intermediate/counts_raw.rds, results/sample_table.csv

######### Libraries ##########
library(data.table)

######### Configuration ##########
RAW_TAR <- file.path("data", "GSE133329_RAW.tar")
SOFT_FILE <- file.path("data", "GSE133329_family.soft")
SUPERVISOR_TABLE <- file.path("data", "reference", "course_count_table.csv")
EXTRACT_DIR <- file.path("results", "intermediate", "geo_raw")
COUNTS_OUT <- file.path("results", "intermediate", "counts_raw.rds")
SAMPLE_TABLE_OUT <- file.path("results", "sample_table.csv")

dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

######### Sample table from the GEO series file ##########
# Each sample block in the SOFT file starts with "^SAMPLE = GSM..."; title and characteristics
# lines follow within the block.
soft_lines <- readLines(SOFT_FILE, encoding = "UTF-8")
sample_starts <- grep("^\\^SAMPLE = ", soft_lines)
block_ends <- c(sample_starts[-1] - 1, length(soft_lines))

read_soft_field <- function(block, prefix) {
  hit <- grep(prefix, block, fixed = TRUE, value = TRUE)
  if (length(hit) != 1) stop("Expected exactly one '", prefix, "' line per sample block")
  sub(prefix, "", hit, fixed = TRUE)
}

sample_table <- rbindlist(lapply(seq_along(sample_starts), function(i) {
  block <- soft_lines[sample_starts[i]:block_ends[i]]
  data.table(
    gsm = read_soft_field(block, "^SAMPLE = "),
    title = read_soft_field(block, "!Sample_title = "),
    genotype_text = read_soft_field(block, "!Sample_characteristics_ch1 = genotype/variation: "),
    infection_text = read_soft_field(block, "!Sample_characteristics_ch1 = infection: ")
  )
}))

# Translate the GEO free text into short factor labels; stop on anything unexpected.
sample_table[, genotype := fcase(
  grepl("TRIM28 wt$", genotype_text), "wt",
  grepl("TRIM28 6KR", genotype_text), "6KR"
)]
sample_table[, infection := fcase(
  grepl("mock", infection_text), "mock",
  grepl("Influenza A", infection_text), "infected"
)]
sample_table[, replicate := match(sub(".*_", "", title), c("I", "II", "III"))]
if (anyNA(sample_table[, .(genotype, infection, replicate)])) {
  stop("Could not parse genotype, infection or replicate for some samples")
}

sample_table[, group := paste(genotype, infection, sep = "_")]
sample_table[, sample := paste(group, replicate, sep = "_")]
sample_table[, genotype := factor(genotype, levels = c("wt", "6KR"))]
sample_table[, infection := factor(infection, levels = c("mock", "infected"))]
sample_table[, group := factor(
  group, levels = c("wt_mock", "wt_infected", "6KR_mock", "6KR_infected")
)]
setorder(sample_table, genotype, infection, replicate)

######### Count matrix from the per-sample files ##########
untar(RAW_TAR, exdir = EXTRACT_DIR)
count_files <- list.files(EXTRACT_DIR, pattern = "^GSM\\d+_.*\\.txt\\.gz$", full.names = TRUE)
names(count_files) <- sub("_.*", "", basename(count_files))

# Every GEO sample needs exactly one count file and vice versa.
if (!setequal(names(count_files), sample_table$gsm) || length(count_files) != nrow(sample_table)) {
  stop("Count files and GEO samples do not match one to one")
}

count_tables <- lapply(count_files[sample_table$gsm], fread)
for (gsm in names(count_tables)) {
  if (!identical(names(count_tables[[gsm]]), c("Identifier", "matchCounts"))) {
    stop("Unexpected columns in the count file of ", gsm)
  }
}

# Bind by position only after confirming identical gene IDs in identical order.
gene_ids <- count_tables[[1]]$Identifier
same_genes <- vapply(count_tables, function(tab) identical(tab$Identifier, gene_ids), logical(1))
if (!all(same_genes)) stop("Gene IDs differ between count files: ", names(which(!same_genes)))
if (anyDuplicated(gene_ids)) stop("Duplicated gene IDs in the count files")

counts_raw <- sapply(count_tables, function(tab) tab$matchCounts)
dimnames(counts_raw) <- list(gene_ids, sample_table$sample)

if (anyNA(counts_raw) || any(counts_raw < 0) || any(counts_raw != round(counts_raw))) {
  stop("Counts must be non-negative integers without NA")
}
storage.mode(counts_raw) <- "integer"

######### Cross-check against the supervisor table ##########
# The course table is not public; without it the check is skipped and reported as NA.
cross_check <- list(same_genes_same_order = NA, same_samples = NA, identical_counts = NA)
if (file.exists(SUPERVISOR_TABLE)) {
  supervisor <- fread(SUPERVISOR_TABLE, sep = ";")
  supervisor_counts <- as.matrix(supervisor[, -1])
  rownames(supervisor_counts) <- supervisor[[1]]
  colnames(supervisor_counts) <- sub("_.*", "", colnames(supervisor_counts))

  cross_check$same_genes_same_order <- identical(rownames(supervisor_counts), gene_ids)
  cross_check$same_samples <- setequal(colnames(supervisor_counts), sample_table$gsm)
  if (cross_check$same_genes_same_order && cross_check$same_samples) {
    supervisor_counts <- supervisor_counts[, sample_table$gsm]
    cross_check$identical_counts <- all(supervisor_counts == counts_raw)
  }
}

######### Output and summary ##########
saveRDS(counts_raw, COUNTS_OUT)
fwrite(
  sample_table[, .(sample, gsm, title, genotype, infection, group, replicate)],
  SAMPLE_TABLE_OUT
)

cat("Samples:", ncol(counts_raw), "| genes:", nrow(counts_raw), "\n")
cat("Genes with zero counts in all samples:", sum(rowSums(counts_raw) == 0), "\n")
cat("Library sizes (million reads):\n")
print(round(colSums(counts_raw) / 1e6, 2))
cat("\nCross-check with the supervisor table:\n")
print(unlist(cross_check))
cat("\nWritten:", COUNTS_OUT, "and", SAMPLE_TABLE_OUT, "\n")
