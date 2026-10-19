# Input data

The real input is not tracked in git. Place these files here before running the pipeline:

| File | Source | Used by |
|---|---|---|
| `GSE133329_RAW.tar` | GEO series GSE133329, "Download family → Supplementary file" (<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133329>) | `scripts/01_build_count_matrix.R` |
| `GSE133329_family.soft` | same GEO page, SOFT formatted family file (unpacked from `GSE133329_family.soft.gz`) | `scripts/01_build_count_matrix.R` (sample metadata) |
| `reference/course_count_table.csv` | count table handed out in the practical course | cross-check only; script 01 skips the check if the file is missing |

md5 checksums and sizes are in `CHECKSUMS.txt`; check them from inside `data/` with
`md5sum -c CHECKSUMS.txt`.

## Format of the real data

- `GSE133329_RAW.tar` holds 12 gzipped featureCounts tables, one per sample
  (`GSM3905824_wt_TR28_mock_I.txt.gz` … `GSM3905835_TR28_6KR_MOI10_III.txt.gz`). Each is
  tab-separated with the columns `Identifier` (Ensembl gene ID without version, GRCh38 release
  91) and `matchCounts` (integer read count); 21,386 genes in the same order in every file.
- `GSE133329_family.soft` provides title, genotype and infection per GSM sample.

Script 01 turns these into a gene × sample count matrix (`results/intermediate/counts_raw.rds`)
and a sample table (`results/sample_table.csv`).

## Example files

`example_counts.csv` and `example_sample_table.csv` show the layout of these two intermediate
tables with **simulated** data:

- `example_counts.csv`: column `gene_id` (Ensembl gene ID) plus one integer count column per
  sample; 200 real gene IDs with negative-binomial random counts, no biological signal.
- `example_sample_table.csv`: one row per sample with `sample`, `gsm`, `title`, `genotype`
  (`wt`, `6KR`), `infection` (`mock`, `infected`), `group` and `replicate`; identical to the
  real sample table.

The pipeline reads the GEO files directly and is not wired to run on the example files; they
document the expected format.
