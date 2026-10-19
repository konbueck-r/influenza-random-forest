# influenza-random-forest

So this is a project from my bioinformatics practical course at the Friedrich Schiller University
Jena (spring 2025), where I looked at which genes change when lung epithelial cells get infected
with influenza A, and whether a Random Forest can tell the different groups apart only from their
gene expression. After the course I rebuilt the whole analysis as a proper R pipeline, because I
wanted to be sure that every number is right and that also someone else can run it again.

**If you just want to see what came out, have a look at the website:**
<https://konbueck-r.github.io/influenza-random-forest/>. There are the main results with
interactive plots, and the full report as PDF.

## What I did

The data come from a published experiment (GEO
[GSE133329](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133329), Schmidt et al., PNAS
2019): A549 cells with either the normal host protein TRIM28 (wt) or a mutant that cannot be
modified by SUMO anymore (6KR), mock treated or infected with influenza A, three replicates each,
so 12 samples in total.

- **Differential expression with DESeq2**, with one model for infection, genotype and the
  interaction of both, so I could ask which genes respond to infection, which to the genotype, and
  which respond to infection differently depending on the genotype.
- **Over-representation analysis with g:Profiler**, to see which processes are behind these gene
  lists.
- **A Random Forest classifier for the four groups**, evaluated with nested leave-one-out
  cross-validation. The important part here is that the gene selection happens inside every fold
  on the 11 training samples only, because otherwise the test sample would sneak into the gene
  selection and the accuracy would look better than it is.

The short version of the results: infection changes thousands of genes and is easy to detect, the
genotype is much harder. My planned model got 8 of 12 samples right and mixed up only the
genotypes, because its gene selection was dominated by infection genes. A model that selects
genotype and infection genes separately got 12 of 12, but I chose that rule after I saw the first
result, so I'd say this number is optimistic and would need an independent dataset to be
confirmed. The details and all limitations are in the report and on the website.

## How to run it

You need R 4.5.3 with Bioconductor 3.22. The package versions are in `renv.lock`, so
`renv::restore()` should give you the same setup (the project does not activate renv by itself).
The enrichment script also needs internet, because g:Profiler is a web service, and since they
update their databases from time to time, a later run can give slightly different terms (the
version I used is written to `results/enrichment/gprofiler_version.txt`).

The raw data are not in this repository. Where to download them and how they should look is
described in [`data/README.md`](data/README.md). Then, from the project folder:

```bash
Rscript run_all.R
```

This runs the scripts `01` to `07` in order and takes about 5 minutes on my laptop, most of it is
the Random Forest cross-validation. For the report you also need Quarto, and then you can build it
from inside `report/` with:

```bash
quarto render report.qmd
```

Every number in the report and on the website is read from `results/` when they are built, so
nothing is typed in by hand. How to build the website is explained in
[`website/README.md`](website/README.md).

## What is where

| Script | What it does | Where the output goes |
|---|---|---|
| `01_build_count_matrix.R` | count matrix and sample table from the GEO files | `results/intermediate/`, `results/sample_table.csv` |
| `02_qc_exploration.R` | gene filter, rlog, PCA, sample distances | `results/`, `figures/qc/` |
| `03_differential_expression.R` | DESeq2 interaction model, five contrasts | `results/deg/`, `figures/deg/` |
| `04_enrichment.R` | over-representation analysis with g:Profiler | `results/enrichment/`, `figures/enrichment/` |
| `05_random_forest.R` | four-class Random Forest, nested leave-one-out CV | `results/random_forest/`, `figures/random_forest/` |
| `06_biomarker_summary.R` | candidate biomarker table | `results/biomarkers/`, `figures/biomarkers/` |
| `07_report_figures.R` | figures for the report and the website | `figures/report/` |

The helper functions (gene filter, colours, plot theme) are in `functions.R`, the report source is
in `report/` and the website source in `website/`.

## Licence

The code is under the MIT licence (`LICENSE`), so feel free to reuse it. The photos in
`website/images/` are not covered by that and please don't reuse them. The RNA-seq data belong to
GEO series GSE133329 and are not part of this repository.
