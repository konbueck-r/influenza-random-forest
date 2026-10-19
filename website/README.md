# Project website

Multi-page Quarto website that presents the analysis: overview, data and QC, differential
expression, Random Forest, methods and reproducibility. Every number is read from `../results/`
when the site is built, as in the report.

## Build

Requires the pipeline outputs (`Rscript run_all.R` in the project root), Quarto, and Microsoft Word
for the PDF conversion. From this folder:

```bash
powershell -File make_report_pdf.ps1   # public report PDF -> assets/
quarto render                          # site -> _site/
```

`quarto render` first runs `copy_assets.R`, which copies the report figures into `assets/`.
`_site/` and `assets/` are build output and not tracked.

The public report PDF is rendered with `-P public:true`: it shows only the author's name and the
course, not the matriculation number or the supervisors, which are kept in the untracked file
`../report/private_details.txt` for the submission version.

## Publish

`quarto publish gh-pages` renders the site and pushes it to the `gh-pages` branch, which GitHub
Pages serves.
