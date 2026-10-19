# Build the public report PDF for the website: render report.qmd with -P public:true (no
# matriculation number, no supervisor names) and convert the docx to PDF with Microsoft Word.
# Run from website/:  powershell -File make_report_pdf.ps1
# Output: assets/influenza_rf_biomarkers_report.pdf

$ErrorActionPreference = "Stop"
$website = $PSScriptRoot
$reportDir = Join-Path $website "..\report"
$docx = Join-Path $reportDir "influenza_rf_biomarkers_report_public.docx"
$pdf = Join-Path $website "assets\influenza_rf_biomarkers_report.pdf"

Push-Location $reportDir
try {
  quarto render report.qmd -P public:true -o influenza_rf_biomarkers_report_public.docx
  if ($LASTEXITCODE -ne 0) { throw "quarto render failed" }
} finally {
  Pop-Location
}

New-Item -ItemType Directory -Force (Split-Path $pdf) | Out-Null
$word = New-Object -ComObject Word.Application
try {
  $document = $word.Documents.Open((Resolve-Path $docx).Path, $false, $true)
  $document.SaveAs2([System.IO.Path]::GetFullPath($pdf), 17)  # 17 = wdFormatPDF
  $document.Close(0)
} finally {
  $word.Quit()
}
Write-Output "Written $pdf"
