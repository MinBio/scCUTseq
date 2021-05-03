# Render markdown file with custom params

params = list(basedir = "/mnt/AchTeraD/data/",
              run = "BICRO284",
              extra_id = "",
              run_description = "scCUTseq on Prostate (P6) and BRCA1",
              libraries = c("NZ231", "NZ249", "MS159"),
              library_description = c("scCUTseq on Prostate P6L3C3", "scCUTseq on Breast cancer", "scCUTseq on Prostate P6L3C1"),
              binsizes = c("500000"))


rmarkdown::render(
  input = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/Sequencing-reports/sequencing-qc.Rmd",
  output_file = paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/sequencing-reports/", 
                       params$run, params$extra_id, "-sequence_report.html"),
  params = params
)
