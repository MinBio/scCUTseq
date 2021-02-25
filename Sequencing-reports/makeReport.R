# Render markdown file with custom params

params = list(basedir = "/mnt/AchTeraD/data/",
              run = "BICRO253",
              extra_id = "",
              run_description = "scCUTseq on Mouse brain E11.5",
              libraries = c("NZ184"),
              library_description = c("scCUTseq on Mouse brain E11.5"),
              binsizes = c("500000"))


rmarkdown::render(
  input = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/Sequencing-reports/sequencing-qc.Rmd",
  output_file = paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/sequencing-reports/", 
                       params$run, params$extra_id, "-sequence_report.html"),
  params = params
)
