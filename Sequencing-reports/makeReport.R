# Render markdown file with custom params

params = list(basedir = "/mnt/AchTeraD/data/",
              run = "scCUTseq_turin1",
              extra_id = "",
              run_description = "scCUTseq performed in Turin",
              libraries = c("EB1", "EB2"),
              library_description = c("scCUTseq performed in Turin lib1", "scCUTseq performed in Turin lib2"),
              binsizes = c("1000000", "500000"))


rmarkdown::render(
  input = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/Sequencing-reports/sequencing-qc.Rmd",
  output_file = paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/sequencing-reports/", 
                       params$run, params$extra_id, "-sequence_report.html"),
  params = params
)
  