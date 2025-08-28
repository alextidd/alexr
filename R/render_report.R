render_report <- function(report, rerun = FALSE) {
  subdir <- strsplit(report, "/")[[1]][2]
  name <- sub("^[^_]+_", "", subdir)
  date <- format(Sys.Date(), "%Y%m%d")
  outdir <- file.path("reports", subdir, date)
  basename <- paste0(date, "_", name)
  rmarkdown::render(
    report,
    output_file = paste0(basename, ".html"),
    output_dir = outdir,
    params = list(title = name,
                  cache_path = file.path(outdir, paste0(basename, "_cache/")),
                  rerun = rerun))
}