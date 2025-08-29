#' Render R Markdown report with organized output structure
#'
#' Renders an R Markdown report to HTML with automatic directory organization
#' and naming conventions. The function creates a dated output directory structure
#' and passes parameters for title, cache path, and rerun options.
#'
#' @param report Character string specifying the path to the R Markdown file to render.
#'   Expected to be in format \code{reports/date_name/report.Rmd} where the subdirectory
#'   name will be used for organizing outputs.
#' @param rerun Logical indicating whether to rerun cached chunks. Default is FALSE.
#'   When TRUE, forces re-execution of cached code chunks.
#'
#' @return Invisibly returns the path to the rendered HTML file. The function is
#'   called primarily for its side effect of creating the rendered report.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Extracts subdirectory name from the report path (second path component)
#'   \item Removes prefix from subdirectory name to create a clean report name
#'   \item Creates output directory: \code{reports/date_name/YYYYMMDD/}
#'   \item Generates output filename: \code{YYYYMMDD_name.html}
#'   \item Sets up cache directory: \code{reports/date_name/YYYYMMDD/YYYYMMDD_name_cache/}
#'   \item Passes parameters to the R Markdown document: title, cache_path, rerun
#' }
#'
#' @examples
#' \dontrun{
#' # Render a report located at "reports/20250820_mutation_analysis/report.Rmd"
#' render_report("reports/20250820_mutation_analysis/report.Rmd")
#' 
#' # Force rerun of cached chunks
#' render_report("reports/20250820_mutation_analysis/report.Rmd", rerun = TRUE)
#' 
#' # This would create:
#' # - Output directory: reports/20250820_mutation_analysis/20250828/
#' # - HTML file: 20250828_mutation_analysis.html
#' # - Cache directory: reports/20250820_mutation_analysis/20250828/20250828_mutation_analysis_cache/
#' }
#'
#' @importFrom rmarkdown render
#' @export
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