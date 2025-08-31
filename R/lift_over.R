#' Lift over genomic coordinates between genome builds
#'
#' Converts genomic coordinates from one genome build to another using a liftOver
#' chain file. The function handles both single position variants and genomic ranges,
#' automatically detecting the input format and maintaining all original data.
#'
#' @param variants A data frame containing variant/genomic position data with columns:
#'   \describe{
#'     \item{chr}{Chromosome}
#'     \item{pos}{Position (for single positions)}
#'     \item{start, end}{Start and end positions (for genomic ranges)}
#'     \item{...}{Additional columns will be preserved in the output}
#'   }
#'   Must contain either \code{pos} column OR both \code{start} and \code{end} columns.
#' @param lift_over_chain Path to the liftOver chain file for coordinate conversion
#'   between genome builds (e.g., hg19 to hg38).
#'
#' @return A data frame with lifted coordinates and all original columns preserved.
#'   For range data: \code{start} and \code{end} columns contain the new coordinates.
#'   For position data: \code{pos} column contains the new coordinates.
#'   Only successfully lifted coordinates are returned.
#'
#' @details 
#' The function performs the following operations:
#' \itemize{
#'   \item Automatically detects input format (single positions vs ranges)
#'   \item Converts data to GenomicRanges object for liftOver processing
#'   \item Applies coordinate conversion using the specified chain file
#'   \item Filters to successfully lifted coordinates only
#'   \item Rejoins lifted coordinates with original variant data
#'   \item Returns data with updated coordinates and all original columns
#' }
#'
#' The function uses an inner join, so only variants that successfully lift over
#' will be included in the output. Variants that fail to lift over are excluded.
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' library(rtracklayer)
#' library(dplyr)
#' 
#' # Example with single positions
#' variants_pos <- data.frame(
#'   chr = c("chr1", "chr2", "chr3"),
#'   pos = c(1000000, 2000000, 3000000),
#'   ref = c("A", "G", "C"),
#'   alt = c("T", "C", "A")
#' )
#' 
#' # Lift over from hg19 to hg38
#' lifted_pos <- lift_over(variants_pos, "hg19ToHg38.over.chain")
#'
#' # Example with genomic ranges
#' variants_range <- data.frame(
#'   chr = c("chr1", "chr2"),
#'   start = c(1000000, 2000000),
#'   end = c(1000100, 2000100),
#'   gene = c("GENE1", "GENE2")
#' )
#'
#' lifted_range <- lift_over(variants_range, "hg19ToHg38.over.chain")
#' }
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer liftOver import.chain
#' @importFrom tibble as_tibble
#' @importFrom dplyr transmute inner_join mutate select distinct
#' @export
lift_over <- function(variants, lift_over_chain) {
  # check position column type
  if (all(c("start", "end") %in% colnames(variants))) {
    pos_type <- "range"
  } else if ("pos" %in% colnames(variants)) {
    pos_type <- "pos"
  } else {
    stop("Input `variants` must contain either 'start' and 'end' columns or 'pos' column.")
  }

  # convert to granges object
  gr <-
    GRanges(
      seqnames = variants$chr,
      ranges = IRanges(
        start = ifelse(pos_type == "range", variants$start, variants$pos),
        end = ifelse(pos_type == "range", variants$end, variants$pos)))

  # set old positions as metadata
  if (pos_type == "range") {
    gr$start_old <- variants$start
    gr$end_old <- variants$end
  } else {
    gr$pos <- variants$pos
  }

  # liftover
  lifted <-
    gr %>%
    rtracklayer::liftOver(rtracklayer::import.chain(lift_over_chain)) %>%
    unlist() %>%
    tibble::as_tibble() %>%
    dplyr::distinct()

  # add all variant columns back
  if (pos_type == "range") {
    lifted %>%
      dplyr::transmute(chr = seqnames, start_new = start, end_new = end,
                       start = start_old, end = end_old) %>%
      dplyr::inner_join(variants) %>%
      dplyr::mutate(start = start_new, end = end_new) %>%
      dplyr::select(-start_new, -end_new)
  } else {
    lifted %>%
      dplyr::transmute(chr = seqnames, pos_new = start,
                       pos = pos) %>%
      dplyr::inner_join(variants) %>%
      dplyr::mutate(pos = pos_new) %>%
      dplyr::select(-pos_new)
  }
}