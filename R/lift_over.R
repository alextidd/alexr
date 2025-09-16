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
#' @param keep_old_pos Logical indicating whether to keep the original coordinates
#'   in the output. Default is FALSE. When TRUE, original coordinates are preserved
#'   as \code{pos_old} (for positions) or \code{start_old}/\code{end_old} (for ranges).
#'
#' @return A data frame with lifted coordinates and all original columns preserved.
#'   For range data: \code{start} and \code{end} columns contain the new coordinates.
#'   For position data: \code{pos} column contains the new coordinates.
#'   When \code{keep_old_pos = TRUE}, original coordinates are preserved as additional columns.
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
#'   \item Optionally preserves original coordinates when \code{keep_old_pos = TRUE}
#' }
#'
#' The function uses an inner join, so only variants that successfully lift over
#' will be included in the output. Variants that fail to lift over are excluded.
#' 
#' When \code{keep_old_pos = TRUE}, the original coordinates are preserved as:
#' \itemize{
#'   \item \code{pos_old} for single position data
#'   \item \code{start_old} and \code{end_old} for range data
#' }
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
#' # Lift over and keep original positions
#' lifted_with_old <- lift_over(variants_pos, "hg19ToHg38.over.chain", keep_old_pos = TRUE)
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
lift_over <- function(variants, lift_over_chain, keep_old_pos = FALSE) {

  # check position column type, create granges object
  if (all(c("start", "end") %in% colnames(variants))) {

    pos_type <- "range"
    gr <-
      GRanges(seqnames = variants$chr,
              ranges = IRanges(start = variants$start, end = variants$end),
              start_old = variants$start, end_old = variants$end)

  } else if ("pos" %in% colnames(variants)) {

    pos_type <- "pos"
    gr <-
      GRanges(seqnames = variants$chr,
              ranges = IRanges(start = variants$pos, end = variants$pos),
              pos_old = variants$pos)

  } else {

    stop("`variants` must contain either 'start' and 'end' columns or 'pos' column.")

  }

  # liftover
  lifted_gr <-
    gr %>%
    rtracklayer::liftOver(rtracklayer::import.chain(lift_over_chain)) %>%
    unlist()

  # warn if variants lost
  if (length(gr) > length(lifted_gr)) {
    warning(paste0(length(gr) - length(lifted_gr), " / ", length(gr),
                   " variant(s) failed to lift over and will be removed."))
  }

  # convert to tibble
  lifted <-
    lifted_gr %>%
    tibble::as_tibble() %>%
    dplyr::distinct()

  # add all variant columns back, optionally remove old positions
  if (pos_type == "range") {

    lifted <-
      lifted %>%
      dplyr::transmute(chr = seqnames, start, end,
                       start_old, end_old) %>%
      dplyr::inner_join(variants %>% dplyr::rename(start_old = start,
                                                   end_old = end))
    if (keep_old_pos == FALSE) {
      lifted <- lifted %>% dplyr::select(-start_old, -end_old)
    }

  } else {

    lifted <-
      lifted %>%
      dplyr::transmute(chr = seqnames, pos = start, pos_old) %>%
      dplyr::inner_join(variants %>% dplyr::rename(pos_old = pos))

    if (keep_old_pos == FALSE) {
      lifted <- lifted %>% dplyr::select(-pos_old)
    }

  }

  # return
  lifted

}