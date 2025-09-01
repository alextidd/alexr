#' Type mutations based on reference and alternative alleles
#'
#' Classifies mutations into different types (SNV, insertion, deletion, DNV, MNV)
#' based on the length of reference and alternative alleles. Optionally splits
#' complex mutations (DNVs and MNVs) into individual single nucleotide changes.
#'
#' @param variants A data frame containing variants with columns:
#'   \describe{
#'     \item{ref}{Reference allele sequence}
#'     \item{alt}{Alternative allele sequence}
#'     \item{chr}{Chromosome (required if split_complex = TRUE)}
#'     \item{pos}{Position (required if split_complex = TRUE)}
#'     \item{...}{Additional columns will be preserved in the output}
#'   }
#' @param split_complex Logical indicating whether to split DNVs and MNVs into
#'   individual nucleotide changes. Default is FALSE. When TRUE, multi-nucleotide
#'   variants are expanded to multiple rows, one for each position, but still 
#'   retain the "dnv" or "mnv" label.
#'
#' @return A data frame with the original data plus a new \code{type} column
#'   containing mutation classifications:
#'   \describe{
#'     \item{snv}{Single nucleotide variant (ref and alt both length 1)}
#'     \item{ins}{Insertion (ref length 1, alt length > 1)}
#'     \item{del}{Deletion (ref length > 1, alt length 1)}
#'     \item{dnv}{Dinucleotide variant (ref and alt both length 2)}
#'     \item{mnv}{Multi-nucleotide variant (ref and alt both length > 2)}
#'   }
#'   When \code{split_complex = TRUE}, DNVs and MNVs are converted to multiple
#'   SNV rows with updated positions and alleles.
#'
#' @details 
#' The function classifies mutations based on allele lengths:
#' \itemize{
#'   \item \strong{SNV}: Single nucleotide changes (1 base -> 1 base)
#'   \item \strong{Insertion}: Reference is 1 base, alternative is longer
#'   \item \strong{Deletion}: Reference is longer, alternative is 1 base
#'   \item \strong{DNV}: Both reference and alternative are 2 bases
#'   \item \strong{MNV}: Both reference and alternative are >2 bases
#' }
#'
#' When \code{split_complex = TRUE}, the function:
#' \itemize{
#'   \item Identifies DNVs and MNVs in the input data
#'   \item Splits each complex variant into constituent single nucleotide changes
#'   \item Updates positions incrementally for each nucleotide
#'   \item Creates separate rows for each position with corresponding ref/alt bases
#'   \item Preserves all other columns from the original data
#' }
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' 
#' # Example mutation data
#' mutations <- data.frame(
#'   chr = c("1", "2", "3", "4", "5"),
#'   pos = c(1000, 2000, 3000, 4000, 5000),
#'   ref = c("A", "T", "ATG", "CG", "ATCG"),
#'   alt = c("G", "TCC", "A", "GA", "GCTA"),
#'   sample = c("S1", "S2", "S3", "S4", "S5")
#' )
#' 
#' # Basic typing without splitting complex variants
#' typed <- type_mutations(mutations)
#' 
#' # Typing with complex variant splitting
#' typed_split <- type_mutations(mutations, split_complex = TRUE)
#' }
#'
#' @importFrom dplyr mutate case_when filter group_by across reframe select bind_rows
#' @export
type_mutations <- function(variants, split_complex = FALSE) {
  typed_variants <-
    variants %>%
    dplyr::mutate(type = dplyr::case_when(
      nchar(ref) == 1 & nchar(alt) == 1 ~ "snv",
      nchar(ref) == 1 & nchar(alt) > 1 ~ "ins",
      nchar(ref) > 1 & nchar(alt) == 1 ~ "del",
      nchar(ref) == 2 & nchar(alt) == 2 ~ "dnv",
      nchar(ref) > 2 & nchar(alt) > 2 ~ "mnv"
    ))

  # expand dnv/mnv mutations to all positions
  if (split_complex == TRUE) {
    if (any(typed_variants$type %in% c("dnv", "mnv"))) {
      typed_variants_dnv_mnv <-
        typed_variants %>%
        dplyr::filter(type %in% c("dnv", "mnv")) %>%
        dplyr::mutate(mut_id_x = paste(chr, pos, ref, alt, type)) %>%
        dplyr::group_by(dplyr::across(-c(pos, alt, ref))) %>%
        dplyr::reframe(pos = pos:(pos + nchar(alt) - 1),
                       ref = strsplit(ref, split = "") %>% unlist(),
                       alt = strsplit(alt, split = "") %>% unlist()) %>%
        dplyr::select(-mut_id_x)
      typed_variants <-
        typed_variants %>%
        dplyr::filter(!(type %in% c("dnv", "mnv"))) %>%
        dplyr::bind_rows(typed_variants_dnv_mnv)
    }
  }

  # return
  typed_variants
}