#' Annotate variants using dNdScv
#'
#' Annotates a set of variants using the dNdScv package to identify coding mutations
#' and their functional impact. The function handles multiple annotations per mutation
#' and maintains the original variant structure.
#'
#' @param variants A data frame containing variant data with columns:
#'   \describe{
#'     \item{sample_id_col}{Column containing sample/donor identifiers (default: donor_id)}
#'     \item{chr}{Chromosome}
#'     \item{pos}{Position}
#'     \item{ref}{Reference allele}
#'     \item{alt}{Alternative allele}
#'     \item{...}{Additional columns will be preserved in the output}
#'   }
#' @param refcds Reference database for dNdScv annotation. Should be a valid RefCDS
#'   object or path to reference data compatible with dNdScv.
#' @param sample_id_col Character string specifying the column name to use as sampleID
#'   for dNdScv. Default is "donor_id".
#'
#' @return A data frame with the original variant data joined with dNdScv annotations.
#'   When multiple annotations exist for the same mutation, they are collapsed into
#'   comma-separated strings. The output maintains one row per mutation-cell combination.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Prepares variant data in dNdScv input format
#'   \item Runs dNdScv annotation with no limits on mutations per gene or sample
#'   \item Collapses multiple annotations for the same mutation into comma-separated values
#'   \item Joins annotations back to the original variant data
#'   \item Preserves all original columns from the input variants data frame
#' }
#'
#' The dNdScv parameters used are:
#' \itemize{
#'   \item \code{max_muts_per_gene_per_sample = Inf} (no gene-level limits)
#'   \item \code{max_coding_muts_per_sample = Inf} (no sample-level limits)
#'   \item \code{outp = 1} (return annotated mutations)
#' }
#'
#' @examples
#' \dontrun{
#' # Example variant data with default donor_id column
#' variants <- data.frame(
#'   donor_id = c("sample1", "sample1", "sample2"),
#'   chr = c("17", "12", "17"),
#'   pos = c(7577120, 25398281, 7577121),
#'   ref = c("C", "G", "A"),
#'   alt = c("T", "A", "G"),
#'   cell_id = c("cell1", "cell2", "cell3")
#' )
#'
#' # Annotate variants using default donor_id column
#' annotated <- annotate_variants(variants, "path/to/refcds")
#' 
#' # Annotate variants using custom sample column
#' variants2 <- data.frame(
#'   sample_name = c("sample1", "sample1", "sample2"),
#'   chr = c("17", "12", "17"),
#'   pos = c(7577120, 25398281, 7577121),
#'   ref = c("C", "G", "A"),
#'   alt = c("T", "A", "G")
#' )
#' annotated2 <- annotate_variants(variants2, "path/to/refcds", sample_id_col = "sample_name")
#' }
#'
#' @importFrom dndscv dndscv
#' @importFrom dplyr select distinct rename group_by summarise across right_join everything
#' @importFrom rlang sym :=
#' @export
annotate_variants <- function(variants, refcds, sample_id_col = "donor_id") {

  # annotate mutations with dndscv
  dndscv_in <-
    variants %>%
    dplyr::select(sampleID = !!rlang::sym(sample_id_col), chr, pos, ref, mut = alt) %>%
    dplyr::distinct()
  dndscv_out <-
    dndscv::dndscv(dndscv_in, max_muts_per_gene_per_sample = Inf,
                   max_coding_muts_per_sample = Inf, outp = 1,
                   refdb = refcds)

  # collapse multiple annotations of the same mutation
  # (we want to maintain 1 mutation-x-cell per row)
  annots <-
    dndscv_out$annotmuts %>%
    dplyr::rename(!!rlang::sym(sample_id_col) := sampleID) %>%
    dplyr::rename(alt = mut) %>%
    dplyr::group_by(!!rlang::sym(sample_id_col), chr, pos, ref, alt) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(.x, collapse = ",")),
                     .groups = "drop")

  # join variants with dndscv output, collapse multiple annotations
  dplyr::right_join(annots, variants)

}