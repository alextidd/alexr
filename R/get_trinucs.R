#' Get trinucleotide contexts for SNV mutations
#'
#' This function extracts trinucleotide contexts for single nucleotide variants (SNVs)
#' and converts mutations to pyrimidine base representation for standardized analysis.
#' Only variants with single base reference and alternative alleles are processed.
#'
#' @param variants A data frame containing variant information with columns: chr, pos, ref, alt.
#'   Only SNVs (single nucleotide variants) will be processed.
#' @param fasta A Rsamtools FaFile object or path to a FASTA file containing the reference genome.
#'
#' @return A data frame with the original variant information plus additional columns:
#'   \describe{
#'     \item{trinuc_ref}{The trinucleotide context (3 bases) centered on the variant position}
#'     \item{sub}{The substitution in format "REF>ALT"}
#'     \item{sub_py}{The substitution converted to pyrimidine base representation}
#'     \item{trinuc_ref_py}{The trinucleotide context converted to pyrimidine base representation}
#'   }
#'   
#' @details 
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters variants to include only SNVs (single base ref and alt)
#'   \item Extracts trinucleotide context (position-1 to position+1) from the reference genome
#'   \item Converts purine bases (A, G) to their pyrimidine complements (T, C) for standardized representation
#' }
#'
#' @examples
#' \dontrun{
#' library(Rsamtools)
#' variants <- data.frame(
#'   chr = "chr1", 
#'   pos = c(1000, 2000), 
#'   ref = c("A", "G"), 
#'   alt = c("T", "C")
#' )
#' fasta <- FaFile("reference.fa")
#' result <- get_trinucs(variants, fasta)
#' }
#'
#' @importFrom dplyr filter rowwise mutate ungroup case_when
#' @importFrom Rsamtools scanFa
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
get_trinucs <- function(variants, fasta) {

  variants %>%
    # get snvs only
    dplyr::filter(nchar(ref) == 1, nchar(alt) == 1) %>%
    # get the trinucleotide context
    dplyr::rowwise() %>%
    dplyr::mutate(
      trinuc_ref = Rsamtools::scanFa(fasta,
                                     GRanges(chr, IRanges(pos - 1,
                                                          pos + 1))) %>%
                    as.vector()) %>%
    dplyr::ungroup() %>%
    # annotate the mutation from the pyrimidine base
    dplyr::mutate(
      sub = paste(ref, alt, sep = ">"),
      sub_py = dplyr::case_when(
        ref %in% c("A", "G") ~ chartr("TCGA", "AGCT", sub),
        TRUE ~ sub),
      trinuc_ref_py = dplyr::case_when(
        ref %in% c("A", "G") ~ chartr("TCGA", "AGCT", trinuc_ref),
        TRUE ~ trinuc_ref))

}
