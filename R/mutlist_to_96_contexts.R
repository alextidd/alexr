#' Convert mutation list to 96-context trinucleotide mutation matrix
#'
#' Converts a list of mutations into a 96-trinucleotide context mutation matrix
#' for mutational signature analysis. The function extracts trinucleotide contexts
#' from a reference genome and converts all mutations to pyrimidine base representation.
#'
#' @param muts A data frame containing mutation data with columns:
#'   \describe{
#'     \item{sample_id}{Sample identifier}
#'     \item{chr}{Chromosome (1-22, X, Y)}
#'     \item{pos}{Position on chromosome}
#'     \item{ref}{Reference allele (A, C, G, T)}
#'     \item{alt}{Alternative allele (A, C, G, T)}
#'   }
#' @param fasta A Rsamtools FaFile object or path to indexed FASTA file containing
#'   the reference genome sequence for trinucleotide context extraction.
#'
#' @return A numeric matrix with 96 columns and one row per sample. Columns represent
#'   the 96 possible trinucleotide contexts (6 substitution types × 16 trinucleotide
#'   contexts). Row names are sample IDs and column names follow the format
#'   "substitution,flanking-bases" (e.g., "C>A,A-A").
#'
#' @details 
#' The function performs the following steps for each sample:
#' \itemize{
#'   \item Filters mutations to standard chromosomes (1-22, X, Y) and valid bases (A,C,G,T)
#'   \item Extracts trinucleotide context (position ± 1) from the reference genome
#'   \item Converts all mutations to pyrimidine base representation (C>X or T>X)
#'   \item For purine reference bases (A,G), converts to complement strand
#'   \item Counts mutations for each of the 96 possible trinucleotide contexts
#'   \item Returns a matrix suitable for mutational signature analysis
#' }
#'
#' The 96 trinucleotide contexts consist of:
#' \itemize{
#'   \item 6 substitution types: C>A, C>G, C>T, T>A, T>C, T>G
#'   \item 16 trinucleotide contexts: all combinations of flanking bases (A,C,G,T)
#' }
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' library(Rsamtools)
#' 
#' # Example mutation data
#' mutations <- data.frame(
#'   sample_id = c("sample1", "sample1", "sample2", "sample2"),
#'   chr = c("1", "2", "1", "3"),
#'   pos = c(1000000, 2000000, 1500000, 3000000),
#'   ref = c("C", "G", "A", "T"),
#'   alt = c("T", "A", "G", "C")
#' )
#' 
#' # Generate 96-context matrix
#' mut_matrix <- muts_to_96_contexts(mutations, "reference.fa")
#' }
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools scanFa
#' @export
muts_to_96_contexts <- function(muts, fasta) {

  # libraries
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")

  # list of samples
  samples <- unique(muts$sample_id)

  # initiate matrix
  trinuc_mut_mat <- matrix(0, ncol = 96, nrow = length(samples))

  for (n in seq_along(samples)) {
    s <- samples[n]
    mutations <-
      as.data.frame(muts[muts$sample_id == s, c("chr", "pos", "ref", "alt")])
    colnames(mutations) <- c("chr", "pos", "ref", "mut")
    mutations$pos <- as.numeric(mutations$pos)
    mutations <- mutations[(mutations$ref %in% c("A", "C", "G", "T")) &
                             (mutations$mut %in% c("A", "C", "G", "T")) &
                             mutations$chr %in% paste0("chr", c(1:22, "X", "Y")), ]

    mutations$trinuc_ref <-
      as.vector(scanFa(fasta, GRanges(mutations$chr,
                                      IRanges(mutations$pos - 1,
                                              mutations$pos + 1))))

    ntcomp <- c(T = "A", G = "C", C = "G", A = "T")
    mutations$sub <- paste(mutations$ref, mutations$mut, sep = ">")
    mutations$trinuc_ref_py <- mutations$trinuc_ref

    for (j in seq_along(mutations)) {
      if (mutations$ref[j] %in% c("A", "G")) {
        mutations$sub[j] <-
          paste(ntcomp[mutations$ref[j]], ntcomp[mutations$mut[j]], sep = ">")
        mutations$trinuc_ref_py[j] <-
          paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j], split = "")[[1]])],
                collapse = "")
      }
    }

    freqs <-
      table(paste0(mutations$sub, ",", substr(mutations$trinuc_ref_py, 1, 1),
                   "-", substr(mutations$trinuc_ref_py, 3, 3)))
    sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    ctx_vec <- paste(rep(c("A", "C", "G", "T"), each = 4),
                     rep(c("A", "C", "G", "T"), times = 4),
                     sep = "-")
    full_vec <- paste(rep(sub_vec, each = 16),
                      rep(ctx_vec, times = 6),
                      sep = ",")

    freqs_full <- freqs[full_vec]
    freqs_full[is.na(freqs_full)] <- 0
    names(freqs_full) <- full_vec
    trinuc_mut_mat[n, ] <- freqs_full
    print(s)
  }

  colnames(trinuc_mut_mat) <- full_vec
  rownames(trinuc_mut_mat) <- samples
  return(trinuc_mut_mat)
}
