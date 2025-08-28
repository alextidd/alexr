#' Genotype variants in a BAM file using deepSNV::bam2R.
#'
#' This function genotypes a set of variants by querying a BAM file and
#' calculating variant allele frequencies (VAF) for each variant. It supports
#' SNVs, insertions, and deletions, but filters out complex variants.
#'
#' @param variants A data frame containing variant information with columns: chr, pos, ref, alt.
#' @param bam Path to the BAM file to query for genotyping.
#' @param min_bq Minimum base quality score for reads.
#' @param min_mq Minimum mapping quality score for reads.
#' @param mask Genomic regions to mask during analysis.
#'
#' @return A tibble with columns: chr, pos, ref, alt, total_depth, alt_depth, alt_vaf.
#'   The alt_vaf column contains the variant allele frequency calculated as alt_depth / total_depth.
#' @examples
#' \dontrun{
#' variants <- data.frame(chr = "1", pos = 1000, ref = "A", alt = "T")
#' genotype_variants(variants, "sample.bam", min_bq = 30, min_mq = 30, mask = 3844)
#' }
#'
#' @export
genotype_variants <- function(variants, bam, min_bq, min_mq, mask = 0) {
  
  # get and type variants
  vars <-
    variants %>%
    dplyr::mutate(
      type = dplyr::case_when(nchar(ref) == 1 & nchar(alt) == 1 ~ "snv",
                              nchar(ref) == 1 & nchar(alt) > 1 ~ "ins",
                              nchar(ref) > 1 & nchar(alt) == 1 ~ "del",
                              TRUE ~ "complex"))

  # check for variants that are not snv / ins / del
  if ("complex" %in% vars$type) {
    message("Complex variants are not supported!")
    message(paste(vars %>% dplyr::filter(type == "complex") %>% nrow(),
                  "complex variant(s) were found and will be removed."))
    vars <- vars %>% dplyr::filter(type != "complex")
  }

  # genotype all sites
  geno <-
    vars %>%
    dplyr::distinct(chr, pos, ref, alt, type) %>%
    purrr::pmap(function(chr, pos, ref, alt, type) {

      paste(chr, pos, ref, alt, type, "\n") %>% cat()

      # check type
      if (!(type %in% c("snv", "dnv", "mnv", "ins", "del"))) {
        stop(paste0("Variant type '", type, "' not recognized!"))
      }

      # look ahead if deletion
      if (type == "del") {
        pos_i <- pos + 1
      } else {
        pos_i <- pos
      }

      # query bam
      calls <- deepSNV::bam2R(bam, chr, pos_i, pos_i, mask = opts$mask,
                              q = opts$min_bq, mq = opts$min_mq)

      # calculate total depth
      total_depth <-
        sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t", "-", "_")],
            na.rm = TRUE)

      # calculate alt depth
      if (type == "del") {
        alt_depth <- sum(calls[, c("-", "_")])
      } else if(type == "ins") {
        alt_depth <- sum(calls[, c("INS", "ins")])
      } else {
        alt_i <- unlist(strsplit(alt, ""))[1]
        alt_depth <- calls[1, alt_i] + calls[1, tolower(alt_i)]
      }

      # return
      tibble::tibble(chr = chr, pos = pos, ref = ref, alt = alt,
                     total_depth = total_depth, alt_depth = alt_depth) %>%
        dplyr::mutate(alt_vaf = alt_depth / total_depth)

    }) %>%
    dplyr::bind_rows()

}