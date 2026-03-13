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
#' @param pileup Logical; if TRUE, returns full pileup counts for each base and indel at each site. If FALSE, returns summary statistics for ref and alt.
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
#'
#' @return If pileup = FALSE: a tibble with columns chr, pos, ref, alt, total_depth, alt_depth, ref_depth, alt_vaf (variant allele frequency, alt_depth / total_depth). If pileup = TRUE: a tibble with columns chr, pos, A, T, C, G, N, ins, del, total_depth (base and indel counts per site).
#'   Complex variants are automatically filtered out and not genotyped.
#'
#' @examples
#' \dontrun{
#' variants <- data.frame(chr = "1", pos = 1000, ref = "A", alt = "T")
#' genotype_variants(variants, "sample.bam", min_bq = 30, min_mq = 30, mask = 3844)
#' genotype_variants(variants, "sample.bam", min_bq = 30, min_mq = 30, pileup = TRUE)
#' }
#'
#' @export
genotype_variants <- function(variants, bam, min_bq, min_mq, mask = 0, pileup = FALSE) {

  if (pileup == FALSE) {
    # get and type variants
    vars <- alexr::type_variants(variants)

    # check for variants that are not snv / ins / del
    if ("complex" %in% vars$type) {
      message("Complex variants are not supported!")
      message(paste(vars %>% dplyr::filter(type == "complex") %>% nrow(),
                    "complex variant(s) were found and will be removed."))
      vars <- vars %>% dplyr::filter(type != "complex")
    }
  } else {
    vars <- variants
  }

  # genotype all sites
  geno <-
    vars %>%
    dplyr::select(chr, pos, dplyr::any_of(if (pileup) character(0) else c("ref", "alt", "type"))) %>%
    dplyr::distinct() %>%
    purrr::pmap(function(chr, pos, type = "", ref = NULL, alt = NULL) {
      # testing: # vars %>% dplyr::distinct(chr, pos, ref, alt, type) %>% as.list() %>% list2env(envir = globalenv()); min_bq <- 30; min_mq <- 30; mask <- 3844

      paste(chr, pos, ref, alt, type, "\n") %>% cat()

      # look ahead if deletion
      if (type == "del") {
        pos_i <- pos + 1
      } else {
        pos_i <- pos
      }

      # query bam
      calls <- deepSNV::bam2R(bam, chr, pos_i, pos_i, mask = mask,
                              q = min_bq, mq = min_mq)

      # calculate total depth
      total_depth <-
        sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t", "-", "_",
                      "N", "n")],
            na.rm = TRUE)
      
      # return full pileup
      if (pileup == TRUE) {

        tibble::tibble(
          chr = chr, pos = pos,
          A   = calls[1, "A"]   + calls[1, "a"],
          T   = calls[1, "T"]   + calls[1, "t"],
          C   = calls[1, "C"]   + calls[1, "c"],
          G   = calls[1, "G"]   + calls[1, "g"],
          N   = calls[1, "N"]   + calls[1, "n"],
          ins = calls[1, "INS"] + calls[1, "ins"],
          del = calls[1, "-"]   + calls[1, "_"],
          total_depth = total_depth)

      } else {

        # check type
        if (!(type %in% c("snv", "dnv", "mnv", "ins", "del"))) {
          stop(paste0("Variant type '", type, "' not recognized!"))
        }

        # calculate ref depth
        if (type == "ins") {
          # ref = no insertion; everything that isn't an insertion read
          ref_depth <- total_depth - alt_depth
        } else {
          # snv / del: ref = reads matching the first base of ref allele
          # (for del, bam2R is already queried at pos+1, so ref_i is correct)
          ref_i <- unlist(strsplit(ref, ""))[1]
          ref_depth <- calls[1, ref_i] + calls[1, tolower(ref_i)]
        }

        # calculate alt depth
        if (type == "del") {
          alt_depth <- sum(calls[, c("-", "_")])
        } else if (type == "ins") {
          alt_depth <- sum(calls[, c("INS", "ins")])
        } else {
          alt_i <- unlist(strsplit(alt, ""))[1]
          alt_depth <- calls[1, alt_i] + calls[1, tolower(alt_i)]
        }

        # return
        tibble::tibble(chr = chr, pos = pos, ref = ref, alt = alt,
                       total_depth = total_depth,
                       alt_depth = alt_depth, ref_depth = ref_depth) %>%
          dplyr::mutate(alt_vaf = alt_depth / total_depth)

      }

    }) %>%
    dplyr::bind_rows()
  
  # return
  return(geno)

}
