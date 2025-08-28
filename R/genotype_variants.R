genotype_variants <- function(variants, bam, min_bq, min_mq, mask) {
  
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