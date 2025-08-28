get_mut_trinucs <- function(dat, fasta) {
  dat %>%
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