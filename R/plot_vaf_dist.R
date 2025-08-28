plot_vaf_dist <- function(p_dat, p_title = "") {
  n_unique_muts <-
    p_dat %>%
    dplyr::distinct(chr, pos, ref, alt) %>%
    nrow() %>%
    format(big.mark = ",", scientific = FALSE)
  n_total_muts <-
    p_dat %>%
    nrow() %>%
    format(big.mark = ",", scientific = FALSE)

  p_dat_binned <-
    p_dat %>%
    dplyr::mutate(
      alt_vaf = ifelse(total_depth == 0, 0, alt_vaf),
      total_depth_bin = cut(
        total_depth, breaks = c(-Inf, 0, 10, 20, 50, 100, Inf),
        labels = c("0", "<10", "10–20", "20–50", "50–100", ">100")),
      alt_vaf_bin = cut(alt_vaf, breaks = seq(0, 1, length.out = 30 + 1),
                        include.lowest = TRUE),
      alt_vaf_min = stringr::str_extract(alt_vaf_bin, "[\\d\\.]+") %>% as.numeric(),
      alt_vaf_high = as.numeric(stringr::str_extract(alt_vaf_bin, "(?<=,)\\s*[\\d\\.]+")),
      alt_vaf_mid = (alt_vaf_min + alt_vaf_high) / 2) %>%
    dplyr::count(alt_vaf_bin, alt_vaf_mid, total_depth_bin)
  p <-
    p_dat_binned %>%
    ggplot2::ggplot(ggplot2::aes(x = alt_vaf_mid, y = n, fill = total_depth_bin)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0(p_title, "\nVAF distribution"),
                  subtitle = paste(n_total_muts, "total mutations,",
                                   n_unique_muts, "unique mutations"),
         x = "VAF bin", y = "n mutations") +
    ggplot2::lims(x = c(0, 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  print(p)
}