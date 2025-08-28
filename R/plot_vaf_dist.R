#' Plot variant allele frequency (VAF) distribution
#'
#' Creates a histogram showing the distribution of variant allele frequencies (VAF)
#' stratified by total sequencing depth. The plot displays VAF values binned into
#' 30 equal intervals and colored by depth categories.
#'
#' @param p_dat A data frame containing variant data with columns:
#'   \describe{
#'     \item{chr}{Chromosome}
#'     \item{pos}{Position}
#'     \item{ref}{Reference allele}
#'     \item{alt}{Alternative allele}
#'     \item{alt_vaf}{Variant allele frequency (0-1)}
#'     \item{total_depth}{Total sequencing depth at the position}
#'   }
#' @param p_title Character string for the plot title prefix. Default is an empty string.
#'   The final title will be "p_title\nVAF distribution".
#'
#' @return Invisibly returns NULL. The function prints the plot directly using print().
#'
#' @details 
#' The function creates a VAF distribution plot with the following features:
#' \itemize{
#'   \item VAF values are binned into 30 equal intervals from 0 to 1
#'   \item Total depth is categorized into 6 bins: 0, <10, 10-20, 20-50, 50-100, >100
#'   \item Colors represent different depth categories using viridis color scale
#'   \item Subtitle shows total number of mutations and unique mutations
#'   \item X-axis ranges from 0 to 1 (VAF)
#'   \item Y-axis shows mutation counts with no expansion at zero
#' }
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#' 
#' # Example variant data
#' variants <- data.frame(
#'   chr = c("chr1", "chr2", "chr1"),
#'   pos = c(1000, 2000, 3000),
#'   ref = c("A", "G", "C"),
#'   alt = c("T", "C", "A"),
#'   alt_vaf = c(0.25, 0.45, 0.15),
#'   total_depth = c(50, 80, 25)
#' )
#' 
#' plot_vaf_dist(variants, "Sample 1")
#' }
#'
#' @importFrom dplyr distinct mutate count
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_viridis_d theme_classic
#' @importFrom ggplot2 labs lims scale_y_continuous
#' @importFrom stringr str_extract
#' @export
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