#' Plot B-allele frequency (BAF) across chromosomes
#'
#' Creates a genome-wide plot of B-allele frequencies showing both VAF and BAF
#' (1-VAF) values across all chromosomes. The plot includes alternating chromosome
#' backgrounds, gene annotations, and a reference line at 0.5.
#'
#' @param p_dat A data frame containing variant data with columns:
#'   \describe{
#'     \item{chr}{Chromosome (will be converted to factor with levels 1-22, X, Y)}
#'     \item{pos}{Position on chromosome}
#'     \item{mut_vaf}{Mutation variant allele frequency}
#'     \item{total_depth}{Total sequencing depth (filtered to >20)}
#'   }
#' @param genes A character vector of gene names for annotation, or NULL.
#'   If provided, gene coordinates will be looked up using the refcds file.
#'   Default is NULL.
#' @param refcds Path to RefCDS reference file for gene coordinate lookup.
#'   Required if genes are specified. Default is NULL.
#' @param p_alpha Transparency level for points (0-1). Default is 0.05.
#' @param p_size Size of points. Default is 0.01.
#'
#' @return A ggplot object showing BAF values across chromosomes.
#'
#' @details 
#' The function creates a BAF plot with the following features:
#' \itemize{
#'   \item Filters variants to those with total_depth > 20
#'   \item Calculates BAF as 1 - VAF (B-allele frequency)
#'   \item Displays both VAF and BAF values as points
#'   \item Uses alternating gray/white backgrounds for chromosomes
#'   \item Adds a red reference line at 0.5 (expected heterozygous frequency)
#'   \item Facets by chromosome with free scaling
#'   \item Optional gene annotation as vertical lines
#'   \item Removes x-axis labels for cleaner display
#' }
#'
#' @note This function has dependencies on undefined variables (opts, p_source)
#'   that may need to be defined in the calling environment or passed as parameters.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#' library(ggh4x)
#' 
#' # Example variant data
#' variants <- data.frame(
#'   chr = c("1", "2", "1", "X"),
#'   pos = c(1000000, 2000000, 1500000, 500000),
#'   mut_vaf = c(0.3, 0.7, 0.45, 0.6),
#'   total_depth = c(50, 80, 60, 40)
#' )
#' 
#' # Basic BAF plot
#' plot_baf(variants)
#' 
#' # With gene annotations
#' plot_baf(variants, genes = c("TP53", "KRAS"), refcds = "path/to/refcds.RData")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_rect geom_hline geom_point geom_vline
#' @importFrom ggplot2 scale_x_continuous scale_fill_manual theme_classic theme
#' @importFrom ggplot2 element_rect element_blank labs unit
#' @importFrom dplyr filter mutate group_by cur_group_id bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggh4x facet_grid2
#' @importFrom purrr set_names map map_lgl
#' @importFrom tibble tibble
#' @export
plot_baf <- function(p_dat, genes = NULL, refcds = NULL, p_alpha = 0.05, p_size = 0.01) {

  # libraries
  library(ggplot2)

  # get gene coordinates
  if (!is.null(genes)) {

    # read refcds
    load(refcds)

    # get gene coordinates
    p_genes <-
      genes %>%
      purrr::set_names() %>%
      purrr::map(function(g) {
        i <- RefCDS[[which(purrr::map_lgl(RefCDS, ~ .x$gene_name == g))]]
        tibble::tibble(chr = i$chr,
                      pos = (min(i$intervals_cds) + max(i$intervals_cds)) / 2)
      }) %>%
      dplyr::bind_rows(.id = "gene")

  } else {

    p_genes <- NULL

  }

  # prep data
  p_dat2 <-
    p_dat %>%
    dplyr::filter(total_depth > 20) %>%
    dplyr::mutate(mut_baf = 1 - mut_vaf,
                  chr = factor(sub("^chr", "", chr),
                               levels = c(as.character(1:22), "X", "Y"))) %>%
    # get alternating colours and chr bounds
    dplyr::group_by(chr) %>%
    dplyr::mutate(chr_alternating = dplyr::cur_group_id() %% 2,
                  min_pos = min(pos), max_pos = max(pos)) %>%
    # plot
    tidyr::pivot_longer(cols = c("mut_vaf", "mut_baf"), names_to = "vaf_type")

  # plot
  p <-
    p_dat2 %>%
    ggplot(aes(x = pos, y = value)) +
    # add alternating coloured chromosomes
    geom_rect(aes(ymin = 0, ymax = 1, xmin = min_pos, xmax = max_pos,
                  fill = as.factor(chr_alternating)),
              show.legend = FALSE) +
    # add line at vaf = 0.5
    geom_hline(yintercept = 0.5, colour = "red") +
    # add baf points
    geom_point(size = p_size, alpha = p_alpha) +
    # add baf bands
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("white", "#e3e3e3")) +
    ggh4x::facet_grid2(. ~ chr, scales = "free_x", space = "free_x") +
    theme_classic() +
    theme(panel.spacing = unit(0, "lines"),
          panel.border = element_rect(color = "grey", fill = NA,
                                      linewidth = 0),
          strip.background = element_rect(color = "grey", fill = NA,
                                          linewidth = 0, linetype = "solid"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = paste(opts$id, "-", p_source))

  # add gene annotations
  if (!is.null(p_genes)) {
    p <-
      p +
      geom_vline(
        data = p_genes %>% dplyr::mutate(chr = factor(chr, levels = levels(p_dat2$chr))),
        aes(xintercept = pos, colour = gene))
  }

  # return
  p

}