#' Plot trinucleotide mutation signature
#'
#' Creates a bar plot showing the distribution of mutations across different
#' trinucleotide contexts, grouped by substitution type. The plot displays
#' mutation counts for each trinucleotide context in pyrimidine base representation.
#'
#' @param p_dat A data frame containing mutation data with columns:
#'   \describe{
#'     \item{sub_py}{Substitution type in pyrimidine representation (e.g., "C>A", "T>G")}
#'     \item{trinuc_ref_py}{Trinucleotide context in pyrimidine representation}
#'   }
#' @param p_title Character string for the plot title. Default is an empty string.
#'
#' @return A ggplot object showing mutation counts by trinucleotide context,
#'   faceted by substitution type with color-coded bars.
#'
#' @details 
#' The function creates a mutation signature plot with the following features:
#' \itemize{
#'   \item Six substitution types: C>A, C>G, C>T, T>A, T>C, T>G
#'   \item Each substitution type is displayed in a separate facet
#'   \item X-axis shows all possible trinucleotide contexts for each substitution
#'   \item Y-axis shows the number of mutations
#'   \item Colors follow standard mutation signature conventions
#'   \item Missing trinucleotide contexts are filled with zero counts
#' }
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#' 
#' # Example mutation data
#' mutations <- data.frame(
#'   sub_py = c("C>A", "C>T", "T>G"),
#'   trinuc_ref_py = c("ACA", "TCG", "ATG")
#' )
#' 
#' plot_mut_trinucs(mutations, "Mutation Signature")
#' }
#'
#' @importFrom dplyr mutate n_distinct add_count count right_join filter
#' @importFrom ggplot2 ggplot aes geom_col facet_grid guides guide_axis theme_minimal
#' @importFrom ggplot2 theme element_text scale_fill_manual scale_y_continuous labs
#' @export
plot_mut_trinucs <- function(p_dat, p_title = "") {
  sub_colours <-
    c("C>A" = "dodgerblue", "C>G" = "black", "C>T" = "red",
      "T>A" = "grey70", "T>C" = "olivedrab3", "T>G" = "plum2")
  nucs <- c("T", "C", "G", "A")
  py_nucs <- c("C", "T")
  all_subs <-
    expand.grid(do.call(paste0, expand.grid(nucs, py_nucs, nucs)), nucs) %>%
    dplyr::transmute(trinuc_ref_py = Var1, ref = substr(Var1, 2, 2), alt = Var2,
                     sub_py = paste0(ref, ">", alt)) %>%
    dplyr::filter(alt != ref)

  p_dat %>%
    dplyr::count(sub_py, trinuc_ref_py) %>%
    dplyr::right_join(all_subs) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    ggplot(aes(x = trinuc_ref_py, y = n, fill = sub_py)) +
    geom_col() +
    facet_grid(~ sub_py, scales = "free_x") +
    guides(x = guide_axis(angle = 90)) +
    theme_minimal() +
    theme(axis.text.x = element_text(family = "mono"),
          legend.position = "none") +
    scale_fill_manual(values = sub_colours) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = p_title,
         subtitle = paste(format(nrow(p_dat), big.mark = ","), "mutations"),
         x = "", y = "number of mutations")
}