#' Plot variant allele frequency (VAF) heatmap
#'
#' Creates a clustered heatmap showing variant allele frequencies across identifiers
#' and mutations. The heatmap displays VAF values with hierarchical clustering
#' of both rows (mutations) and columns (identifiers), with optional annotations.
#'
#' @param p_dat A data frame containing variant data with columns:
#'   \describe{
#'     \item{gene}{Gene name}
#'     \item{chr}{Chromosome}
#'     \item{pos}{Position}
#'     \item{ref}{Reference allele}
#'     \item{alt}{Alternative allele}
#'     \item{alt_vaf}{Variant allele frequency (0-1)}
#'     \item{id}{Identifier}
#'     \item{...}{Additional columns for annotations (specified in annotations parameter)}
#'   }
#' @param p_title Character string for the plot title prefix. Default is an empty string.
#'   The final title will be "p_title VAF heatmap n IDs, m mutations".
#' @param annotations Character vector of column names in p_dat to use for column annotations.
#'   Default is an empty vector (no annotations).
#' @param show_rownames Logical indicating whether to show row names (mutation IDs) on the heatmap.
#'   Default is TRUE.
#'
#' @return A pheatmap object. The function displays the heatmap directly.
#'
#' @details 
#' The function creates a VAF heatmap with the following features:
#' \itemize{
#'   \item Filters to mutations with VAF > 0
#'   \item Creates mutation IDs in format: gene_chr_pos_ref_alt
#'   \item Hierarchical clustering of both rows (mutations) and columns (IDs)
#'   \item Row annotations show the number of IDs with each mutation
#'   \item Column annotations can be specified for additional ID metadata
#'   \item Uses reversed magma color palette (viridis package)
#'   \item NA values are replaced with 0
#'   \item Title includes total ID and mutation counts
#' }
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(pheatmap)
#' 
#' # Example variant data
#' variants <- data.frame(
#'   gene = c("TP53", "KRAS", "TP53"),
#'   chr = c("17", "12", "17"),
#'   pos = c(7577120, 25398281, 7577121),
#'   ref = c("C", "G", "A"),
#'   alt = c("T", "A", "G"),
#'   alt_vaf = c(0.45, 0.32, 0.28),
#'   id = c("cell1", "cell2", "cell1"),
#'   sample_type = c("tumor", "tumor", "normal")
#' )
#' 
#' plot_vaf_heatmap(variants, "Sample 1", annotations = "sample_type")
#' }
#'
#' @importFrom dplyr ungroup filter mutate add_count select n_distinct count
#' @importFrom tibble column_to_rownames
#' @importFrom reshape2 dcast
#' @importFrom pheatmap pheatmap
#' @importFrom viridis magma
#' @importFrom rlang syms
#' @export
plot_vaf_heatmap <- function(p_dat, p_title = "", annotations = c(), show_rownames = TRUE) {

  # prepare data
  p_dat2 <-
    p_dat %>%
    dplyr::ungroup() %>%
    # only keep mutations with VAF > 0
    dplyr::filter(alt_vaf > 0) %>%
    dplyr::mutate(mut_id = paste(gene, chr, pos, ref, alt, sep = "_")) %>%
    # count number of cells with each mutation
    dplyr::add_count(mut_id, name = "n_cells_w_mut")

  # reshape data for heatmap
  heatmap_data <-
    p_dat2 %>%
    reshape2::dcast(mut_id + n_cells_w_mut ~ id, value.var = "alt_vaf")
  rownames(heatmap_data) <- heatmap_data$mut_id
  heatmap_matrix <-
    heatmap_data %>%
    dplyr::select(-mut_id, -n_cells_w_mut) %>%
    as.matrix()
  n_total_muts <- dplyr::n_distinct(p_dat2$mut_id)
  n_total_cells <- dplyr::n_distinct(p_dat2$id)

  # calculate the number of mutations per column (id) + other annots
  annotations_col <-
    p_dat2 %>%
    dplyr::count(id, !!!syms(annotations)) %>%
    tibble::column_to_rownames("id")

  # replace NAs with 0
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # prepare annotations for columns
  annotation_row <- data.frame(n_cells_w_mut = heatmap_data$n_cells_w_mut)
  rownames(annotation_row) <- heatmap_data$mut_id

  # plot
  pheatmap::pheatmap(
    mat = heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    annotation_row = annotation_row,
    annotation_col = annotations_col,
    color = rev(viridis::magma(100)),
    main = paste0(p_title, "\nVAF heatmap\n",
                  n_total_cells, " cells, ", n_total_muts, " mutations"))

}
