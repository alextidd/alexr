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
    reshape2::dcast(mut_id + n_cells_w_mut ~ cell_id, value.var = "alt_vaf")
  rownames(heatmap_data) <- heatmap_data$mut_id
  heatmap_matrix <-
    heatmap_data %>%
    dplyr::select(-mut_id, -n_cells_w_mut) %>%
    as.matrix()
  n_total_muts <- dplyr::n_distinct(p_dat2$mut_id)
  n_total_cells <- dplyr::n_distinct(p_dat2$cell_id)

  # calculate the number of mutations per column (id) + other annots
  annotations_col <-
    p_dat2 %>%
    dplyr::count(cell_id, !!!syms(annotations)) %>%
    tibble::column_to_rownames("cell_id")

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
