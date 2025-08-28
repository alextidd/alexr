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
    dplyr::mutate(n_cells = dplyr::n_distinct(cell_id)) %>%
    dplyr::add_count(mut_id, name = "n_cells_w_mut") %>%
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