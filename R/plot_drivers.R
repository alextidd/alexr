plot_drivers <-
  function(genes2plot,
           annot_muts,
           dndsout,
           scores,
           target_genes,
           newkc,
           RefCDS,
           gr_genes,
           max_genes = Inf,
           only_sig_dnds = TRUE,
           sig_p = 0.01,
           driver_density_sample_ids = NULL,
           gene2dc,
           logscale = TRUE,
           use_indel_sites = TRUE) {

    # libraries
    library(ggplot2)
    library(patchwork)

    # set the theme
    theme_set(theme_classic())

    # impact palettes
    impact_pal <-
      c("indels" = "chocolate3", "splice" = "darkorchid2",
        "nonsense" = "darkorchid4", "missense" = "cadetblue",
        "synonymous" = "grey70")
    impact2_pal <-
      c("missense" = "cadetblue", "nonsense+splice" = "darkorchid2",
        "indels" = "chocolate3")

    # 1) mutations observed
    p_dat_1 <-
      dndsout$sel_cv %>%
      tibble::as_tibble() %>%
      dplyr::select(gene = gene_name, dplyr::starts_with("n_")) %>%
      # filter genes
      dplyr::filter(gene %in% genes2plot) %>%
      tidyr::pivot_longer(cols = -gene, names_to = "impact",
                          values_to = "n") %>%
      # sort genes in decreasing order of number of mutations
      dplyr::group_by(gene) %>%
      dplyr::mutate(n_total = sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(gene = forcats::fct_reorder(gene, n_total,
                                                .desc = TRUE)) %>%
      dplyr::arrange(gene) %>%
      # select first n genes
      dplyr::group_by(gene) %>%
      dplyr::filter(dplyr::cur_group_id() <= max_genes) %>%
      dplyr::ungroup() %>%
      # recode impacts
      dplyr::mutate(impact = dplyr::recode(
                      impact, n_syn = "synonymous",
                      n_mis = "missense", n_spl = "splice", n_non = "nonsense",
                      n_ind = "indels") %>% factor(levels = names(impact_pal)))

    # save gene order based on n mutations
    genes2plot <- levels(p_dat_1$gene)

    # plot
    p_1 <-
      p_dat_1 %>%
      ggplot(aes(x = gene, y = n, fill = impact)) +
      geom_col() +
      scale_fill_manual(values = impact_pal) +
      # second axis with mean n muts per donor
      scale_y_continuous(
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / dplyr::n_distinct(dndsout$annotmuts$sampleID),
        name = "mean n muts per donor")) +
      scale_x_discrete(guide = guide_axis(angle = -90)) +
      labs(y = "total mutations") +
      theme(axis.title.x = element_blank())

    # 2) dN/dS ratios from the dNdScv model
    p_dat_2 <-
      dndsout$sel_cv %>%
      tibble::as_tibble() %>%
      dplyr::select(gene = gene_name, dnds_missense = wmis_cv,
                    `dnds_nonsense+splice` = wnon_cv, dnds_indels = wind_cv,
                    p_missense = pmis_cv, `p_nonsense+splice` = ptrunc_cv,
                    p_indels = pindpos_cv) %>%
      # convert to long format
      tidyr::pivot_longer(cols = -gene, names_to = c("name", "impact"),
                          names_pattern = "^(dnds|p)_(.*)$") %>%
      tidyr::pivot_wider() %>%
      # filter genes
      dplyr::filter(gene %in% genes2plot) %>%
      # factor impacts and genes, recode and factor impacts
      dplyr::mutate(
        impact = factor(impact, levels = names(impact2_pal)),
        gene = factor(gene, levels = genes2plot),
        `dN/dS ratio` = case_when(only_sig_dnds == TRUE & p > sig_p ~ NA,
                                 TRUE ~ dnds))

    # plot
    p_2 <-
      p_dat_2 %>%
      ggplot(aes(x = gene, y = `dN/dS ratio`, fill = impact)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = impact2_pal) +
      geom_hline(yintercept = 1, color = "grey") +
      scale_x_discrete(guide = guide_axis(angle = -90)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.title.x = element_blank())

    # 3) % of mutant cells
    #    This is based on some assumptions:
    #    1. The calculation will be restricted to samples listed in the
    #       driver_density_sample_ids argument.
    #    2. The values are averages across all donors weighted by duplex
    #       coverage.
    #    3. Values are bounded between sum(duplexVAFs) and sum(2*duplexVAFs), as
    #       explained in Martincorena et al 2018 for standard VAFs.
    #    4. To account for VAFs (clone sizes), dN/dS ratios for significant
    #       genes will be recalculated weighting mutations by their frequency.
    #       This should provide an unbiased point estimate for the purpose of
    #       this analysis, but CI95% cannot be calculated with geneci.

    # mutation table with each mutation represented N times (N = times_called)
    annot_muts_n <-
      annot_muts %>%
      tidyr::uncount(times_called) %>%
      # optionally filter samples based on driver density + filter genes
      dplyr::filter(length(driver_density_sample_ids) == 0 |
                      sampleID %in% driver_density_sample_ids) %>%
      # add a numeric suffix to sampleID to avoid removal of redundant mutations
      dplyr::ungroup() %>%
      dplyr::mutate(sampleID = paste(sampleID, dplyr::row_number()))

    # calculate VAF-weighted dN/dS ratios to estimate the density of driver
    # mutations per cell
    dndsaux <-
      annot_muts_n %>%
      dndscv::dndscv(gene_list = target_genes,
                     max_muts_per_gene_per_sample = Inf,
                     max_coding_muts_per_sample = Inf,
                     constrain_wnon_wspl = TRUE,
                     mingenecovs = 0,
                     dc = gene2dc / mean(gene2dc),
                     kc = newkc,
                     outmats = TRUE,
                     onesided = TRUE,
                     cv = scores,
                     refdb = RefCDS,
                     maxcovs = 10,
                     use_indel_sites = use_indel_sites)

    # cell fractions
    cell_fractions <-
      dndsaux$annotmuts %>%
      # remove numeric suffix
      dplyr::mutate(sampleID = gsub(" [0-9]+$", "", sampleID)) %>%
      # add duplex vaf and cell fraction info
      dplyr::left_join(
        annot_muts %>% dplyr::distinct(sampleID, chr, pos, ref, mut,
                                       duplex_vaf, cellfraction)) %>%
      # unique mutations
      dplyr::distinct() %>%
      # get number of samples
      dplyr::ungroup() %>%
      dplyr::mutate(n_samples = dplyr::n_distinct(sampleID)) %>%
      # filter
      dplyr::filter(gene %in% genes2plot,
                    impact %in% c("Missense", "Nonsense", "Essential_Splice",
                                  "no-SNV")) %>%
      dplyr::group_by(gene, impact) %>%
      dplyr::summarise(hi = sum(cellfraction) / unique(n_samples),
                       lo = sum(duplex_vaf) / unique(n_samples),
                       .groups = "drop") %>%
      tidyr::pivot_longer(cols = c("hi", "lo"), names_to = "cf_bound",
                          values_to = "cf") %>%
      tidyr::complete(gene, impact, cf_bound, fill = list(cf = 0))

    # driver fractions
    driver_fractions <-
      dndsaux$sel_cv %>%
      # filter genes
      dplyr::filter(gene_name %in% genes2plot) %>%
      tidyr::pivot_longer(cols = c("wmis_cv", "wnon_cv", "wspl_cv", "wind_cv"),
                          values_to = "w_cv", names_to = "impact") %>%
      dplyr::select(gene = gene_name, impact, w_cv) %>%
      dplyr::mutate(impact = recode(impact,
                                    wmis_cv = "Missense",
                                    wnon_cv = "Nonsense",
                                    wind_cv = "no-SNV",
                                    wspl_cv = "Essential_Splice")) %>%
      # approximate fraction per driver gene
      dplyr::mutate(
        driver_fraction = (w_cv - 1) / w_cv,
        driver_fraction = ifelse(driver_fraction < 0, 0, driver_fraction))

    # combine cell and driver fractions
    p_dat_3 <-
      # estimated driver density per gene (% cells w driver mutation)
      dplyr::full_join(driver_fractions, cell_fractions) %>%
      dplyr::mutate(driver_density = cf * driver_fraction * 100) %>%
      # sum per gene
      dplyr::group_by(gene, cf_bound) %>%
      dplyr::summarise(driver_density = sum(driver_density),
                       .groups = "drop") %>%
      tidyr::pivot_wider(names_from = cf_bound,
                         values_from = driver_density) %>%
      # factor genes
      dplyr::mutate(gene = factor(gene, levels = genes2plot))

    # plot
    p_3 <-
      p_dat_3 %>%
      ggplot(aes(x = gene)) +
      geom_col(aes(y = hi), fill = "indianred3") +
      geom_col(aes(y = lo), fill = "white", width = 1) +
      scale_x_discrete(guide = guide_axis(angle = -90)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(y = "% cells w/ driver mutation") +
      theme(axis.title.x = element_blank())

    # 4) Unbiased BAM VAFs for the non-synonymous mutations in each gene
    p_dat_4 <-
      annot_muts %>%
      # filter gene
      dplyr::filter(gene %in% genes2plot) %>%
      dplyr::mutate(
        gene = factor(gene, levels = genes2plot),
        bam_vaf_corr = (bam_mut - times_called) / (bam_cov - duplex_cov),
        bam_vaf_corr = ifelse(bam_vaf_corr < 0, NA, bam_vaf_corr)) %>%
      dplyr::filter(gene %in% genes2plot,
                    impact %in% c("Missense", "Nonsense",
                                  "Essential_Splice"),
                    !is.na(bam_vaf_corr)) %>%
      dplyr::mutate(
        bam_vaf_corr_cut = cut(bam_vaf_corr, include.lowest = TRUE,
                               breaks = c(0, 1e-4, 1e-3, 1e-2, 0.05, 1)))

    # create vaf bin palette
    vaf_bin_colours <-
      c("#482677FF", "#2D708EFF", "#29AF7FFF", "#B8DE29FF", "red")
    vaf_bin_pal <-
      p_dat_4 %>%
      dplyr::filter(!is.na(bam_vaf_corr_cut)) %>%
      dplyr::pull(bam_vaf_corr_cut) %>%
      unique() %>%
      sort() %>%
      {setNames(vaf_bin_colours, .)}

    # plot
    p_4 <-
      p_dat_4 %>%
      dplyr::mutate(
        bam_vaf_corr_cut = factor(bam_vaf_corr_cut,
                                  levels = rev(names(vaf_bin_pal)))) %>%
      ggplot(aes(x = gene, fill = bam_vaf_corr_cut)) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = vaf_bin_pal) +
      scale_x_discrete(guide = guide_axis(angle = -90)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(fill = "BAM VAF bin")

    # optionally apply log scale
    if (logscale) {
      p_3 <- p_3 + scale_y_log10()
    }

    # return full plot
    p_1 / p_2 / p_3 / p_4

  }