#' Make Gene Set Enrichment Plot
#'
#' @param GSE_obj go object, result of either go_GSEA
#' @param plot_title string
#' @param plot_type string, one of "dot", "ridge", "cluster", "heat", "upset","tree", "traditional"
#' @param categories_to_show number; defaults to 30
#' @param label_length number; defaults to 75, how long to let a label be before wrapping
#' @param fontsize number; defaults to 9
#' @param foldChangeList list, defaults to null, fold changes with names attached - needed for coloring by FC
#' @param legend_loc string; defaults to "bottom" where the legend should be placed
#' @param seed number; defaults to 5, for reproducibility
#'
#' @return a plot or a string
#' @export
#'
#' @examples
#' makeGSEplot(onevstwo_GSEA_BP, "One Vs Two Dotplot", "dot")
makeGSEplot <-
  function(GSE_obj,
           plot_title,
           plot_type,
           categories_to_show = 30,
           fontsize = 9,
           label_length = 75,
           foldChangeList = NULL,
           legend_loc = "bottom",
           seed = 5) {
    set.seed(seed = seed)
    if (stringi::stri_cmp_equiv(plot_type, "dot", strength = 1)) {
      dotplot <-
        enrichplot::dotplot(
          GSE_obj,
          showCategory = categories_to_show,
          label_format = label_length,
          font.size = fontsize,
          title = plot_title
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize * 1.35),
          legend.position = legend_loc
        )

      return(dotplot)
    }
    else if (stringi::stri_cmp_equiv(plot_type, "ridge", strength = 1)) {
      ridge <-
        enrichplot::ridgeplot(
          GSE_obj,
          label_format = label_length,
          showCategory = categories_to_show,
          fill = "p.adjust"
        ) +
        ggplot2::labs(x = "Enrichment Distribution", title = plot_title) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = fontsize),
          plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize * 1.35),
          legend.position = legend_loc
        )

      return(ridge)
    }
    else if (stringi::stri_cmp_equiv(plot_type, "heat", strength = 1)) {
      GSE_read <- DOSE::setReadable(GSE_obj, "org.Mm.eg.db", "ENSEMBL")
      heat <-
        enrichplot::heatplot(
          GSE_read,
          showCategory = categories_to_show,
          label_format = label_length,
          foldChange = foldChangeList
        ) +
        ggplot2::labs(title = plot_title) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize * 1.35),
          legend.position = legend_loc
        )

      return(heat)
    }
    else if (stringi::stri_cmp_equiv(plot_type, "upset", strength = 1)) {
      upset <-
        enrichplot::upsetplot(GSE_obj, n = categories_to_show) +
        ggplot2::labs(title = plot_title) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize * 1.35),
          legend.position = legend_loc
        )

      return(upset)
    }
    else if (stringi::stri_cmp_equiv(plot_type, "cluster", strength = 1)) {
      GSE_obj <- enrichplot::pairwise_termsim(GSE_obj)
      emap <-
        enrichplot::emapplot(
          GSE_obj,
          showCategory = categories_to_show,
          node_label = "category",
          shadowtext = FALSE,
          cluster.params = list(
            cluster = TRUE,
            method = stats::kmeans,
            n = ceiling(sqrt(categories_to_show)),
            legend = FALSE,
            label_style = "ggforce",
            label_words_n = 4,
            label_format = label_length
          ),
          cex.params = list(
            category_node = 0.85,
            category_label = 0.85,
            label_group = 1.2,
            line = 0.4
          )
        ) +
        ggplot2::labs(title = plot_title) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize * 1.35),
          legend.position = legend_loc
        )
      return(emap)
    }
    else if (stringi::stri_cmp_equiv(plot_type, "tree", strength = 1)) {
      GSE_obj <- enrichplot::pairwise_termsim(GSE_obj)
      tree <-
        enrichplot::treeplot(
          GSE_obj,
          showCategory = categories_to_show,
          label_format = label_length,
          offset.params = list(
            bar_tree = ggplot2::rel(1.5),
            tiplab = ggplot2::rel(0.75),
            hexpand = 0.35
          ),
          cluster.params = list(
            method = "ward.D",
            n = ceiling(sqrt(categories_to_show)),
            label_words_n = 4,
            label_format = label_length
          )
        ) +
        ggplot2::labs(title = plot_title) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize * 1.35),
          legend.position = legend_loc
        )

      return(tree)
    }
    else if (stringi::stri_cmp_equiv(plot_type, "traditional", strength = 1)) {
      gse <-
        enrichplot::gseaplot2(
          GSE_obj,
          geneSetID = 1:categories_to_show,
          base_size = fontsize,
          title = plot_title,
          pvalue_table = TRUE
        )
      return(gse)
    }
    else {
      return(
        paste(
          "Please set plot_type to be one of the following:",
          "dot",
          "ridge",
          "cluster",
          "heat",
          "upset",
          "tree",
          "or traditional",
          sep = " "
        )
      )
    }
  }
