#' Plot Gene Counts
#'
#' @param mgi string, MGI symbol for the gene of interest
#' @param norm_df data frame, transformed DESeq object as a data frame
#' @param label_samples boolean, defaults to FALSE, if true, adds labels for samples while avoiding overlaps
#' @param overlaps_allowed number, defaults to 10, increase to allow more overlaps in labeling
#' @param yaxistozero boolean, defaults to FALSE, if true the Y axis will
#' @param comparisons list of length 2 lists, each length 2 list is the comparisons to use for stats, defaults to an empty list
#' @param kw_hjust integer, dictates the horizontal position of the kruskal wallis value
#' @param kw_vjust integer, dictates the vertical position of the kruskal wallis value
#' @param yaxis string, label for the Y axis
#' @param metadat dataframe, meta data
#' @param plot_aes dataframe, containing labels, colors, shapes, and fills for the data points being plotted
#'
#' @return plot of normalized gene counts
#' @export
#'
#' @examples
#' plot_genecounts("Tcf7", vst_df, metadata, treatment_aes, comparisons = stat_comparisons)
#'
plot_genecounts <- function(mgi, norm_df, metadat, plot_aes, comparisons = list(),
                            label_samples = FALSE, overlaps_allowed = 10, yaxistozero = FALSE, kw_hjust = 0.2, kw_vjust = 1,
                            yaxis = "rLog Normalized Reads"){

  norm_plotdata <- data.frame(t(norm_df[norm_df$MGI_Symbol == as.character(mgi), ]))
  norm_plotdata <- dplyr::select(norm_plotdata, dplyr::contains("ENS"))
  BiocGenerics::colnames(norm_plotdata) <- NULL
  BiocGenerics::colnames(norm_plotdata) <- "norm"
  actual_data <- as.numeric(dplyr::tally(norm_plotdata) - 2)
  gene_info <- data.frame(t(dplyr::slice_head(norm_plotdata, n = 2)))
  norm_plotdata <- dplyr::slice_tail(norm_plotdata, n = actual_data)
  norm_plotdata <- dplyr::bind_cols(metadat, norm_plotdata)
  norm_plotdata$norm <- as.numeric(norm_plotdata$norm)
  yaxismax <- max(norm_plotdata$norm)
  yaxismin <- min(norm_plotdata$norm)

  # create ggplot
  countPlot <-
    ggplot2::ggplot(norm_plotdata, ggplot2::aes(x = Group, y = norm, color = Group, shape = Group, fill = Group)) +
    ggplot2::ggtitle(gene_info$MGI_Symbol, subtitle = gene_info$MGI_Desc) +
    ggplot2::scale_color_manual(labels = plot_aes$labels, values = plot_aes$outlineID, name = "Group: ") +
    ggplot2::scale_shape_manual(labels = plot_aes$labels, values = plot_aes$shapeID, name = "Group: ") +
    ggplot2::scale_fill_manual(labels = plot_aes$labels, values = plot_aes$fillID, name = "Group: ") +
    ggplot2::xlab(BiocGenerics::paste("Group")) +
    ggplot2::ylab(yaxis) +
    ggplot2::scale_y_continuous(n.breaks = 8)

  # display plot with global theme and labels
  if (label_samples) {
    countPlot <- countPlot +
      global_theme() +
      ggplot2::geom_point(size = 3.5, stroke = 2) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
      add_labels(norm_plotdata, overlaps_allowed)
  }
  else {
    countPlot <- countPlot +
      global_theme() +
      ggplot2::geom_jitter(size = 3.5, width = 0.1, height = 0, stroke = 2) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_line(color = "gray95", size = 0.25, linetype = 2),
                     panel.grid.minor.y = ggplot2::element_line(color = "gray95", size = 0.25,linetype = 2))
  }

  # edit y axis
  if (yaxistozero) {countPlot <- countPlot + ggplot2::expand_limits(y = 0)}
  else if ((yaxismin - 5) < 0) {countPlot <- countPlot + ggplot2::expand_limits(y = c(0, yaxismax + 3))}
  else {countPlot <- countPlot + ggplot2::expand_limits(y = c(yaxismin - 5, yaxismax + 3))}

  # add comparison stats
  countPlot <- countPlot +
    ggpubr::stat_compare_means(comparisons = comparisons, label.y.npc = "bottom", vjust = 0.5, label = "p.signif",  size = 5.5, step.increase = 0.04) +
    ggpubr::stat_compare_means(label.y.npc = "bottom", show.legend = FALSE, size = 3.2, vjust = kw_vjust, hjust = kw_hjust)

  return(countPlot)
}
