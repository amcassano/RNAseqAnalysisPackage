#' Volcano plot
#'
#' creates a volcano plot showing the relationship between log2foldchange and p value for all genes in a pairwise comparison
#'
#' @param deg_df data frame, pairwise comparison results, should contain at minimum Log2FC MGI symbol, and adjusted p value
#' @param dircolors list of strings, colors for significant up, down, and insignificant, defaults to red, blue and grey
#' @param dirshapes list of numbers, shapes for points up, down, and insignificant, defaults to circles for all
#' @param l2fc_cutoff number, log2 fold change threshold for significance, defaults to 1
#' @param pval_cutoff number, adjusted p value threshold for significance, defaults to 0.01
#' @param plotTitle string, plot title
#' @param cond1 string, name of condition 1
#' @param cond2 string, name of condition 2
#'
#' @return ggplot volcano plot
#' @export
#'
#' @examples
#' volcano_plot(tol_vs_naive, "Tolerant vs. Naive", "Tolerant", "Naive")
#' volcano_plot(tol_vs_naive, "Tolerant vs. Naive", "Tolerant", "Naive", pval_cutoff = 0.005)
volcano_plot <- function(deg_df, plotTitle, cond1, cond2, l2fc_cutoff = 1, pval_cutoff = 0.01,
                         dircolors = c("#a50000", "#00009c", "gray70"), dirshapes = c(19, 19, 19)) {
  # make labels for genes
  dirUp <- paste("Upregulated in", cond1, "", sep = " ")
  dirDown <- paste("Upregulated in", cond2, "", sep = " ")
  notsig <- paste("Not significantly", "differentially expressed", sep = "\n")
  names(dircolors) <- c(dirUp, dirDown, notsig)
  names(dirshapes) <- c(dirUp, dirDown, notsig)

  # save dataframe of differentially expressed genes as a new data frame for manipulation, filter out unannotated genes
  volcanoData <- dplyr::distinct(deg_df, MGI_Symbol, .keep_all = TRUE)

  # keep only padj & log2 fold change columns, remove NA rows
  volcanoData <- dplyr::select(volcanoData, c(Adj_P_Value, Log2FoldChange, MGI_Symbol))
  volcanoData <- dplyr::filter(volcanoData, !is.na(Adj_P_Value))
  volcanoData <- dplyr::filter(volcanoData, !is.na(Log2FoldChange))
  volcanoData <- dplyr::filter(volcanoData, !is.na(MGI_Symbol))

  # add change direction to significantly up or significantly down genes
  volcanoData$ChangeDir <- notsig
  volcanoData$ChangeDir[volcanoData$Log2FoldChange > l2fc_cutoff & volcanoData$Adj_P_Value < pval_cutoff] <- dirUp
  volcanoData$ChangeDir[volcanoData$Log2FoldChange < (-1 * l2fc_cutoff) & volcanoData$Adj_P_Value < pval_cutoff] <- dirDown

  # add column with MGI symbol only if significant
  volcanoData$signifGeneLabel <- NA
  volcanoData$signifGeneLabel[volcanoData$ChangeDir != notsig] <- volcanoData$MGI_Symbol[volcanoData$ChangeDir != notsig]

  # make volcano plot
  volcanoPlot <-
    ggplot2::ggplot(
      volcanoData,
      ggplot2::aes(
        x = Log2FoldChange, y = -log10(Adj_P_Value), color = ChangeDir,
        shape = ChangeDir, fill = ChangeDir, label = signifGeneLabel
      )
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(name = "", values = dircolors) +
    ggplot2::scale_fill_manual(name = "", values = dircolors) +
    ggplot2::scale_shape_manual(name = "", values = dirshapes) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::labs(
      title = plotTitle, subtitle = paste(cond1, "vs.", cond2),
      caption = paste("Log2 FC cutoff: ", l2fc_cutoff, "\n Adj. P val. cutoff: ", pval_cutoff, sep = ""),
      x = "Log2 Fold Change", y = "-Log10(Adj. P Value)"
    ) +
    ggplot2::geom_vline(xintercept = c(-1 * l2fc_cutoff, l2fc_cutoff), col = "#767676") +
    ggplot2::geom_hline(yintercept = -log10(pval_cutoff), col = "#767676") +
    ggrepel::geom_label_repel(
      min.segment.length = 0, size = 2.25, point.padding = 0.25, label.padding = 0.2, box.padding = 0.25,
      max.overlaps = 70, show.legend = FALSE, fill = "white"
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      axis.line = ggplot2::element_line("black", 1),
      aspect.ratio = (10 / 12),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, face = "italic"),
      axis.text = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9)
    )

  return(volcanoPlot)
}
