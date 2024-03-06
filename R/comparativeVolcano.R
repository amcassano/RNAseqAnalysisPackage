#' Make comparative volcano plot
#'
#' @param deg_df1 dataframe of pairwise comparison results, must be annotated
#' @param deg_df2 dataframe of pairwise comparison results, must be annotated
#' @param plotTitle string, plot title
#' @param conditions list of 4 strings, c("numerator for df1", "denominator for df1", "numerator for df2", "denominator for df2")
#' @param l2fc_cutoff number, where to draw the cutoff lines for log2 fold changes and to determine which points to color/label
#' @param pval_cutoff number, points with P value lower than this in either dataframe will not be plotted
#' @param dircolors list of 5 colors (strings) for both up, both down, up in 1 & down in 2, up in 2 & down in 1, and not changed
#' @param dirshapes list of 5 shapes (numbers) for both up, both down, up in 1 & down in 2, up in 2 & down in 1, and not changed
#' @param overlaps number, indicates max overlaps for the labeling
#' @param boxlabels boolean, defaults to FALSE, if true draws boxes arond labels
#'
#' @return plot of log2 fold change vs log 2 fold change
#' @export
#'
#' @examples comp_volcano(Rej_vs_Tol, Naive_vs_Tol, "Naive and Rejecting comparison", "Rejecting", "Naive", "tolerant")
comp_volcano <- function (deg_df1, deg_df2, plotTitle,
                          conditions = c("numerator1", "denominator1", "num2", "denom2"),
                          l2fc_cutoff = 0.8, pval_cutoff = 0.45,
                          dircolors = c("#a50000", "#00009c", "darkgreen", "purple4", "gray70"),
                          dirshapes = c(24, 25, 23, 23, 21), overlaps = 50, boxlabels = FALSE)
{
  #remove any data points with a Pvalue of NA or of greater than p val cutoff
  deg_df1 <- dplyr::filter(deg_df1, !is.na(Adj_P_Value), Adj_P_Value < pval_cutoff)
  deg_df2 <- dplyr::filter(deg_df2, !is.na(Adj_P_Value), Adj_P_Value < pval_cutoff)

  #filter data sets to get rid of any duplicates
  genelist <- dplyr::select(deg_df1, MGI_Symbol)
  genelist$dupe <- duplicated(genelist)
  genelist <- dplyr::select(genelist, dupe)

  genelist <- tibble::rownames_to_column(genelist, var = "ENSMBL")
  deg_df1 <- tibble::rownames_to_column(deg_df1, var = "ENSMBL")
  deg_df2 <- tibble::rownames_to_column(deg_df2, var = "ENSMBL")

  deg_df1 <- dplyr::left_join(deg_df1, genelist, by = "ENSMBL")
  deg_df2 <- dplyr::left_join(deg_df2, genelist, by = "ENSMBL")

  deg_df1 <- tibble::column_to_rownames(deg_df1, var = "ENSMBL")
  deg_df2 <- tibble::column_to_rownames(deg_df2, var = "ENSMBL")

  deg_df1 <- dplyr::filter(deg_df1, MGI_Symbol != "", !dupe, !is.na(MGI_Symbol))
  deg_df2 <- dplyr::filter(deg_df2, MGI_Symbol != "", !dupe, !is.na(MGI_Symbol))

  #condition labels
  cond1 <- paste(conditions[1], "/", conditions[2], sep = "")
  cond2 <- paste(conditions[3], "/", conditions[4], sep = "")

  #set labels
  if (conditions[2] == conditions[4]) {
    bothUp <- paste("Up in both ", conditions[1], " & ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    oneUpOnly <- paste("Up in ", conditions[1], " & down in ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    twoUpOnly <- paste("Up in ", conditions[3], " & down in ", conditions[1], " (vs. ", conditions[2], ")", sep = "")
    bothDown <- paste("Down in both ", conditions[1], " & ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    notsig <- paste("Not up or down in either ", conditions[1], " or ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    full_title <- paste(conditions[1], "and", conditions[3], "vs", conditions[2])
  }
  else{
    bothUp <- paste("Up in both", cond1, "&", cond2, sep = " ")
    oneUpOnly <- paste("Up in", cond1, "& down in", cond2, sep = " ")
    twoUpOnly <- paste("Up in", cond2, "& down in", cond1, sep = " ")
    bothDown <- paste("Down in both", cond1, "&", cond2, sep = " ")
    notsig <- paste("Not up or down", "in either", cond1, "or", cond2, sep = " ")
    full_title <- paste(cond1, "and", cond2)
  }
  names(dircolors) <- c(bothUp, bothDown, oneUpOnly, twoUpOnly, notsig)
  names(dirshapes) <- c(bothUp, bothDown, oneUpOnly, twoUpOnly, notsig)
  legorder <- c(bothUp, bothDown, oneUpOnly, twoUpOnly, notsig)

  #combine the dataframes into 1 with both fold changes
  combinedData <- dplyr::inner_join(x = deg_df1, y = deg_df2, by = "MGI_Symbol", suffix = c(".c1", ".c2"))
  combinedData <- dplyr::select(combinedData, c(MGI_Symbol, Log2FoldChange.c1, Log2FoldChange.c2))

  #assign labels for if the data point is up/down in both/one condition
  combinedData$ChangeDir <- notsig
  combinedData$ChangeDir[combinedData$Log2FoldChange.c1 > l2fc_cutoff & combinedData$Log2FoldChange.c2 > l2fc_cutoff] <- bothUp
  combinedData$ChangeDir[combinedData$Log2FoldChange.c1 < -l2fc_cutoff & combinedData$Log2FoldChange.c2 < -l2fc_cutoff] <- bothDown
  combinedData$ChangeDir[combinedData$Log2FoldChange.c1 > l2fc_cutoff & combinedData$Log2FoldChange.c2 < -l2fc_cutoff] <- oneUpOnly
  combinedData$ChangeDir[combinedData$Log2FoldChange.c1 < -l2fc_cutoff & combinedData$Log2FoldChange.c2 > l2fc_cutoff] <- twoUpOnly

  #assign gene labels for those up/down in both
  combinedData$signifGeneLabel <- NA
  combinedData$signifGeneLabel[combinedData$ChangeDir == bothUp] <- combinedData$MGI_Symbol[combinedData$ChangeDir == bothUp]
  combinedData$signifGeneLabel[combinedData$ChangeDir == bothDown] <- combinedData$MGI_Symbol[combinedData$ChangeDir == bothDown]
  combinedData$signifGeneLabel[combinedData$ChangeDir == oneUpOnly] <- combinedData$MGI_Symbol[combinedData$ChangeDir == oneUpOnly]
  combinedData$signifGeneLabel[combinedData$ChangeDir == twoUpOnly] <- combinedData$MGI_Symbol[combinedData$ChangeDir == twoUpOnly]

  comboPlot <- ggplot2::ggplot(combinedData,
                               ggplot2::aes(x = Log2FoldChange.c1, y = Log2FoldChange.c2,
                                            color = ChangeDir, shape = ChangeDir, fill = ChangeDir,
                                            label = signifGeneLabel)) +
    ggplot2::geom_point(size = 1.75) +
    ggplot2::scale_color_manual(name = "", values = dircolors, breaks = legorder) +
    ggplot2::scale_fill_manual(name = "", values = dircolors, breaks = legorder) +
    ggplot2::scale_shape_manual(name = "", values = dirshapes, breaks = legorder) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::labs(title = plotTitle,
                  subtitle = paste(full_title," DEGs", sep = ""),
                  caption = paste("\n Log2 FC cutoff: ", l2fc_cutoff,
                                  "\n", "\n",
                                  "Not showing any points with adj. P value of greater than ", pval_cutoff,
                                  sep = ""),
                  x = paste("Log2 Fold Change ", cond1, sep = ""),
                  y = paste("Log2 Fold Change ", cond2, sep = "")) +
    ggplot2::geom_vline(xintercept = c(-1 * l2fc_cutoff,l2fc_cutoff), col = "gray85") +
    ggplot2::geom_hline(yintercept = c(-1 * l2fc_cutoff,l2fc_cutoff), col = "gray85")
  if (boxlabels) {
    comboPlot <- comboPlot + ggrepel::geom_label_repel(min.segment.length = 0,
                                                       size = 2.75,
                                                       point.padding = 0.25,
                                                       box.padding = 0.7,
                                                       label.padding = 0.3,
                                                       color = "white",
                                                       force = 1.35,
                                                       force_pull = 1.2,
                                                       max.overlaps = overlaps,
                                                       show.legend = FALSE,
                                                       na.rm =  TRUE,
                                                       segment.size = 0.35) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent"),
                     plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                     axis.line = ggplot2::element_line("black", 1),
                     aspect.ratio = (10/12), panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, face = "italic"),
                     axis.text = ggplot2::element_text(size = 9),
                     axis.title = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 11),
                     legend.text = ggplot2::element_text(size = 9))
  }
  else{
    comboPlot <- comboPlot + ggrepel::geom_text_repel(min.segment.length = 0,
                                                      size = 2.75,
                                                      point.padding = 0.25,
                                                      box.padding = 0.7,
                                                      force = 1.35,
                                                      force_pull = 1.2,
                                                      max.overlaps = overlaps,
                                                      show.legend = FALSE,
                                                      na.rm =  TRUE,
                                                      segment.size = 0.35) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent"),
                     plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
                     axis.line = ggplot2::element_line("black", 1),
                     aspect.ratio = (10/12), panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, face = "italic"),
                     axis.text = ggplot2::element_text(size = 9),
                     axis.title = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 11),
                     legend.text = ggplot2::element_text(size = 9))
  }
  return(comboPlot)
}
