#' Heatmap from Data Frame
#'
#' @param plottitle string, title of the plot
#' @param genelist data frame, contains mgi symbols for genes to be plotted
#' @param met data frame, meta data
#' @param counts_df data frame, contains normalized counts
#' @param clusterRows boolean, defaults to TRUE, should the rows cluster by similarity
#' @param gaps list of numbers, defaults to empty list, dictates where there would be gaps if desired between conditions
#' @param height number, defaults to 10, dictates how tall the cells in the heatmap will be
#' @param annocolors list of strings, list of colors for use in labeling
#' @param colorscale colorRampPalette funciton call - sets the color scale, defaults to red (hi) white blue (low)
#'
#' @return heatmap plot
#' @export
#'
#' @examples
#' heatmapFromDF("DEGs Naive vs Tol", ranked_naive_vs_tol$MGI_Symbol,
#'               metadata, rlog_df, c("blue", "red"))
heatmapFromDF <- function(plottitle,
                          genelist,
                          met,
                          counts_df,
                          annocolors,
                          clusterRows = TRUE,
                          gaps = c(),
                          height = 10,
                          colorscale = grDevices::colorRampPalette(c("blue3", "white", "red2")) (100)) {
  #create data frame for heatmap that combines the genes to plot and counts
  heatmap_data <- dplyr::select(counts_df, -MGI_Desc)
  heatmap_data <-  dplyr::filter(heatmap_data, MGI_Symbol %in% genelist$MGI_Symbol)
  heatmap_data <-  tibble::rownames_to_column(heatmap_data, var = "GeneID")

  create_heatmap(title = plottitle,
                 genes_and_counts = heatmap_data,
                 met = met,
                 annocolors = annocolors,
                 clusterRows = clusterRows,
                 gaps = gaps,
                 height = height,
                 colorscale = colorscale)
}

