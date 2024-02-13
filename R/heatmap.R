#' Basic heatmap
#'
#' @param title string, title of the plot
#' @param met data frame, meta data for the experiment
#' @param clusterRows boolean, defaults to TRUE, should the heatmap cluster rows by similarity
#' @param gaps list of numbers, defaults to empty list, dictates where there would be gaps if desired between conditions
#' @param height number, defaults to 10, dictates how tall the cells in the heatmap will be
#' @param annocolors list of strings, list of colors for use in labeling
#' @param genes_and_counts data frame, contains genes of interest and the corresponding normalized read counts
#' @param colorscale list of 3 colors sets the color scale, defaults to red (hi) white blue (low)
#'
#' @return pheatmap plot
#' @export
#'
#' @examples
#' create_heatmap("Oxphos Genes", oxhposgenes_andcounts, metadata, c("Tol" = "blue", "Rej" = "green"))
create_heatmap <- function(title,
                           genes_and_counts,
                           met,
                           annocolors,
                           clusterRows = TRUE,
                           gaps = c(),
                           height = 10,
                           colorscale = c("blue4", "white", "red3")) {
  if (length(genes_and_counts) == 0) {return("No genes match")}
  else{
    heatmap <-
      pheatmap::pheatmap(
        genes_and_counts[3:ncol(genes_and_counts)],
        color = grDevices::colorRampPalette(colorscale) (100),
        cluster_rows = clusterRows, cluster_cols = FALSE, clustering_distance_rows = "correlation",
        gaps_col = gaps, annotation_col = met, annotation_colors = annocolors,
        scale = "row", labels_row = genes_and_counts[, 2], labels_col = c(rep("", ncol(genes_and_counts))),
        fontsize = height + 1.5, fontsize_row = height + 0.25, fontsize_col = 10,
        cellheight = height, cellwidth = height + 3,
        main = title, border_color = NA, annotation_legend = TRUE)
    return(heatmap)
  }
}
