#' Sample to Sample Distances
#'
#' @param normobj DESeqTransform object
#' @param annotationcolors list, list of colors for use in heat map
#' @param metadat data frame, meta-data data frame of 1 variable
#' @param clusterrows boolean, defaults to False, should there be row clustering in heatmap
#' @param clustercols boolean, defaults to False,should there be column clustering in heatmap
#' @param fontsize numeric, defaults to 9, font size for heatmap
#'
#' @return heatmap
#' @export
#'
#' @examples
#' sample_distances(rlog_norm,
#'                  list(Condition = c("Acute Rejection" = "#0404b9", Tolerant = "#389538")),
#'                  meta)
sample_distances <- function(normobj, annotationcolors, metadat, clusterrows = FALSE, clustercols = FALSE, fontsize = 9) {
  distances <- stats::dist(t(SummarizedExperiment::assay(normobj)))
  print(distances)
  samp_dist_matrix <- as.matrix(distances)

  # plot the sample distances as a heatmap
  hmap <- pheatmap::pheatmap(samp_dist_matrix, clustering_distance_rows = distances, clustering_distance_cols = distances,
                             cluster_cols = F, cluster_rows = F,
                             annotation_col = metadat, annotation_colors = annotationcolors,
                             main = "Sample to Sample Distances", fontsize_row = fontsize, fontsize_col = fontsize)
  return(hmap)
}
