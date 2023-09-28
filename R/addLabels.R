#' Add Labels
#'
#' adds labels to each point on a ggplot
#'
#' @param dat data frame
#' @param overlaps numeric, number of overlaps allowed
#'
#' @return geom label repel object
#' @export
#'
#' @examples
#' pca_data <- DESeq2::plotPCA(rlognorm, intgroup = Group, returnData = TRUE)
#' add_labels(pca_data)
add_labels <- function(dat, overlaps = 6) {
  return(ggrepel::geom_label_repel(data = dat,ggplot2::aes(label = BiocGenerics::rownames(dat)), color = "white",segment.color = "black",
                                   min.segment.length = 0, size = 3, point.padding = 1, label.padding = 0.2,
                                   box.padding = 0.8, max.overlaps = overlaps, show.legend = FALSE))
}
