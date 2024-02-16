#' Add Labels
#'
#' adds labels to each point on a ggplot
#'
#' @param dat data frame
#' @param overlaps numeric, number of overlaps allowed
#' @param boxed boolean, if there should be a box drawn under the label
#'
#' @return geom label repel object
#' @export
#'
#' @examples
#' pca_data <- DESeq2::plotPCA(rlognorm, intgroup = Group, returnData = TRUE)
#' add_labels(pca_data)
add_labels <- function(dat, overlaps = 6, boxed = TRUE) {
  if (boxed) {
    return(ggrepel::geom_label_repel(
      data = dat,
      ggplot2::aes(label = BiocGenerics::rownames(dat)),
      color = "white", segment.color = "black", segment.size = 0.35,
      min.segment.length = 0,
      size = 3,
      point.padding = 0.25, label.padding = 0.5, box.padding = 0.75,
      max.overlaps = overlaps,
      show.legend = FALSE,
      force_pull = 1.2, force = 1.35, na.rm = TRUE
    ))
  }
  else {
    return(ggrepel::geom_text_repel(
      data = dat,
      ggplot2::aes(label = BiocGenerics::rownames(dat)),
      color = "white", segment.color = "black",  segment.size = 0.35,
      min.segment.length = 0,
      size = 3,
      point.padding = 0.25, box.padding = 0.75,
      max.overlaps = overlaps,
      show.legend = FALSE,
      force = 1.35, force_pull = 1.2,
      na.rm = TRUE
    ))
      }
  }
