#' GSEA ridgeplot
#'
#' @param gsea_obj gene set object
#' @param plottitle string
#' @param cats_to_show  number
#' @param font_size number
#'
#' @return ridgeplot
#' @export
#'
#' @examples
#' makeGSEAridgeplot(tolvsnaiveGSEA, "")
makeGSEAridgeplot <- function(gsea_obj, plottitle, cats_to_show = 25, font_size = 8){
  ridge <-
    enrichplot::ridgeplot(gsea_obj, label_format = 75) +
    ggplot2::labs(x = "Enrichment Distribution", title = plottitle) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = font_size), plot.title = ggplot2::element_text(hjust = 0.5))

  return(ridge)

}
