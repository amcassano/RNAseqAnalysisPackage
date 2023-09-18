#' GSEA Dotplot
#'
#' @param gsea_obj gene set object
#' @param plottitle string, plot title
#' @param categories_to_show number, defaults to 25
#' @param font_size number, defaults to 8
#' @param split boolean, defaults to FALSE
#'
#' @return dotplot
#' @export
#'
#' @examples
#' makeGSEAdotplot(tolvsnaiveGSEA, "")
makeGSEAdotplot <- function(gsea_obj, plottitle, categories_to_show = 25, font_size = 8, split = FALSE){

  if (split) {
    dotplot <- enrichplot::dotplot(gsea_obj,
                       showCategory = categories_to_show,
                       split = ".sign",
                       label_format = 75, font.size = font_size,
                       title = plottitle)  + DOSE::facet_grid(.~.sign)
  }
  else{
    dotplot <- enrichplot::dotplot(gsea_obj,
                       showCategory = categories_to_show, split = NULL,
                       label_format = 75, font.size = font_size,
                       title = plottitle) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  return(dotplot)
}
