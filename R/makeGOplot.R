#' Make GO Plot
#'
#' @param GO_obj go object, result of either getOE_analysis or get_GSEA
#' @param plot_title string
#' @param plot_type string, one of "dot", "ridge", "cluster",
#' @param categories_to_show number; defaults to 30
#' @param fontsize number; defaults to 9
#' @param splitby string; defaults to NULL
#'
#' @return a plot or a string
#' @export
#'
#' @examples
#' makeGOplot(onevstwo_GSEA_BP, "One Vs Two Dotplot", "dot")
makeGOplot <- function(GO_obj, plot_title, plot_type, categories_to_show = 30, fontsize = 9){
  if (stringi::stri_cmp_equiv(plot_type, "dot", strength = 1)) {
      dotplot <-
        enrichplot::dotplot(GO_obj, color = "p.adjust",  size = "Percentage",
                                     showCategory = categories_to_show,
                                     label_format = 75, font.size = fontsize,
                                     title = plot_title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    return(dotplot)
  }
  else if(stringi::stri_cmp_equiv(plot_type, "ridge", strength = 1)) {
    ridge <-
      enrichplot::ridgeplot(GO_obj, label_format = 75, showCategory = categories_to_show, fill = "p.adjust") +
      ggplot2::labs(x = "Enrichment Distribution", title = plot_title) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize), plot.title = ggplot2::element_text(hjust = 0.5))

    return(ridge)
  }
  else if(stringi::stri_cmp_equiv(plot_type, "cluster", strength = 1)) {
    emap <-
      enrichplot::emapplot(GO_obj, showCategory = categories_to_show, label_format = 75, font.size = fontsize) +
      ggplot2::labs(title = plot_title)

    return(emap)
  }
  else {
    return(paste("Please set plot_type to be one of the following", "dot", "ridge", "cluster", sep = "\n"))
  }
}
