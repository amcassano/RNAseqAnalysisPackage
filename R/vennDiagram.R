#' Plot Venn Diagram
#'
#' Plot a venn diagram showing the number of shared DEGs between 2+ comparisons
#' Can be any number of lists but for best Venn Diagram
#'
#' @param degs_list list of named lists, lists of MGI symbols from significant DEGs of pairwise comparisons, should be named with the label
#' @param plt_title string, name of the plot
#' @param edgecolor list of strings, colors for the circles of the venn diagram, should be in the same order as the lists
#'
#' @return ggplot venn diagram
#' @export
#'
#' @examples
#' plot_venndiagram(list("Treg" = ranked_TregNaive$MGI_Symbol, "Tol" = ranked_TolNaive$MGI_Symbol))
#'
#'
plot_venndiagram <- function(degs_list, plt_title = "", edgecolor = c("#000000")) {
  venn <-
    ggVennDiagram::ggVennDiagram(degs_list, lty = 1, edge_size = 1.25, set_size = 4.25, label_size = 4, label = "count", label_alpha = 0.4) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#8c9fbf") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2)) +
    ggplot2::scale_color_manual(values = edgecolor)

  venn <- venn + ggplot2::labs(title = plt_title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14))

  return(venn)
}
