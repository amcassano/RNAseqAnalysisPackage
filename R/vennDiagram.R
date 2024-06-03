#' Plot Venn Diagram
#'
#' Plot a venn diagram showing the number of shared DEGs between 2+ comparisons
#' Can be any number of lists but for best Venn Diagram
#'
#' @param degs_list list of named lists, lists of MGI symbols from significant DEGs of pairwise comparisons, should be named with the label
#' @param plt_title string, name of the plot
#' @param edgecolor list of strings, colors for the circles of the venn diagram, should be in the same order as the lists
#' @param showintersects boolean, if true the plot will be interactive and will show the genes in the intersects
#' @param colorscale list of two string colors
#'
#' @return ggplot venn diagram
#' @export
#'
#' @examples
#' plot_venndiagram(list("Treg" = ranked_TregNaive$MGI_Symbol, "Tol" = ranked_TolNaive$MGI_Symbol))
#' plot_venndiagram(list("Treg" = ranked_TregNaive$MGI_Symbol, "Tol" = ranked_TolNaive$MGI_Symbol, "Rej" = ranked_RejNaive$MGI_Symbol),
#'  edgecolor = c("#ffe280", "#000000", "#b0001b"),
#'  show intersects = TRUE,
#'  colorscale = c("#fde9ec", "#540f0f"))
#'
plot_venndiagram <- function(degs_list, plt_title = "", edgecolor = c("#000000"), showintersects = FALSE, colorscale = c("#FFFFFF", "#8c9fbf")) {
  venn <-
    ggVennDiagram::ggVennDiagram(degs_list,
                                 lty = 1,
                                 edge_size = 1.25,
                                 set_size = 4.25,
                                 show_intersect = showintersects,
                                 label_size = 4,
                                 label = "count",
                                 label_alpha = 0.4) +
    ggplot2::scale_fill_gradient(low = colorscale[1], high = colorscale[2]) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2)) +
    ggplot2::scale_color_manual(values = edgecolor)

  venn <- venn + ggplot2::labs(title = plt_title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14))

  return(venn)
}

#' Get Venn Diagram Overlaps
#'
#' @param degs_list list of named lists, lists of MGI symbols from significant DEGs of pairwise comparisons, should be named with the label
#'
#' @return dataframe of overlaps between comparisions
#' @export
#'
#' @examples
#' venn_overlaps(list("Treg" = ranked_TregNaive$MGI_Symbol, "Tol" = ranked_TolNaive$MGI_Symbol))
#' venn_overlaps(list("Treg" = ranked_TregNaive$MGI_Symbol, "Tol" = ranked_TolNaive$MGI_Symbol, "Rej" = ranked_RejNaive$MGI_Symbol))
venn_overlaps <- function(degs_list){
  venn_obj <- ggVennDiagram::Venn(degs_list)
  overlaps <- ggVennDiagram::process_region_data(venn_obj)
  overlaps <- tidyr::unnest(overlaps, item)

  return(overlaps)
}
