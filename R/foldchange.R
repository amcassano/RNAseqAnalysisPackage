#' Get Fold Change List
#'
#' @param degs df of pairwise comparison genes, not filtered
#'
#' @return list of fold changes with names of gene ids
#' @export
#'
#' @examples
#' getFoldChanges(one_vs_two)
getFoldChanges <- function(degs){
  geneset <- degs$Log2FoldChange
  degs <- tibble::rownames_to_column(degs, var = "GeneID")
  names(geneset) <- degs$GeneID
  geneset <- stats::na.omit(geneset)
  geneset <- sort(geneset, decreasing = TRUE)

  return(geneset)
}
