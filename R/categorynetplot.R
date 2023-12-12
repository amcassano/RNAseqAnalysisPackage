#' Create Category Netplot
#'
#' @param go_obj go object, result of either getOE_analysis or get_GSEA
#' @param deg_df dataframe
#' @param pval_cutoff number
#' @param l2fc_cutoff number
#' @param categories_to_show number; defaults to 30
#' @param seed number; defaults to 2, for set seed for reproducibility
#' @param shape string; one of "star", "circle", "gem", "dh", 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'.
#'
#' @return category netplot
#' @export
#'
#' @examples
#' categorynetplot(goobj, one_vs_two)
categorynetplot <- function(go_obj, deg_df, pval_cutoff = 0.05, l2fc_cutoff = 0.5, categories_to_show = 30, seed = 2, shape = "circle"){
  #remove any rows with no Gene ID or with duplicates (there shouldn't be any but it'd mess everything up)
  deg_df <- tibble::rownames_to_column(deg_df, var = "GeneID")
  deg_df <- dplyr::filter(deg_df, !is.na(GeneID))
  deg_df <- dplyr::distinct(deg_df, GeneID, .keep_all = TRUE)

  sigDEGs <- dplyr::filter(deg_df, Adj_P_Value <= pval_cutoff)
  sigDEGs <- dplyr::filter(sigDEGs, abs(Log2FoldChange) >= l2fc_cutoff)

  sigGenesFC <- sigDEGs$Log2FoldChange
  names(sigGenesFC) <- sigDEGs$GeneID
  sigGenesFC <-  sort(sigGenesFC, decreasing = TRUE)
  go_read <- DOSE::setReadable(go_obj, "org.Mm.eg.db", "ENSEMBL")
  set.seed(seed)
  cnet <-
    enrichplot::cnetplot(go_read,
                         showCategory = categories_to_show,
                         color.params = list(foldChange = sigGenesFC, edge = TRUE, category = "#a71919"),
                         cex.params = list(category_label = 0.9, gene_label = 0.85),
                         layout = shape)

  return(cnet)
}
