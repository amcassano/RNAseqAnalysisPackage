#' Create Category Netplot
#'
#' @param go_obj go object, result of either getOE_analysis or get_GSEA
#' @param deg_df dataframe
#' @param pval_cutoff number
#' @param l2fc_cutoff number
#' @param ontol string; either "ALL", "BP", "MF", or "CC"
#' @param categories_to_show number; defaults to 30
#' @param seed number; defaults to 2, for set seed for reproducibility
#'
#' @return category netplot
#' @export
#'
#' @examples
#' categorynetplot(goobj, one_vs_two)
categorynetplot <- function(go_obj, deg_df, pval_cutoff = 0.05, l2fc_cutoff = 0.5, categories_to_show = 30, ontol = "ALL", seed = 2){
  #remove any rows with no Gene ID or with duplicates (there shouldn't be any but it'd mess everything up)
  deg_df <- tibble::rownames_to_column(deg_df, var = "GeneID")
  deg_df <- dplyr::filter(deg_df, !is.na(GeneID))
  deg_df <- dplyr::distinct(deg_df, GeneID, .keep_all = TRUE)

  allGenes <- deg_df$GeneID

  sigDEGs <- dplyr::filter(deg_df, Adj_P_Value <= pval_cutoff)
  sigDEGs <- dplyr::filter(sigDEGs, abs(Log2FoldChange) >= l2fc_cutoff)

  sigGenesFC <- sigDEGs$Log2FoldChange
  names(sigGenesFC) <- sigDEGs$GeneID
  set.seed(seed)
  cnet <-
    enrichplot::cnetplot(go_obj, showCategory = categories_to_show, foldChange = sigGenesFC)

  return(cnet)
}
