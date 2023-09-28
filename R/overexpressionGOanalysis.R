#' Over-expression GO Term Analysis
#'
#' @param deg_df dataframe
#' @param pval_cutoff number
#' @param l2fc_cutoff number
#' @param ontol string; either "ALL", "BP", "MF", or "CC"
#'
#' @return data frame
#' @export
#'
#' @examples
#' getOE_analysis(one_vs_two)
getOE_analysis <- function(deg_df, pval_cutoff = 0.05, l2fc_cutoff = 0.5, ontol = "ALL"){
  #remove any rows with no Gene ID or with duplicates (there shouldn't be any but it'd mess everything up)
  deg_df <- tibble::rownames_to_column(deg_df, var = "GeneID")
  deg_df <- dplyr::distinct(deg_df, GeneID, .keep_all = TRUE)

  allGenes <- deg_df$GeneID

  sigDEGs <- dplyr::filter(deg_df, Adj_P_Value <= pval_cutoff)
  sigDEGs <- dplyr::filter(sigDEGs, abs(Log2FoldChange) >= l2fc_cutoff)

  sigGenes <- sigDEGs$GeneID

  OE_GO <- clusterProfiler::enrichGO(gene = sigGenes,
                                   universe = allGenes,
                                   keyType = "ENSEMBL",
                                   OrgDb = "org.Mm.eg.db",
                                   ont = ontol,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)


  return(OE_GO)
}
