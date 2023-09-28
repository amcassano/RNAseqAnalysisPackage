#' Get GSEA
#'
#' @param degs annotated, pairwise comparison dataframe, *not* the significant genes
#' @param ontol string, either "BP", "MF", "CC", or "ALL", defaults to ALL
#'
#' @return GSEA object
#' @export
#'
#' @examples
#' getGSEAobj(tol_vs_naive, ontol = "BP")
getGSEAobj <- function(degs, ontol = "ALL"){
  #set up the gene set - get fold change and gene names
  geneset <- degs$Log2FoldChange
  degs <- tibble::rownames_to_column(degs, var = "GeneID")
  names(geneset) <- degs$GeneID
  geneset <- stats::na.omit(geneset)
  geneset <- sort(geneset, decreasing = TRUE)

  #call gene set function
  gsea <- clusterProfiler::gseGO(geneList = geneset,
                      ont = ontol,
                      keyType = "ENSEMBL",
                      minGSSize = 1,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      verbose = F,
                      OrgDb = "org.Mm.eg.db",
                      pAdjustMethod = "fdr")

  return(gsea)
}

