#' GSEA Gene ONtology
#'
#' @param degs annotated, pairwise comparison dataframe, *not* the significant genes
#' @param ontol string, either "BP", "MF", "CC", or "ALL", defaults to ALL
#'
#' @return GSEA object
#' @export
#'
#' @examples
#' go_GSEA(tol_vs_naive, ontol = "BP")
go_GSEA <- function(degs, ontol = "ALL"){
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
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      verbose = F,
                      OrgDb = "org.Mm.eg.db",
                      pAdjustMethod = "fdr")

  return(gsea)
}

#' KEGG GSEA
#'
#' @param degs annotated, pairwise comparison dataframe, *not* the significant genes
#'
#' @return GSEA object
#' @export
#'
#' @examples
#' kegg_GSEA(tol_vs_naive)
kegg_GSEA <- function(degs){
  degs <- stats::na.omit(degs)
  degs <- dplyr::distinct(degs, EntrezID, .keep_all = TRUE)

  #set up the gene set - get fold change and gene names
  geneset <- degs$Log2FoldChange
  names(geneset) <- degs$EntrezID
  geneset <- stats::na.omit(geneset)
  geneset <- sort(geneset, decreasing = TRUE)

  #call gene set function
  gsea <- clusterProfiler::gseKEGG(geneList = geneset,
                                   organism = "mmu",
                                   minGSSize = 10,
                                   maxGSSize = 1000,
                                   pvalueCutoff = 0.05)

  return(gsea)
}
