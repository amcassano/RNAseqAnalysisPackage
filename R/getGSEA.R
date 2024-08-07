#' GSEA Gene Ontology
#'
#' @param degs annotated, pairwise comparison dataframe, *not* the significant genes
#' @param ontol string, either "BP", "MF", "CC", or "ALL", defaults to ALL. note must be ALL if used for REVIGO, can be any choice for use with other plots
#' @param pval_cutoff number, defaults to 0.05, p value cutoff value
#'
#' @return GSEA object
#' @export
#'
#' @examples
#' go_GSEA(tol_vs_naive, ontol = "BP")
go_GSEA <- function(degs, ontol = "ALL", pval_cutoff = 0.05) {
  FCs <- getFoldChanges(degs)
  # call gene set function
  gsea <- clusterProfiler::gseGO(
    geneList = FCs,
    ont = ontol,
    keyType = "ENSEMBL",
    minGSSize = 10,
    maxGSSize = 750,
    pvalueCutoff = pval_cutoff,
    verbose = F,
    OrgDb = "org.Mm.eg.db",
    pAdjustMethod = "fdr"
  )

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
kegg_GSEA <- function(degs) {
  degs <- stats::na.omit(degs)
  degs <- dplyr::distinct(degs, EntrezID, .keep_all = TRUE)

  # set up the gene set - get fold change and gene names
  geneset <- degs$Log2FoldChange
  names(geneset) <- degs$EntrezID
  geneset <- stats::na.omit(geneset)
  geneset <- sort(geneset, decreasing = TRUE)

  # call gene set function
  gsea <- clusterProfiler::gseKEGG(
    geneList = geneset,
    organism = "mmu",
    minGSSize = 10,
    maxGSSize = 1000,
    pvalueCutoff = 0.05
  )

  return(gsea)
}
