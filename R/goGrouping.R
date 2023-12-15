#' Group Go
#' Get a list of GO terms at a designated level that are represented by a list of genes input
#' @param genes dataframe, list of DEGs, must be annotated and have Log2FoldChange and Adj_P_Value columns
#' @param pval_cutoff number, defaults to 0.05
#' @param l2fc_cutoff number, defaults to 0.5
#' @param ontol string, "BP", "CC", or "MF" only
#' @param lvl number, defaults to 2, the level of GO to query
#'
#' @return gg, dataframe
#' @export
#'
#' @examples
#' goGrouping(one_vs_two, 0.01, 0.8, "BP", 4)
goGrouping <- function(genes, pval_cutoff = 0.05, l2fc_cutoff = 0.5, ontol = "BP", lvl = 2){
  genes <- tibble::rownames_to_column(genes, var = "GeneID")
  genes <- dplyr::filter(genes, Adj_P_Value <= pval_cutoff)
  genes <- dplyr::filter(genes, abs(Log2FoldChange) >= l2fc_cutoff)
  genes <- dplyr::distinct(genes, GeneID)

  gg <- clusterProfiler::groupGO(gene = genes$GeneID,
                                 OrgDb = "org.Mm.eg.db",
                                 ont = ontol,
                                 keyType = "ENSEMBL",
                                 level = lvl,
                                 readable = TRUE)

  gg <- as.data.frame(gg)
  gg <- dplyr::filter(gg, Count > 0)
  gg <- dplyr::arrange(gg, desc(Count))
  return(gg)
}

#' Gene List from Go Grouping
#'
#' @param gogroup dataframe, output of goGrouping function
#' @param rownum number, which row to return genes from
#'
#' @return dataframe, MGI symbols
#' @export
#'
#' @examples
#' genelistFromGOGroup(goGrouping(one_vs_two), 2)
genelistFromGOGroup <- function(gogroup, rownum){
  gogroup <- dplyr::slice(gogroup, rownum)
  allgenes <- gogroup$geneID
  listgenes <- stringr::str_split(allgenes, "/")
  listgenes <- listgenes[[1]]
  listgenes <- as.data.frame(listgenes)
  listgenes <- dplyr::rename(listgenes, MGI_Symbol = listgenes)
  return(listgenes)
}
