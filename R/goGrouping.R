#' Group Go
#' Get a list of GO terms at a designated level that are represented by a list of genes input
#' @param genes dataframe, list of DEGs, must be annotated and have Log2FoldChange and Adj_P_Value columns
#' @param pval_cutoff number, defaults to 0.05
#' @param l2fc_cutoff number, defaults to 0.5
#' @param ontol string, "BP", "CC",  "MF", "ALL", or "BP_MF"
#' @param lvl number, defaults to 2, the level of GO to query
#'
#' @return gg, dataframe
#' @export
#'
#' @examples
#' goGrouping(one_vs_two, 0.01, 0.8, "BP", 4)
goGrouping <- function(genes, pval_cutoff = 0.05, l2fc_cutoff = 0.5, ontol = c("BP", "MF", "CC", "ALL", "BP_MF"), lvl = 2) {
  genes <- tibble::rownames_to_column(genes, var = "GeneID")
  genes <- dplyr::filter(genes, Adj_P_Value <= pval_cutoff)
  genes <- dplyr::filter(genes, abs(Log2FoldChange) >= l2fc_cutoff)
  genes <- dplyr::distinct(genes, GeneID)
  if (ontol == "ALL") {
    gg <- dplyr::bind_rows(
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "BP",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      ),
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "MF",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      ),
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "CC",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      )
    )
  }
  else if(ontol == "BP_MF"){
    gg <- dplyr::bind_rows(
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "BP",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      ),
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "MF",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      )
    )
  }
  else{
    gg <- clusterProfiler::groupGO(
      gene = genes$GeneID,
      OrgDb = "org.Mm.eg.db",
      ont = ontol,
      keyType = "ENSEMBL",
      level = lvl,
      readable = TRUE
    )
  }
  gg <- as.data.frame(gg)
  gg <- dplyr::filter(gg, Count > 0)
  gg <- dplyr::arrange(gg, desc(Count))
  gg <- tibble::remove_rownames(gg)
  return(gg)
}

#' Group GO from rLog DF
#'
#' @param rlogdf rlog dataframe - ideally filtered to only include significantly changed/relevant genes
#' @param ontol string, "BP", "CC",  "MF", "ALL", or "BP_MF"
#' @param lvl number, defaults to 2, the level of GO to query
#'
#' @return gg, dataframe
#' @export
#'
#' @examples
#' goGrouping_rlog(sig_rlog, "BP", 3)
goGrouping_rlog <- function(rlogdf, ontol = c("BP", "MF", "CC", "ALL", "BP_MF"), lvl = 2) {
  genes <- tibble::rownames_to_column(rlogdf, var = "GeneID")
  genes <- dplyr::distinct(genes, GeneID)
  if (ontol == "ALL") {
    gg <- dplyr::bind_rows(
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "BP",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      ),
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "MF",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      ),
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = "CC",
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      )
    )
  }
  else{
    gg <- clusterProfiler::groupGO(
      gene = genes$GeneID,
      OrgDb = "org.Mm.eg.db",
      ont = ontol,
      keyType = "ENSEMBL",
      level = lvl,
      readable = TRUE
    )
  }

  gg <- as.data.frame(gg)
  gg <- dplyr::filter(gg, Count > 0)
  gg <- dplyr::arrange(gg, desc(Count))
  gg <- tibble::remove_rownames(gg)
  return(gg)
}

#' Get Go Groups from all levels (2-6)
#'
#' @param genes dataframe, list of DEGs, must be annotated and have Log2FoldChange and Adj_P_Value columns
#' @param pval_cutoff number, defaults to 0.05
#' @param l2fc_cutoff number, defaults to 0.5
#' @param ontol string, "BP", "CC",  "MF", "ALL", or "BP_MF"
#'
#' @return gg, dataframe
#' @export
#'
#' @examples
#' goGrouping_all_levels(one_vs_two, 0.01, 0.8, "BP")
goGrouping_all_levels <- function(genes, pval_cutoff = 0.05, l2fc_cutoff = 0.5, ontol = c("BP", "MF", "CC", "ALL", "BP_MF")){
  gogroups <- dplyr::bind_rows(
    goGrouping(genes = genes,
               pval_cutoff = pval_cutoff,
               l2fc_cutoff = l2fc_cutoff,
               ontol = ontol,
               lvl = 2),
    goGrouping(genes = genes,
               pval_cutoff = pval_cutoff,
               l2fc_cutoff = l2fc_cutoff,
               ontol = ontol,
               lvl = 3),
    goGrouping(genes = genes,
               pval_cutoff = pval_cutoff,
               l2fc_cutoff = l2fc_cutoff,
               ontol = ontol,
               lvl = 4),
    goGrouping(genes = genes,
               pval_cutoff = pval_cutoff,
               l2fc_cutoff = l2fc_cutoff,
               ontol = ontol,
               lvl = 5),
    goGrouping(genes = genes,
               pval_cutoff = pval_cutoff,
               l2fc_cutoff = l2fc_cutoff,
               ontol = ontol,
               lvl = 6)
  )
  return(gogroups)
}

#' Get Go Groups from all levels (2-6) from rLogDF
#'
#' @param rlogdf rlog dataframe - ideally filtered to only include significantly changed/relevant genes
#' @param ontol string, "BP", "CC",  "MF", "ALL", or "BP_MF"
#' @return gg, dataframe
#' @export
#'
#' @examples
#' goGrouping_rlog_all_levels(one_vs_two, 0.01, 0.8, "BP")
goGrouping_rlog_all_levels <- function(rlogdf, ontol = c("BP", "MF", "CC", "ALL", "BP_MF")){
  gogroups <- dplyr::bind_rows(
    goGrouping_rlog(rlogdf = rlogdf,
                    ontol = ontol,
                    lvl = 2),
    goGrouping_rlog(rlogdf = rlogdf,
                    ontol = ontol,
                    lvl = 3),
    goGrouping_rlog(rlogdf = rlogdf,
                    ontol = ontol,
                    lvl = 4),
    goGrouping_rlog(rlogdf = rlogdf,
                    ontol = ontol,
                    lvl = 5),
    goGrouping_rlog(rlogdf = rlogdf,
                    ontol = ontol,
                    lvl = 6)
  )
  return(gogroups)
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
genelistFromGOGroup <- function(gogroup, rownum = NULL, goid = NULL) {
  if (!is.null(rownum)) {
    gogroup <- dplyr::slice(gogroup, rownum)
  } else if (!is.null(goid)) {
    gogroup <- dplyr::filter(gogroup, stringr::str_detect(ID, goid))
  } else {
    gogroup <- dplyr::slice(gogroup, 1)
  }
  allgenes <- gogroup$geneID
  listgenes <- stringr::str_split(allgenes, "/")
  listgenes <- listgenes[[1]]
  listgenes <- as.data.frame(listgenes)
  listgenes <- dplyr::rename(listgenes, MGI_Symbol = listgenes)
  return(listgenes)
}
