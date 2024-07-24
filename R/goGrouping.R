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
goGrouping <- function(genes,
                       pval_cutoff = 0.05,
                       l2fc_cutoff = 0.5,
                       ontol = c("BP", "MF", "CC", "ALL", "BP_MF"),
                       lvl = 2) {
  if (ontol != "ALL" | ontol != "BP" | ontol != "MF" | ontol != "CC" | ontol != "BP_MF") {
    return("Please make sure ontol is one of the following: 'BP', 'CC', 'MF', 'ALL' or 'BP_MF'")
  }
  else if (ontol == "ALL") {
    bp <- goGrouping(genes, pval_cutoff, l2fc_cutoff, "BP", lvl)
    mf <- goGrouping(genes, pval_cutoff, l2fc_cutoff, "MF", lvl)
    cc <- goGrouping(genes, pval_cutoff, l2fc_cutoff, "CC", lvl)
    gg <- dplyr::bind_rows(bp, mf, cc, .id = "tbl")
    gg$Ontol <- NULL
    gg$Ontol[gg$tbl == "1"] <- "BP"
    gg$Ontol[gg$tbl == "2"] <- "MF"
    gg$Ontol[gg$tbl == "3"] <- "CC"
    gg <- dplyr::select(gg, -tbl)
    return(gg)
  }
  else if (ontol == "BP_MF") {
    bp <- goGrouping(genes, pval_cutoff, l2fc_cutoff, "BP", lvl)
    mf <- goGrouping(genes, pval_cutoff, l2fc_cutoff, "MF", lvl)
    gg <- dplyr::bind_rows(bp, mf,.id = "tbl")
    gg$Ontol <- NULL
    gg$Ontol[gg$tbl == "1"] <- "BP"
    gg$Ontol[gg$tbl == "2"] <- "MF"
    gg <- dplyr::select(gg, -tbl)
    return(gg)
  }
  else{
    genes <- tibble::rownames_to_column(genes, var = "GeneID")
    genes <- dplyr::filter(genes, Adj_P_Value <= pval_cutoff)
    genes <- dplyr::filter(genes, abs(Log2FoldChange) >= l2fc_cutoff)
    genes <- dplyr::distinct(genes, GeneID)
    gg <- as.data.frame(
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = ontol,
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      )
    )
    gg <- dplyr::filter(gg, Count > 0)
    gg <- dplyr::arrange(gg, desc(Count))
    gg <- tibble::remove_rownames(gg)
    gg <- dplyr::select(gg, -GeneRatio)
    return(gg)
  }
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
goGrouping_rlog <- function(rlogdf,
                            ontol = c("BP", "MF", "CC", "ALL", "BP_MF"),
                            lvl = 2) {
  if (ontol != "ALL" | ontol != "BP" | ontol != "MF" | ontol != "CC" | ontol != "BP_MF") {
    return("Please make sure ontol is one of the following: 'BP', 'CC', 'MF', 'ALL' or 'BP_MF'")
  }
  else if (ontol == "ALL") {
    bp <-  goGrouping_rlog(rlogdf, "BP", lvl)
    mf <- goGrouping_rlog(rlogdf, "MF", lvl)
    cc <- goGrouping_rlog(rlogdf, "CC", lvl)
    gg <- dplyr::bind_rows(bp, mf, cc, .id = "tbl")
    gg$Ontol <- NULL
    gg$Ontol[gg$tbl == "1"] <- "BP"
    gg$Ontol[gg$tbl == "2"] <- "MF"
    gg$Ontol[gg$tbl == "3"] <- "CC"
    gg <- dplyr::select(gg, -tbl)
    return(gg)
  }
  else if (ontol == "BP_MF") {
    bp <-  goGrouping_rlog(rlogdf, "BP", lvl)
    mf <- goGrouping_rlog(rlogdf, "MF", lvl)
    gg <- dplyr::bind_rows(bp, mf, .id = "tbl")
    gg$Ontol <- NULL
    gg$Ontol[gg$tbl == "1"] <- "BP"
    gg$Ontol[gg$tbl == "2"] <- "MF"
    gg <- dplyr::select(gg, -tbl)
    return(gg)
  }
  else{
    genes <- tibble::rownames_to_column(rlogdf, var = "GeneID")
    genes <- dplyr::distinct(genes, GeneID)
    gg <- as.data.frame(
      clusterProfiler::groupGO(
        gene = genes$GeneID,
        OrgDb = "org.Mm.eg.db",
        ont = ontol,
        keyType = "ENSEMBL",
        level = lvl,
        readable = TRUE
      )
    )
    gg <- dplyr::filter(gg, Count > 0)
    gg <- dplyr::arrange(gg, desc(Count))
    gg <- tibble::remove_rownames(gg)
    gg <- dplyr::select(gg, -GeneRatio)
    return(gg)
  }
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
goGrouping_all_levels <- function(genes,
                                  pval_cutoff = 0.05,
                                  l2fc_cutoff = 0.5,
                                  ontol = c("BP", "MF", "CC", "ALL", "BP_MF")) {
  lvl2 <- goGrouping(genes, pval_cutoff, l2fc_cutoff, ontol, lvl = 2)
  lvl3 <- goGrouping(genes, pval_cutoff, l2fc_cutoff, ontol, lvl = 3)
  lvl4 <- goGrouping(genes, pval_cutoff, l2fc_cutoff, ontol, lvl = 4)
  lvl5 <- goGrouping(genes, pval_cutoff, l2fc_cutoff, ontol, lvl = 5)
  lvl6 <- goGrouping(genes, pval_cutoff, l2fc_cutoff, ontol, lvl = 6)
  gg <- dplyr::bind_rows(lvl2, lvl3, lvl4, lvl5, lvl6, .id = "tbl")
  gg$Level <- NULL
  gg$Level[gg$tbl == "1"] <- "2"
  gg$Level[gg$tbl == "2"] <- "3"
  gg$Level[gg$tbl == "3"] <- "4"
  gg$Level[gg$tbl == "4"] <- "5"
  gg$Level[gg$tbl == "5"] <- "6"
  gg <- dplyr::select(gg, -tbl)
  return(gg)
}

#' Get Go Groups from all levels (2-6) from rLogDF
#'
#' @param rlogdf rlog dataframe - ideally filtered to only include significantly changed/relevant genes
#' @param ontol string, "BP", "CC",  "MF", "ALL", or "BP_MF"
#' @return gg, dataframe
#' @export
#'
#' @examples
#' goGrouping_rlog_all_levels(sigrlog, "BP")
goGrouping_rlog_all_levels <- function(rlogdf,
                                       ontol = c("BP", "MF", "CC", "ALL", "BP_MF")) {
  lvl2 <- goGrouping_rlog(rlogdf, ontol, 2)
  lvl3 <- goGrouping_rlog(rlogdf, ontol, 3)
  lvl4 <- goGrouping_rlog(rlogdf, ontol, 4)
  lvl5 <- goGrouping_rlog(rlogdf, ontol, 5)
  lvl6 <- goGrouping_rlog(rlogdf, ontol, 6)
  gg <- dplyr::bind_rows(lvl2, lvl3, lvl4, lvl5, lvl6, .id = "tbl")
  gg$Level <- NULL
  gg$Level[gg$tbl == "1"] <- "2"
  gg$Level[gg$tbl == "2"] <- "3"
  gg$Level[gg$tbl == "3"] <- "4"
  gg$Level[gg$tbl == "4"] <- "5"
  gg$Level[gg$tbl == "5"] <- "6"
  gg <- dplyr::select(gg, -tbl)
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
genelistFromGOGroup <- function(gogroup,
                                rownum = NULL,
                                goid = NULL) {
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
