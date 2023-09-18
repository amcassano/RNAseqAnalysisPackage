#' Annotate Biomart
#'
#' add human readable gene IDs, descriptions, and gene biotypes to the data
#'
#' @param res dataframe built from DESeq results, must have ensembl gene ID in 'GeneID' column
#' @param gm biomart genemap, contains GeneID and MGI_symbol columns
#'
#' @return results object now with annotations
#' @export
#'
#' @examples
#' annotate_biomart(tol_vs_naive, gmap)
#' annotate_biomart(tol_vs_naive, gmap, c("MGI_Symbol", "MGI_Desc"))
annotate_biomart <- function(res, gm, gmattri = c("MGI_Symbol", "MGI_Desc", "GeneType")) {
  # join the annotations with the results
  res <- dplyr::left_join(res, gm, by = c("GeneID" = "ensembl_gene_id"))

  res <- dplyr::select(res, -`genemap$ensembl_gene_id`)

  # rename the columns
  colList <- BiocGenerics::colnames(res)
  newCols <- purrr::map(colList, rename_columns)
  BiocGenerics::colnames(res) <- NULL
  BiocGenerics::colnames(res) <- newCols

  res <- reorder_columns(res, gmattri)

  # return annotated results table
  return(res)
}
