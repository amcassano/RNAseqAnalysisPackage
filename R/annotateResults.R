#' Annotate Biomart
#'
#' add human readable gene IDs, descriptions, and gene biotypes to the data
#'
#' @param res dataframe built from DESeq results, must have ensembl gene ID in 'GeneID' column
#' @param gm biomart genemap, contains GeneID and MGI_symbol columns
#'
#' @return results object now with annotations
#' @export

annotate_biomart <- function(res, gm) {

  # join the annotations with the results
  res <- res %>% dplyr::left_join(gm, by = c("GeneID" = "ensembl_gene_id"))

  # rename the columns
  colList <- BiocGenerics::colnames(res)
  newCols <- purrr::map(colList, rename_columns)
  BiocGenerics::colnames(res) <- NULL
  BiocGenerics::colnames(res) <- newCols

  # return annotated results table
  return(res)
}
