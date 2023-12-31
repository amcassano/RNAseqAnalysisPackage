#' reorder Columns
#'
#' @param annotated_df dataframe, annotated results data frame, must include GeneID column
#' @param move2front list of strings, which columns to move to the front of the data frame, defaults to mgi symbol, description and entrez ID
#'
#' @return data frame with columns in correct order and geneID set to row name
#' @export
#'
#' @examples
#' reorder_columns(anno_tolvsnaive)
#' reorder_columns(anno_rejvstol, move2front = c("MGI_Symbol", "MGI_Symbol"))
reorder_columns <- function(annotated_df, move2front = c("MGI_Symbol", "MGI_Desc", "EntrezID")) {
  #set GeneID as the row name
  BiocGenerics::rownames(annotated_df) <- annotated_df$GeneID

  #remove gene ID column, move desired columns to the front
  annotated_df <- dplyr::select(annotated_df, -GeneID)
  annotated_df <- dplyr::relocate(annotated_df, move2front)

  return(annotated_df)
}
