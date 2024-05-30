#' Get PCA object
#'
#' Creates the PCA object that is used in the PCA plot as well as other things like the Screeplot
#'
#' @param norm_df dataframe of normalized counts (rlog df)
#' @param metadat metadata dataframe
#' @param varianceRemove 	Remove this percentage of variables based on low variance.
#'
#' @return pca_obj
#' @export
#'
#' @examples
#' get_pca_obj(rlog_df, meta, 0.1)
get_pca_obj <- function(
  norm_df,
  metadat,
  varianceRemove = 0.9){
  #fix up norm_df for the PCA function
  pca_input <- dplyr::select(norm_df, -c(MGI_Desc))
  pca_input <- tibble::remove_rownames(pca_input)
  pca_input <- dplyr::distinct(pca_input, MGI_Symbol, .keep_all = TRUE)
  pca_input <- dplyr::filter(pca_input, MGI_Symbol != "")
  pca_input <- dplyr::filter(pca_input, !is.na(MGI_Symbol))
  pca_input <- tibble::column_to_rownames(pca_input, var = "MGI_Symbol")

  # create PCA object
  pca_obj <- PCAtools::pca(mat = pca_input,
                           metadata = metadat,
                           removeVar = varianceRemove)

  return(pca_obj)
}
