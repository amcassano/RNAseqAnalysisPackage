#' Pairwise Differentitally Expressed Genes comparison results
#'
#' Get the results of a deseq object
#'
#' @param deseq_obj DESeqdata object
#' @param numerator string, name of the conditon of interest for comparison, cannot contain spaces
#' @param denominator string, name of baseline condition for comparison
#' @param contrast_str string, defaults to "Group", the comparison of interest for pairwise comparison
#' @param tidy_result boolean, defaults to TRUE, dictates if the results are returned as a dataframe or obj
#' @param log2FCshrink boolean, defaults to false, dictates if ashr lfc shrinking formula is used
#'
#' @return DESeq Results object, or a dataframe
#' @export
#'
#' @examples
#' pairwiseDEGresults("Tolerant", "Naive", dseqobj)
#' pairwiseDEGresults("Rejecting", "Tolerant", dseqobj, log2FCshrink = TRUE)
#'
pairwiseDEGresults <-
  function(numerator, denominator, deseq_obj, contrast_str = "Group", tidy_result = TRUE, log2FCshrink = FALSE) {
    # create comparison for results
    con <- c(contrast_str, numerator, denominator)
    if (log2FCshrink) {
      res <- DESeq2::lfcShrink(dds = deseq_obj, contrast = con, type = "ashr", quiet = TRUE)
      if (tidy_result) {
        res <- BiocGenerics::as.data.frame(res)
        res <- tibble::rownames_to_column(res, var = "GeneID")
        # res <- stats::na.omit(res)
      }
    } else {
      res <- DESeq2::results(deseq_obj, con, tidy = tidy_result)
      if (tidy_result) {
        res <- dplyr::rename(res, "GeneID" = "row")
        # res <- stats::na.omit(res)
      }
    }
    return(res)
  }
