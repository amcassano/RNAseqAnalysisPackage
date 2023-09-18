#' getDESeq
#'
#' This function will return the DESeq data object from which normalized counts and pairwise comparison results can be generated
#'
#' @param countsmatrix data frame, cleaned up raw counts dataframe from RNA seq experiment
#' @param metadata data frame, cleaned up metadata table from RNA seq experiment
#' @param expt_design formula, experimental design
#'
#' @return DESeqDataSet object
#' @export
#'
#' @examples
#' getDESeq(raw, meta, ~Group)
getDESeq <- function(countsmatrix, metadata, expt_design) {
  # create DESeq Data Set
  deseq_dataset <- DESeq2::DESeqDataSetFromMatrix(countData = countsmatrix, colData = metadata, design = expt_design)

  # run DESeq2 package on data object.
  dseq <- DESeq2::DESeq(deseq_dataset)

  return(dseq)
}
