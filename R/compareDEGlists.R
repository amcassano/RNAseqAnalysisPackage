#' Compare two lists of DEGs
#'
#' @param deg_df1 dataframe of DEGs from a pairwise comparison
#' @param deg_df2 dataframe of DEGs from a pairwise comparison
#' @param conditions list of 4 strings, c("numerator for df1", "denominator for df1", "numerator for df2", "denominator for df2")
#' @param l2fc_cutoff number, the cutoff for up or down based on log2 fold change
#' @param pval_cutoff pvalue cut off for which points to exclude
#' @param export should the results also be exported to CSVs
#'
#' @return dataframe of genes with l2fc, pval for both conditions and the mgi symbol and description
#' @export
#'
#' @examples
#' compareDEGlists(one_vs_two, three_vs_two, conditions = c("one", "two", "three", "two"), l2fc_cutoff = 1.2, pval_cutoff = 0.25)
compareDEGlists<- function(deg_df1, deg_df2, conditions = c("numerator1", "denominator1", "num2", "denom2"),
                           l2fc_cutoff = 0.8, pval_cutoff = 0.45, export = FALSE){
  #filter out NA & insignificant P values
  deg_df1 <- dplyr::filter(deg_df1, !is.na(Adj_P_Value), Adj_P_Value < pval_cutoff)
  deg_df2 <- dplyr::filter(deg_df2, !is.na(Adj_P_Value), Adj_P_Value < pval_cutoff)

  #filter out empty or NA MGI symbols
  deg_df1 <- tibble::rownames_to_column(deg_df1, var = "ENSMBL")
  deg_df2 <- tibble::rownames_to_column(deg_df2, var = "ENSMBL")
  deg_df1 <- dplyr::filter(deg_df1, MGI_Symbol != "", !is.na(MGI_Symbol))
  deg_df2 <- dplyr::filter(deg_df2, MGI_Symbol != "", !is.na(MGI_Symbol))

  #condition labels
  cond1 <- paste(conditions[1], "/", conditions[2], sep = "")
  cond2 <- paste(conditions[3], "/", conditions[4], sep = "")

  #set labels
  if (conditions[2] == conditions[4]) {
    bothUp <- paste("Up in both ", conditions[1], " & ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    oneUpOnly <- paste("Up in ", conditions[1], " & down in ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    twoUpOnly <- paste("Up in ", conditions[3], " & down in ", conditions[1], " (vs. ", conditions[2], ")", sep = "")
    bothDown <- paste("Down in both ", conditions[1], " & ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    notsig <- paste("Not up or down in either ", conditions[1], " or ", conditions[3], " (vs. ", conditions[2], ")", sep = "")
    full_title <- paste(conditions[1], "and", conditions[3], "vs", conditions[2])
  }
  else{
    bothUp <- paste("Up in both", cond1, "&", cond2, sep = " ")
    oneUpOnly <- paste("Up in", cond1, "& down in", cond2, sep = " ")
    twoUpOnly <- paste("Up in", cond2, "& down in", cond1, sep = " ")
    bothDown <- paste("Down in both", cond1, "&", cond2, sep = " ")
    notsig <- paste("Not up or down", "in either", cond1, "or", cond2, sep = " ")
    full_title <- paste(cond1, "and", cond2)
  }
  #combine data frames
  combinedData <- dplyr::inner_join(x = deg_df1, y = deg_df2, by = "ENSMBL", suffix = c(".df1", ".df2"))
  combinedData <- tibble::column_to_rownames(combinedData, var = "ENSMBL")
  combinedData <- dplyr::select(combinedData, c(MGI_Symbol = MGI_Symbol.df1, MGI_Desc = MGI_Desc.df1,
                                                Log2FoldChange.df1, Adj_P_Value.df1,
                                                Log2FoldChange.df2, Adj_P_Value.df2))

  #add change directions
  combinedData$ChangeDir <- notsig
  combinedData$ChangeDir[combinedData$Log2FoldChange.df1 > l2fc_cutoff & combinedData$Log2FoldChange.df2 > l2fc_cutoff] <- bothUp
  combinedData$ChangeDir[combinedData$Log2FoldChange.df1 < -l2fc_cutoff & combinedData$Log2FoldChange.df2 < -l2fc_cutoff] <- bothDown
  combinedData$ChangeDir[combinedData$Log2FoldChange.df1 > l2fc_cutoff & combinedData$Log2FoldChange.df2 < -l2fc_cutoff] <- oneUpOnly
  combinedData$ChangeDir[combinedData$Log2FoldChange.df1 < -l2fc_cutoff & combinedData$Log2FoldChange.df2 > l2fc_cutoff] <- twoUpOnly

  #filter out not changed genes
  combinedData <- dplyr::filter(combinedData, ChangeDir != notsig)

  #tidy column names
  combinedData <- dplyr::rename_with(.data = combinedData, .fn = ~ gsub(".df1", paste(" (", cond1, ")", sep = ""), .x))
  combinedData <- dplyr::rename_with(.data = combinedData, .fn = ~ gsub(".df2", paste(" (", cond2, ")", sep = ""), .x))

  if(export) {
    #export
    analyzeRNA::export_genelist(dplyr::filter(combinedData, ChangeDir == bothUp), paste(full_title, "up in both"))
    analyzeRNA::export_genelist(dplyr::filter(combinedData, ChangeDir == bothDown), paste(full_title, "down in both"))
    analyzeRNA::export_genelist(dplyr::filter(combinedData, ChangeDir == oneUpOnly), paste(full_title, "up in", conditions[1], "v", conditions[2]))
    analyzeRNA::export_genelist(dplyr::filter(combinedData, ChangeDir == twoUpOnly), paste(full_title, "up in", conditions[3], "v", conditions[4]))
  }

  #return
  return(combinedData)
}
