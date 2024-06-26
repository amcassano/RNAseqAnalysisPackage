#' Filter and Rank Significant Differentially Expressed Genes
#'
#' @description This function \itemize{
#' \item keeps only differentially expressed genes that are significant
#' \item adds an indicator of if genes in the experimental group are upregulated or downregulated compared to baseline
#' \item and calculates a change metric which will combine direction, Log2 Fold Change, and adjusted Pvalue}
#'
#' @param result dataframe contains \itemize{
#' \item gene annotations
#' \item log2 FC and adjusted p value}
#' @param padj_cutoff double, cutoff level for adjusted P vals to keep, defaults to 0.05
#' @param l2fc_cutoff double, cutoff level for log2 fold change, defaults to 0.8 (abs value)
#' @param entrezid boolean, if entrez id is included
#'
#' @return ranked, a sorted and filtered data frame containing \itemize{
#' \item annotations
#' \item change direction
#' \item fold change
#' \item adjusted p value
#' \item change metric (log2 fold change * - log10 padj)}
#' @export
#'
#' @examples
#' signif_deg(anno_tolvsnaive)
#' signif_deg(anno_rejvstol, padj_cutoff = 0.01, l2fc_cutoff = 0.5)
signif_deg <- function(result, padj_cutoff = 0.05, l2fc_cutoff = 0.8, entrezid = FALSE) {
  ranked <- tibble::rownames_to_column(result, "GeneID")

  # create new df, removing any samples with NA as the value for FC or Pvals
  ranked <- tidyr::drop_na(ranked, c(Log2FoldChange, Adj_P_Value))

  # keep only results that have an adjusted pvalue of less than 0.05 and FC more than 0.8
  ranked <- dplyr::filter(ranked, Adj_P_Value <= padj_cutoff)
  ranked <- dplyr::filter(ranked, abs(Log2FoldChange) >= l2fc_cutoff)
  ranked <- dplyr::filter(ranked, MGI_Symbol != "")

  # add the change direction and change metric
  ranked <- dplyr::mutate(ranked, change_dir = ifelse(Log2FoldChange < 0, "Down", "Up"))
  ranked <- dplyr::mutate(ranked, change_metric = Log2FoldChange * -log10(Adj_P_Value))

  # keep only the relevant columns & arrange in order of change metric
  if(entrezid){
    ranked <- dplyr::select(ranked,GeneID, MGI_Symbol, MGI_Desc, EntrezID, Log2FoldChange, Adj_P_Value,
                            "ChangeDirection" = change_dir,"ChangeMetric" = change_metric)
    ranked <- dplyr::arrange(ranked, dplyr::desc(abs(ChangeMetric)))
  }
  else {
    ranked <- dplyr::select(ranked,GeneID, MGI_Symbol, MGI_Desc, Log2FoldChange, Adj_P_Value,
                          "ChangeDirection" = change_dir,"ChangeMetric" = change_metric)
    ranked <- dplyr::arrange(ranked, dplyr::desc(abs(ChangeMetric)))
  }
  BiocGenerics::rownames(ranked) <- ranked$GeneID
  ranked <- dplyr::select(ranked, -GeneID)
  return(ranked)
}
