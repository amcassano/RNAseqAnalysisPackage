#' Filter rLog DataFrame
#' This function will filter out genes that are not significantly different in any of the pairwise comparisons of interest.
#' This chunk will create a filtered dataframe only containing genes that reach significance in your specified pairwise comparisons.
#' This is helpful for removing genes that are not relevant to your conditions to minimize distracting noise.
#'
#' @param p_thresh p value threshold, can be the same as your threshold for filtering pairwise comps but can be higher to preserve more genes
#' @param fc_thresh l2fc threshold, similar to above
#' @param rlogdf your annotated rlog_df
#' @param ... any number of annotated, unfiltered, pairwise comparison results
#'
#' @return dataframe of rlog counts, filtered
#' @export
#'
#' @examples
#' rlog_filtering(0.05, 0.5, rlog_df, one_vs_two, two_vs_three, three_vs_four)
rlog_filtering <- function(p_thresh, fc_thresh, rlogdf, ...){
  # set up inputs for interior loop
  comparisons_list <- list(...)
  num_comparisons <- length(comparisons_list)

  #interior loop
  interior_loop <- function(total_comps, comps_left, comp_list, genelist_temp = NULL){
    list_index <- total_comps - comps_left + 1
    temp_df <- comp_list[[list_index]]
    temp_df <- analyzeRNA::signif_deg(result = temp_df, padj_cutoff = p_thresh, l2fc_cutoff = fc_thresh)
    mgis <- dplyr::select(temp_df, MGI_Symbol)
    if (comps_left == 1) {
      genelist_temp <- dplyr::bind_rows(mgis, genelist_temp)
      genelist_temp <- dplyr::distinct(genelist_temp, genelist_temp$MGI_Symbol, .keep_all = TRUE)
      return(genelist_temp)
    }
    else if (comps_left == total_comps) {
      genelist_temp <- mgis
      comps_left <- comps_left - 1
      interior_loop(total_comps, comps_left, comp_list, genelist_temp)
    }
    else{
      genelist_temp <- dplyr::bind_rows(mgis, genelist_temp)
      genelist_temp <- dplyr::distinct(genelist_temp, genelist_temp$MGI_Symbol, .keep_all = TRUE)
      comps_left <- comps_left - 1
      interior_loop(total_comps, comps_left, comp_list, genelist_temp)
    }
  }
  genelist <- interior_loop(num_comparisons, num_comparisons, comparisons_list)

  rlogdf <- dplyr::filter(rlogdf, MGI_Symbol %in% genelist$MGI_Symbol)
  return(rlogdf)
}
