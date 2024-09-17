#' Pre Filter Counts CSV
#' Filters out genes with counts below the threshold set, improves signal to noise ratio
#'
#' @param rawcounts raw counts dataframe
#' @param numberofsamples number, the total number of samples in the counts dataframe
#' @param avg_count_minimum number, defaults to 200, higher number will be more stringent filtering, lower number will be more generous
#' @param stringency number, the denominator for the second filtering. defaults to 2 which means that genes where the top half of the counts meet the avg will be kept. if 3 then only 1/3 of counts need to meet it etc
#'
#' @return counts dataframe
#' @export
#'
#' @examples
#' prefilterCounts(raw_counts, 17)
#' prefilterCounts(raw_counts, 13, 150)
prefilterCounts <- function(rawcounts, numberofsamples, avg_count_minimum = 200, stringency = 2){

  frac_samples <- ceiling(numberofsamples / stringency)
  row_min <- avg_count_minimum * (numberofsamples - 1)
  row_min_frac <- avg_count_minimum * frac_samples
  rawcounts <- tibble::column_to_rownames(rawcounts, var = "Geneid")
  rawcounts$rowmax <- apply(rawcounts, 1, max)
  rawcounts$rowsum <- apply(rawcounts, 1, sum)
  rawcounts <- tibble::rownames_to_column(rawcounts, var = "rowid")
  rawcounts <- rawcounts %>% dplyr::rowwise() %>%
    dplyr::mutate(topfrac = sum(sort(dplyr::c_across(
      !tidyselect::starts_with("row")
    ), decreasing = TRUE)[1:frac_samples])) %>%
    dplyr::ungroup()
  rawcounts <- dplyr::filter(rawcounts,
                              rowsum - rowmax * 2 >= row_min |
                                topfrac >= row_min_frac)
  rawcounts <- dplyr::select(rawcounts, -c(rowsum, rowmax, topfrac))
  rawcounts <- tibble::column_to_rownames(rawcounts, var = "rowid")
  return(rawcounts)
}
