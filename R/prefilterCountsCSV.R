#' Pre Filter Counts CSV
#' Filters out genes with counts below the threshold set, improves signal to noise ratio
#'
#' @param rawcounts raw counts dataframe
#' @param numberofsamples number, the total number of samples in the counts dataframe
#' @param avg_count_minimum number, defaults to 200, higher number will be more stringent filtering, lower number will be more generous
#'
#' @return counts dataframe
#' @export
#'
#' @examples
#' prefilterCounts(raw_counts, 17)
#' prefilterCounts(raw_counts, 13, 150)
prefilterCounts <- function(rawcounts, numberofsamples, avg_count_minimum = 200){
  row_min <- avg_count_minimum * numberofsamples
  rawcounts$row_max <- apply(rawcounts, 1, max)
  rawcounts$row_sum <- apply(rawcounts, 1, sum)
  rawcounts <- dplyr::filter(rawcounts, row_sum - row_max * 2 > row_min)
  rawcounts <- dplyr::select(rawcounts, -c(row_sum, row_max))
  return(rawcounts)
}
