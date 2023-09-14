#' Export Gene List
#'
#' @param tab data frame, contains ranked DEGs, the output of signif_deg
#' @param filename string, name to save the exported file as, avoid whitespace
#' @param direction string ("up", "down", or NULL), which genes to export
#'
#' @return CSV file exported to CWD
#' @export
#'
#' @examples
#' export_genelist(ranked_naive_vs_tol, "Naive-vs-tol")
#' export_genelist(ranked_naive_vs_tol, "Naive-vs-tol", direction = "up")
#' export_genelist(ranked_naive_vs_tol, "Naive-vs-tol", direction = "down")
export_genelist <- function(tab, filename, direction = NULL) {

  # if direction is "up", export only those genes with a + change metric
  if(is.null(direction)){
    utils::write.csv(tab, file = BiocGenerics::paste(filename, ".csv", sep = ""))
  }
  else if(stringi::stri_cmp_equiv(direction, "up", strength = 1)){
    upreg <- tab %>%
      dplyr::filter(ChangeMetric > 0) %>%
      dplyr::select(-ChangeDirection)

    utils::write.csv(upreg, file = BiocGenerics::paste(filename, "_upregulated", ".csv", sep = ""))
  }
  # if direction is "down", export only genes with - change metric
  else if (stringi::stri_cmp_equiv(direction, "down", strength = 1)) {
    downreg <- tab %>%
      dplyr::filter(ChangeMetric < 0) %>%
      dplyr::select(-ChangeDirection)

    utils::write.csv(downreg, file = BiocGenerics::paste(filename, "_downregulated", ".csv", sep = ""))
  }
  # if no direction is given export all genes
  else {
    utils::write.csv(tab, file = BiocGenerics::paste(filename, ".csv", sep = ""))
  }
}
