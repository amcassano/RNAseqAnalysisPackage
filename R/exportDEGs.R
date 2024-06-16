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
#' export_genelist(ranked_tolvsnaive, "Naive-vs-tol")
#' export_genelist(ranked_tolvsnaive, "Naive-vs-tol", direction = "up")
#' export_genelist(ranked_tolvsnaive, "Naive-vs-tol", direction = "down")
export_genelist <- function(tab, filename,
                            direction = NULL) {
  # if no direction is given
  if (is.null(direction)) {utils::write.csv(tab, file = base::paste0(filename, ".csv"))}

  # if direction is "up", export only those genes with a + change metric
  else if (stringi::stri_cmp_equiv(direction, "up", strength = 1)) {
    upreg <- dplyr::filter(tab, ChangeMetric > 0)
    upreg <- dplyr::select(upreg, -ChangeDirection)
    utils::write.csv(upreg, file = base::paste0(filename, "_upregulated", ".csv"))
  }

  # if direction is "down", export only genes with - change metric
  else if (stringi::stri_cmp_equiv(direction, "down", strength = 1)) {
    downreg <- dplyr::filter(tab, ChangeMetric < 0)
    downreg <- dplyr::select(downreg, -ChangeDirection)
    utils::write.csv(downreg, file = base::paste0(filename, "_downregulated", ".csv"))
  }

  # if no direction is given or if it is anything other than up or down export all genes
  else {utils::write.csv(tab, file = base::paste0(filename, ".csv"))}
}
