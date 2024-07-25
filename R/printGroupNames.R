#' Print Group Names
#' Prints out group names. Pairwise comparison code must use group names to match exactly as returned here
#'
#' @param dseqobj dseq object; the results of running the dseq code
#'
#' @return prints strings of group names
#' @export
#'
#' @examples
#' print_group_names(dseq)
print_group_names <- function(dseqobj){
  names <- DESeq2::resultsNames(dseqobj)
  names <- names[-1]
  names <- stringr::str_remove_all(names, "Group_")
  names <- stringr::str_split(names, "_vs_")
  names <- purrr::flatten(names)
  names <- purrr::keep(names, !base::duplicated(names))
  names <- base::unlist(names)
  print(names)
}
