#' Consolidate Gene List CSVs
#'
#' @param map dataframe, genemap
#' @param ... strings, file names of CSVs to consolidate
#'
#' @return data frame, all unique genes from each data frame combined
#' @export
#'
#' @examples
#' consolidate_gene_list(genemap, "oxphos.csv", "activation.csv", "proliferation.csv")
consolidate_gene_list <- function(map, ...) {
  #set up the mulitple inputs for the interior loop
  filenames <- list(...)
  num_of_lists <- length(filenames)

  interior_loop <- function(total_lists, lists_left, filename_list, dataframe_so_far = NULL) {
    list_index <- total_lists - lists_left + 1
    fname <- paste(filename_list[[list_index]])
    if (lists_left == 1) {
      temp_df <- utils::read.csv(fname, header = FALSE)
      colnames(temp_df) <- "MGI_Symbol"
      dataframe_so_far <- dplyr::bind_rows(temp_df, dataframe_so_far)
      dataframe_so_far <- dplyr::distinct(dataframe_so_far, dataframe_so_far$MGI_Symbol, .keep_all = TRUE)
      return(dataframe_so_far)
    }
    else if (lists_left == total_lists) {
      dataframe_so_far <- utils::read.csv(fname, header = FALSE)
      colnames(dataframe_so_far) <- "MGI_Symbol"
      lists_left <- lists_left - 1
      interior_loop(total_lists, lists_left, filename_list, dataframe_so_far)
    }
    else {
      temp_df <- utils::read.csv(fname, header = FALSE)
      colnames(temp_df) <- "MGI_Symbol"
      dataframe_so_far <- dplyr::bind_rows(temp_df, dataframe_so_far)
      dataframe_so_far <- dplyr::distinct(dataframe_so_far, dataframe_so_far$MGI_Symbol, .keep_all = TRUE)
      lists_left <- lists_left - 1
      interior_loop(total_lists, lists_left, filename_list, dataframe_so_far)
    }
  }
  genelist_temp <- interior_loop(num_of_lists, num_of_lists, filenames)

  combined_genelist <- dplyr::inner_join(x = genelist_temp, y = map, by = c("MGI_Symbol" = "mgi_symbol"))
  combined_genelist <- dplyr::select(combined_genelist, c("GeneID" = combined_genelist$ensembl_gene_id, combined_genelist$MGI_Symbol))

  return(combined_genelist)
}
