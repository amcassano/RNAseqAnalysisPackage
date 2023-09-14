#' Get Shared Genes
#'
#' @param suffix_list list of strings, one for each data frame in the
#' @param d1 data frame, contains at least MGI symbol and description, log 2 fold change, adjusted p value, and change direction
#' @param d2 data frame, contains at least MGI symbol and description, log 2 fold change, adjusted p value, and change direction
#' @param ... data frames, contains at least MGI symbol and description, log 2 fold change, adjusted p value, and change direction
#' @param samedir boolean, defaults to true, dictates if the changes need to be in the same direction
#'
#' @return data frame, contains all the genes significantly up/down across given data frames including fc, p val,
#' @export
#'
#' @examples
#' get_shared_genes(c(".tol", ".rej"), naive_vs_tol, naive_vs_rej)
#' get_shared_genes(c(".tol", ".rej", ".ex-tol"), naive_vs_tol, naive_vs_rej, naive_vs_extol)
#' get_shared_genes(c(".tol", ".rej", ".ex-tol"), naive_vs_tol, naive_vs_rej, naive_vs_extol, samedir = FALSE)
get_shared_genes <- function(suffix_list, dataframe1, dataframe2, ..., samedir = TRUE){

  dfsToAdd <- list(...)
  dfsToAdd <- c(dataframe2, dfsToAdd)
  dfsLeft <- length(dfsToAdd)

  interiorLoop <- function(total_dfs, remaining_dfs, dfs_to_add, suffixes = suffix_list, shared = dataframe1){
    di <- total_dfs - remaining_dfs + 1
    si <- di + 1
    suffix1 <- suffixes[1]
    if (total_dfs == 1) {
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c(suffix1, suffix2))
      mgi1 <- BiocGenerics::paste("MGI_Symbol", suffix1, sep = "")
      desc1 <- BiocGenerics::paste("MGI_Desc", suffix1, sep = "")
      dir1 <- BiocGenerics::paste("ChangeDirection", suffix1, sep = "")
      fc1 <- BiocGenerics::paste("Log2FoldChange", suffix1, sep = "")
      pv1 <- BiocGenerics::paste("Adj_P_Value", suffix1, sep = "")
      dir2 <- BiocGenerics::paste("ChangeDirection", suffix2, sep = "")
      fc2 <- BiocGenerics::paste("Log2FoldChange", suffix2, sep = "")
      pv2 <- BiocGenerics::paste("Adj_P_Value", suffix2, sep = "")
      shared <- shared %>% dplyr::select("MGI_Symbol" = (!!mgi1), "MGI_Desc" = (!!desc1), "ChangeDirection" = (!!dir1),
                                         (!!fc1), (!!pv1), (!!dir2), (!!fc2), (!!pv2))

      if(samedir) {shared <- shared %>% dplyr::filter(ChangeDirection == (!!dir2))}
      shared <- shared %>% dplyr::select(-(!!dir2))

      return(shared)
    }
    else if(total_dfs == remaining_dfs){
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c(suffix1, suffix2))
      mgi1 <- BiocGenerics::paste("MGI_Symbol", suffix1, sep = "")
      desc1 <- BiocGenerics::paste("MGI_Desc", suffix1, sep = "")
      dir1 <- BiocGenerics::paste("ChangeDirection", suffix1, sep = "")
      fc1 <- BiocGenerics::paste("Log2FoldChange", suffix1, sep = "")
      pv1 <- BiocGenerics::paste("Adj_P_Value", suffix1, sep = "")
      dir2 <- BiocGenerics::paste("ChangeDirection", suffix2, sep = "")
      fc2 <- BiocGenerics::paste("Log2FoldChange", suffix2, sep = "")
      pv2 <- BiocGenerics::paste("Adj_P_Value", suffix2, sep = "")
      shared <- shared %>% dplyr::select("MGI_Symbol" = (!!mgi1), "MGI_Desc" = (!!desc1), "ChangeDirection" = (!!dir1),
                                         (!!fc1), (!!pv1), (!!dir2), (!!fc2), (!!pv2))

      if(samedir) {shared <- shared %>% dplyr::filter(ChangeDirection == (!!dir2))}
      shared <- shared %>% dplyr::select(-(!!dir2))

      remaining_dfs <- remaining_dfs - 1
      interiorLoop(total_dfs, remaining_dfs, dfs_to_add, suffixes, shared)
    }
    else if (remaining_dfs == 1){
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      change_dir_comp <- BiocGenerics::paste("ChangeDirection", suffix1, sep = "")
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c("", suffix2))
      mgi2 <- BiocGenerics::paste("MGI_Symbol", suffix2, sep = "")
      desc2 <- BiocGenerics::paste("MGI_Desc", suffix2, sep = "")
      dir2 <- BiocGenerics::paste("ChangeDirection", suffix2, sep = "")
      shared <- shared %>%  dplyr::select(-(!!mgi2)) %>% dplyr::select(-(!!desc2))
      if(samedir) {shared <- shared %>% dplyr::filter(ChangeDirection == (!!dir2))}
      shared <- shared %>% dplyr::select(-(!!dir2))
      return(shared)
    }
    else{
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c("", suffix2))
      mgi2 <- BiocGenerics::paste("MGI_Symbol", suffix2, sep = "")
      desc2 <- BiocGenerics::paste("MGI_Desc", suffix2, sep = "")
      dir2 <- BiocGenerics::paste("ChangeDirection", suffix2, sep = "")
      shared <- shared %>%  dplyr::select(-(!!mgi2)) %>% dplyr::select(-(!!desc2))
      if(samedir) {shared <- shared %>% dplyr::filter(ChangeDirection == (!!dir2))}
      shared <- shared %>% dplyr::select(-(!!dir2))
      remaining_dfs <- remaining_dfs - 1
      interiorLoop(total_dfs, remaining_dfs, dfs_to_add, suffixes, shared)
    }
  }

  sharedgenes <- interiorLoop(dfsLeft, dfsLeft, dfsToAdd, suffix_list)

  return(sharedgenes)
}
