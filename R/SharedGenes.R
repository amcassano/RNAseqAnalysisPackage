#' Get Shared Genes
#'
#' @param suffix_list list of strings, one for each data frame inputted, in the same order
#' @param dataframe1 data frame, contains at least MGI symbol and description, log 2 fold change, adjusted p value, and change direction
#' @param dataframe2 data frame, contains at least MGI symbol and description, log 2 fold change, adjusted p value, and change direction
#' @param ... data frames, contains at least MGI symbol and description, log 2 fold change, adjusted p value, and change direction
#'
#' @return data frame, contains all the genes significantly up/down in the same direction across given data frames including fc, p val,
#' @export
#'
#' @examples
#' get_shared_genes(c(".naive_tol", ".naive_rej"), sig_NaiveTol, sig_NaiveRej)
#' get_shared_genes(c(".naive_tol", ".naive_rej", ".naive_ex-tol"), sig_NaiveTol, sig_NaiveRej, sig_NaiveExTol)
#' get_shared_genes(c(".naive_tol", ".naive_rej", ".naive_ex-tol", ".naive_memory), sig_NaiveTol, sig_NaiveRej, sig_NaiveExTol, sig_NaiveMemory)
get_shared_genes <- function(suffix_list, dataframe1, dataframe2, ...) {
  # set up the intial variables for the interior loop
  dfsToAdd <- list(...)
  dataframe2 <- list(dataframe2)
  dfsToAdd <- c(dataframe2, dfsToAdd)
  dfsLeft <- as.numeric(length(dfsToAdd))
  totaldfs <- as.numeric(length(dfsToAdd))

  #interior loop to iterate through all of the data frames and find significantly changed shared genes
  interiorLoop <- function(total_dfs, remaining_dfs, dfs_to_add, suffixes, shared) {
    di <- total_dfs - remaining_dfs + 1
    si <- di + 1
    suffix1 <- suffixes[1]
    if (total_dfs == 1) {
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c(suffix1, suffix2))
      mgi1 <- base::paste0("MGI_Symbol", suffix1)
      desc1 <- base::paste0("MGI_Desc", suffix1)
      dir1 <- base::paste0("ChangeDirection", suffix1)
      fc1 <- base::paste0("Log2FoldChange", suffix1)
      pv1 <- base::paste0("Adj_P_Value", suffix1)
      dir2 <- base::paste0("ChangeDirection", suffix2)
      fc2 <- base::paste0("Log2FoldChange", suffix2)
      pv2 <- base::paste0("Adj_P_Value", suffix2)
      shared <- dplyr::select(shared, "MGI_Symbol" = (!!mgi1), "MGI_Desc" = (!!desc1), "ChangeDir" = (!!dir1),
                              (!!fc1), (!!pv1), "ChangeDir2" = (!!dir2), (!!fc2), (!!pv2))
      shared$sameDir <- FALSE
      shared$sameDir[shared$ChangeDir == "Up" & shared$ChangeDir2 == "Up"] <- TRUE
      shared$sameDir[shared$ChangeDir == "Down" & shared$ChangeDir2 == "Down"] <- TRUE
      shared <- dplyr::filter(shared, sameDir == TRUE)
      shared <- dplyr::select(shared, -ChangeDir2)
      return(shared)
    }
    else if (total_dfs == remaining_dfs) {
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c(suffix1, suffix2))
      mgi1 <- base::paste0("MGI_Symbol", suffix1)
      desc1 <- base::paste0("MGI_Desc", suffix1)
      dir1 <- base::paste0("ChangeDirection", suffix1)
      fc1 <- base::paste0("Log2FoldChange", suffix1)
      pv1 <- base::paste0("Adj_P_Value", suffix1)
      dir2 <- base::paste0("ChangeDirection", suffix2)
      fc2 <- base::paste0("Log2FoldChange", suffix2)
      pv2 <- base::paste0("Adj_P_Value", suffix2)
      shared <- dplyr::select(shared, "MGI_Symbol" = (!!mgi1), "MGI_Desc" = (!!desc1), "ChangeDir" = (!!dir1),
                              (!!fc1), (!!pv1), "ChangeDir2" = (!!dir2), (!!fc2), (!!pv2))
      shared$sameDir <- FALSE
      shared$sameDir[shared$ChangeDir == "Up" & shared$ChangeDir2 == "Up"] <- TRUE
      shared$sameDir[shared$ChangeDir == "Down" & shared$ChangeDir2 == "Down"] <- TRUE
      shared <- dplyr::filter(shared, sameDir == TRUE)
      shared <- dplyr::select(shared, -ChangeDir2)
      remaining_dfs <- remaining_dfs - 1
      interiorLoop(total_dfs, remaining_dfs, dfs_to_add, suffixes, shared)
    }
    else if (remaining_dfs == 1) {
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c("", suffix2))
      mgi2 <- base::paste0("MGI_Symbol", suffix2)
      desc2 <- base::paste0("MGI_Desc", suffix2)
      dir2 <- base::paste0("ChangeDirection", suffix2)
      shared <- dplyr::select(shared, -(!!mgi2))
      shared <- dplyr::select(shared, -(!!desc2))
      names(shared)[names(shared) == "Adj_P_Value"] <- base::paste0("Adj_P_Value", suffix2)
      names(shared)[names(shared) == "Log2FoldChange"] <- base::paste0("Log2FoldChange", suffix2)
      shared$sameDir <- FALSE
      shared$sameDir[shared$ChangeDir == "Up" & shared$ChangeDirection == "Up"] <- TRUE
      shared$sameDir[shared$ChangeDir == "Down" & shared$ChangeDirection == "Down"] <- TRUE
      shared <- dplyr::filter(shared, sameDir == TRUE)
      shared <- dplyr::select(shared, -c(ChangeDirection, ChangeMetric, sameDir))
      return(shared)
    }
    else {
      d1 <- shared
      d2 <- dfs_to_add[[di]]
      suffix2 <- suffixes[si]
      shared <- dplyr::inner_join(d1, d2, by = "MGI_Symbol", keep = TRUE, suffix = c("", suffix2))
      mgi2 <- base::paste0("MGI_Symbol", suffix2)
      desc2 <- base::paste0("MGI_Desc", suffix2)
      dir2 <- base::paste0("ChangeDirection", suffix2)
      shared <- dplyr::select(shared, -(!!mgi2))
      shared <- dplyr::select(shared, -(!!desc2))
      names(shared)[names(shared) == "Adj_P_Value"] <- base::paste0("Adj_P_Value", suffix2)
      names(shared)[names(shared) == "Log2FoldChange"] <- base::paste0("Log2FoldChange", suffix2)
      shared$sameDir <- FALSE
      shared$sameDir[shared$ChangeDir == "Up" & shared$ChangeDirection == "Up"] <- TRUE
      shared$sameDir[shared$ChangeDir == "Down" & shared$ChangeDirection == "Down"] <- TRUE
      shared <- dplyr::filter(shared, sameDir == TRUE)
      shared <- dplyr::select(shared, -c(ChangeDirection, ChangeMetric))
      remaining_dfs <- remaining_dfs - 1
      interiorLoop(total_dfs, remaining_dfs, dfs_to_add, suffixes, shared)
    }
  }

  sharedgenes <- interiorLoop(total_dfs = totaldfs,
                              remaining_dfs = dfsLeft,
                              dfs_to_add = dfsToAdd,
                              suffixes = suffix_list,
                              shared = dataframe1)
  return(sharedgenes)
}
