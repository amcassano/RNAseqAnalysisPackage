#' Rename Columns
#'
#' @description  rename columns, this allows for consistent naming and easier access to data
#'
#' @param colstr string, column to be renamed
#'
#' @return string, edited column name
#' @export

rename_columns <- function(colstr){

  colstr <- BiocGenerics::paste(colstr)

  if (stringi::stri_cmp_equiv(colstr, "ensembl_gene_id", strength = 1)) {
    return("GeneID") }

  if(stringi::stri_cmp_equiv(colstr, "mgi_symbol", strength = 1)){
    return("MGI_Symbol") }

  else if(stringi::stri_cmp_equiv(colstr, "mgi_description", strength = 1)){
    return("MGI_Desc") }

  else if(stringi::stri_cmp_equiv(colstr, "gene_biotype", strength = 1)){
    return("GeneType") }

  else if(stringi::stri_cmp_equiv(colstr, "entrezgene_id", strength = 1)){
    return("EntrezID") }

  else if(stringi::stri_cmp_equiv(colstr, "go_id", strength = 1)){
    return("GO_ID") }

  else if(stringi::stri_cmp_equiv(colstr, "log2FoldChange", strength = 1)){
    return("Log2FoldChange")
  }

  else if(stringi::stri_cmp_equiv(colstr, "lfcSE", strength = 1)){
    return("Log2FC_SE")
  }

  else if(stringi::stri_cmp_equiv(colstr, "pvalue", strength = 1)){
    return("P_Value")
  }

  else if(stringi::stri_cmp_equiv(colstr, "padj", strength = 1)){
    return("Adj_P_Value")
  }

  else{return(colstr)}
}
