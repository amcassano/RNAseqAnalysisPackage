#' Create Genemap
#'
#' create a reference dataframe containing ensembl gene IDs and corresponding MGI symbols for use in data annotation
#'
#' @param allcounts data frame that at a minimum includes ensembl Gene ID for all of the genes in the data set in a column called GeneID
#' @param bio_mart, string, defaults to ensembl, dictates what is used to build the genemap
#' @param map_genome string, defaults to "mmusculus_gene_ensembl", dictates which genome to use to build the gene map
#' @param GMattributes list of strings, attributes to include in genemap, defaults to ensmbl ID, mgi symbol, mgi description and gene biotype
#'
#' @return dataframe containing all ensmbl IDs and corresponding MGI symbols and descriptions and gene biotypes
#' @export
#'
#' @examples
#' createGenemap(allcounts = allgenes)
#' createGenemap(allgenes, GMattributes = c("ensembl_gene_id", "mgi_symbol"))
#'
createGenemap <- function(allcounts, bio_mart = "ensembl", map_genome = "mmusculus_gene_ensembl",
                          GMattributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description", "gene_biotype")) {
  # get list of all gene ids from the list of all counts
  all_genes <- as.data.frame(allcounts$GeneID)

  # create biomart object
  ensembl <-
    biomaRt::useMart(biomart = bio_mart, dataset = map_genome)

  genemap <-
    biomaRt::getBM(attributes = GMattributes, filters = "ensembl_gene_id", values = all_genes, mart = ensembl)

  genemap <-
    dplyr::distinct(genemap, genemap$ensembl_gene_id, .keep_all = TRUE)

  return(genemap)
}
