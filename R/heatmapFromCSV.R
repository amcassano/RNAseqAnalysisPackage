#' Heatmap from CSV
#'
#' @param plottitle string, title of the plot
#' @param fname string, full filename for the csv listing genes of interest, including path and .csv
#' @param met data frame, meta data information
#' @param counts_df data frame, normalized counts
#' @param clusterRows boolean, defaults to TRUE, should the rows cluster by similarity
#' @param gaps list of numbers, defaults to empty list, dictates where there would be gaps if desired between conditions
#' @param height number, defaults to 10, dictates how tall the cells in the heatmap will be
#' @param annocolors list of strings, list of colors for use in labeling
#' @param colorscale list of 3 colors sets the color scale, defaults to red (hi) white blue (low)
#' @param scaling boolean, defaults to TRUE - if rows should be scaled or not
#' @param displNumbers boolean, defaults to FALSE - if the value of each cell should be displayed within#'
#' @param savetofile filepath where to save plot Filetype decided by extension. png, pdf, tiff, bmp, jpeg.
#'
#' @return heatmap plot
#' @export
#'
#' @examples
#' heatmapFromCSV("OxPhos", "oxphos.csv", metadata, rlog_df, c("blue", "red"))
heatmapFromCSV <- function(plottitle,
                           fname,
                           met,
                           counts_df,
                           annocolors,
                           clusterRows = TRUE,
                           gaps = c(),
                           height = 10,
                           colorscale = c("blue4", "white", "red3"),
                           scaling = TRUE,
                           displNumbers = FALSE,
                           savetofile = NA) {
    # read in the list of genes to include in the heatmap
    genelist <- utils::read.csv(fname, header = FALSE)
    colnames(genelist) <- "MGI_Symbol"
    genelist <- dplyr::distinct(genelist, MGI_Symbol)

    # create data frame with counts for corresponding genes
    heatmap_data <- dplyr::select(counts_df, -MGI_Desc)
    heatmap_data <- dplyr::filter(heatmap_data, MGI_Symbol %in% genelist$MGI_Symbol)
    heatmap_data <- tibble::rownames_to_column(heatmap_data, var = "GeneID")

    create_heatmap(
        title = plottitle,
        genes_and_counts = heatmap_data,
        met = met,
        annocolors = annocolors,
        clusterRows = clusterRows,
        gaps = gaps,
        height = height,
        colorscale = colorscale,
        scaling = scaling,
        displNumbers = displNumbers,
        savetofile = savetofile
    )
}
