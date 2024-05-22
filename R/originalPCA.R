#' Old PCA function
#' this is the original verison of the PCA plot function. it uses the plot PCA function from the DESeq2 package.
#' there are less customize options in this than the other PCA plot
#'
#' @param dseq_transform deseq transform object, rlog_norm not rlog_df
#' @param plot_aes the dataframe containing your plot aesthetics
#' @param plot_title string, plot title
#' @param label_samples boolean, should each dot be labeled
#'
#' @return pca plot
#' @export
#'
#' @examples
#' originalPCA(rlog_norm, plot_aes, "PCA")
originalPCA <- function(dseq_transform, plot_aes, plot_title = "Principal Component Analysis", label_samples = FALSE) {
  # save the results of the PCA plot as an object so that we can plot the data
  pca_data <-
    DESeq2::plotPCA(dseq_transform, intgroup = "Group", returnData = TRUE)

  # get the percent variation of the PCA data so that we can add that to the plot
  percent_variation <- round(100 * attr(pca_data, "percentVar"))
  labs <- plot_aes$labels
  outlines <- plot_aes$outlineID
  fills <- plot_aes$fillID
  shapes <- plot_aes$shapeID
  pca_plot <-
    ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = Group, shape = Group, fill = Group)) +
    ggplot2::geom_point(size = 3.5, stroke = 2) +
    ggplot2::scale_color_manual(labels = labs, values = outlines, name = "Group: ") +
    ggplot2::scale_fill_manual(labels = labs, values = fills, name = "Group: ") +
    ggplot2::scale_shape_manual(labels = labs, values = shapes, name = "Group: ") +
    ggplot2::xlab(paste0("PC1: ", percent_variation[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_variation[2], "% variance")) +
    ggplot2::labs(title = plot_title)

  # display plot with global theme +/- sample labels
  if (label_samples) {
    return(pca_plot + global_theme() + add_labels(pca_data))
  } else {
    return(pca_plot + global_theme())
  }
}
