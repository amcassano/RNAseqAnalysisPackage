#' PCA Analysis
#'
#' Illustrate PCA plot to show sample relationships
#'
#' @param dseq_transform DESeqTransform object
#' @param plot_title string, defaults to "Principal Component Analysis"
#' @param label_samples boolean, if true, labels for each sample will be added, defaults to FALSE
#' @param plot_aes dataframe, containing labels, colors, shapes, and fills for the data points being plotted
#' @param grouping_var string, defaults to Group, what is the grouping variable for this data
#'
#' @return PCA plot
#' @export
#'
#' @examples
#' pca_analysis(rlognorm, aesthetics)
pca_analysis <-
  function(dseq_transform, plot_aes, grouping_var = "Group", plot_title = "Principal Component Analysis", label_samples = FALSE){
    # save the results of the PCA plot as an object so that we can plot the data
    pca_data <-
      DESeq2::plotPCA(dseq_transform, intgroup = grouping_var, returnData = TRUE)

    # get the percent variation of the PCA data so that we can add that to the plot
    percent_variation <- round(100 * attr(pca_data, "percentVar"))

    labs <- plot_aes$labels
    outlines <- plot_aes$outlineID
    fills <- plot_aes$fillID
    shapes <- plot_aes$shapeID

    # create plot
    pca_plot <-
      ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = grouping_var, shape = grouping_var, fill = grouping_var)) +
      ggplot2::geom_point(size = 3.5, stroke = 2.25) +
      ggplot2::scale_color_manual(labels = labs, values = outlines, name = grouping_var) +
      ggplot2::scale_fill_manual(labels = labs, values = fills, name = grouping_var) +
      ggplot2::scale_shape_manual(labels = labs, values = shapes, name = grouping_var) +
      ggplot2::xlab(paste0("PC1: ", percent_variation[1], "% variance")) +
      ggplot2::ylab(paste0("PC2: ", percent_variation[2], "% variance")) +
      ggplot2::labs(title = plot_title)

    # display plot with global theme +/- sample labels
    if (label_samples) { pca_plot + global_theme() + add_labels(pca_data)}
    else {pca_plot + global_theme()}
  }
