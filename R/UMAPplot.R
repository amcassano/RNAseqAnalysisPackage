#' UMAP Analysis
#'
#' plot umap for all samples, reducing the number of variables and showing relationships between samples
#'
#' @param dseq_transform DESeqTransform object
#' @param neighbors integer, minimum of 2
#' @param m_dist integer
#' @param cond_list list, conditions for each data point
#' @param plot_title string, plot title, defaults to "UMAP"
#' @param label_samples boolean, if sample labels should be added to plot, defaults to FALSE
#' @param plot_aes dataframe, containing labels, colors, shapes, and fills for the data points being plotted
#' @param batch boolean, defaults to false, set to true & add set seed for reproducible results
#'
#' @return UMAP plot
#' @export
#'
#' @examples
#' umap_analysis(rlognorm, 5, 10, plot_aes)

umap_analysis <- function(dseq_transform,  neighbors, m_dist, cond_list, plot_aes, plot_title = "UMAP", label_samples = FALSE, batch = FALSE) {

  # data must be transposed before being able to be used as input for UMAP
  transposed_normDEG <- BiocGenerics::as.data.frame(SummarizedExperiment::assay(dseq_transform))
  transposed_normDEG <- t(transposed_normDEG)

  # store the umap data in its own data frame
  zmap <- BiocGenerics::as.data.frame(uwot::umap(transposed_normDEG, n_neighbors = neighbors, min_dist = m_dist, batch = batch))
  zmap <- dplyr::rename(zmap, "UMAP1" = "V1")
  zmap <- dplyr::rename(zmap, "UMAP2" = "V2")

  sample_names <- BiocGenerics::rownames(zmap)

  # create a umap table combining the condition, sample name, and zmap data
  umap_data <- cbind(cond_list, sample_names, zmap)
  umap_data <- dplyr::rename(umap_data, "Group" = "cond_list")
  umap_data <- dplyr::rename(umap_data, "SampleName" = "sample_names")

  umap_data$Group <- factor(umap_data$Group, levels = plot_aes$labels)
  utils::head(umap_data)

  # plot the UMAP
  umap_plot <- ggplot2::ggplot(umap_data, ggplot2::aes(x = umap_data$UMAP1, y = umap_data$UMAP2,
                                                       color = umap_data$Group, shape = umap_data$Group, fill = umap_data$Group)) +
    global_theme() +
    ggplot2::geom_point(size = 3, stroke = 1.5) +
    ggplot2::scale_color_manual(labels = plot_aes$labels, values = plot_aes$outlineID, name = "Group:") +
    ggplot2::scale_fill_manual(labels = plot_aes$labels, values = plot_aes$fillID, name = "Group:") +
    ggplot2::scale_shape_manual(labels = plot_aes$labels, values = plot_aes$shapeID, name = "Group:") +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    ggplot2::labs(title = plot_title)

  # display plot
  if (label_samples) {return(umap_plot + add_labels(umap_data))}
  else {return(umap_plot)}
}
