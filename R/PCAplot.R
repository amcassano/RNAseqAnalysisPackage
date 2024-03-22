#' PCA Analysis
#'
#' Illustrate PCA plot to show sample relationships
#'
#' @param norm_df data frame with the dseq_transform results (i.e. rlog_df)
#' @param plot_title string, defaults to "Principal Component Analysis"
#' @param label_samples boolean, if true, labels for each sample will be added, defaults to FALSE
#' @param plot_aes dataframe, containing labels, colors, shapes, and fills for the data points being plotted
#' @param metadat metadata dataframe
#' @param varianceRemove defaults to 0.05 number that corresponds to the percent
#' @param plot_subtitle subtitle for plot, defaults to null
#' @param plot_caption plot caption, defaults to null
#' @param draw_polygon boolean; should the groups be circled
#' @param draw_ellipses boolean; should stat ellipses be drawn
#' @param point_size size of points
#' @param point_stroke stroke thickness of points
#' @param label_segment_size width of lines connecting dots and labeles
#' @param label_size size of the labels
#' @param label_point_padding point padding for geom-repel
#' @param label_box_padding box padding for geom repel
#' @param label_force_pull for geom repel
#' @param label_force_push for geom repel
#' @param max_overlaps for label overlaps
#' @param circle_transparency alpha for either ellipse or polygon
#' @param circle_border_width border width for either ellipse or polygon
#' @param ellipse_type for the ellipse either "t" or "norm"
#' @param pcs a list of 2 strings that are the principal components to use (c("PC1", "PC2"))
#'
#' @return PCA plot
#' @export
#'
#' @examples
#' pca_analysis(rlog_df, plot_aes, meta)
pca_analysis <-
  function(norm_df,
           plot_aes,
           metadat,
           varianceRemove = 0.05,
           plot_title = "Principal Component Analysis",
           plot_subtitle = NULL,
           plot_caption = NULL,
           draw_polygon = FALSE,
           draw_ellipses = TRUE,
           label_samples = FALSE,
           point_size = 4,
           point_stroke = 2.5,
           label_segment_size = 0.35,
           label_size = 3.5,
           label_point_padding = 0.25,
           label_box_padding = 0.75,
           label_force_pull = 1.2,
           label_force_push = 1.35,
           max_overlaps = 30,
           circle_transparency = 4/15,
           circle_border_width = 3,
           ellipse_type = "t",
           pcs = c("PC1", "PC2")) {
    pca_obj <- get_pca_obj(norm_df, metadat, varianceRemove)

    labs <- plot_aes$labels
    outlines <- plot_aes$outlineID
    fills <- plot_aes$fillID
    shapes <- plot_aes$shapeID


    if (draw_ellipses) {
      xlim <- c(min(pca_obj$rotated[, pcs[1]]) - abs((min(pca_obj$rotated[, pcs[1]]) / 100) * 35),
                max(pca_obj$rotated[, pcs[1]]) + abs((min(pca_obj$rotated[, pcs[1]]) / 100) * 35))

      ylim <- c(min(pca_obj$rotated[, pcs[2]]) - abs((min(pca_obj$rotated[, pcs[2]]) / 100) * 35),
                max(pca_obj$rotated[, pcs[2]]) + abs((min(pca_obj$rotated[, pcs[2]]) / 100) * 35))
    }
    else{
      xlim <- c(min(pca_obj$rotated[, pcs[1]]) - abs((min(pca_obj$rotated[, pcs[1]]) / 100) * 10),
                max(pca_obj$rotated[, pcs[1]]) + abs((min(pca_obj$rotated[, pcs[1]]) / 100) * 10))

      ylim <- c(min(pca_obj$rotated[, pcs[2]]) - abs((min(pca_obj$rotated[, pcs[2]]) / 100) * 10),
                max(pca_obj$rotated[, pcs[2]]) + abs((min(pca_obj$rotated[, pcs[2]]) / 100) * 10))
    }

    #make the plot
    labelsFunc <- xidx <- yidx <- NULL
    plot_obj <- NULL
    plot_obj$x <- pca_obj$rotated[, pcs[1]]
    plot_obj$y <- pca_obj$rotated[, pcs[2]]
    samplelabel <- rownames(pca_obj$metadata)
    if (label_samples) {
      samplelabel <- rownames(pca_obj$metadata)
      plot_obj$samplelabel <- samplelabel
    }
    else{
      samplelabel <- NULL
    }
    plot_obj <- as.data.frame(plot_obj, stringsAsFactors = FALSE)
    plot_obj$Group <- pca_obj$metadata[, "Group"]

    pca_plot <-
      ggplot2::gpplot(plot_obj, ggplot2::aes(x = x, y = y, color = Group, shape = Group, fill = Group)) +
      ggplot2::geom_point(size = point_size, stroke = point_stroke) +
      ggplot2::scale_color_manual(labels = labs, values = outlines,name = "Group: ") +
      ggplot2::scale_fill_manual(labels = labs, values = fills, name = "Group: ") +
      ggplot2::scale_shape_manual(labels = labs, values = shapes,name = "Group: ") +
      ggplot2::xlab(paste0(pcs[1], ": ", round(pca_obj$variance[pcs[1]]), "% variance")) +
      ggplot2::ylab(paste0(pcs[2], ": ", round(pca_obj$variance[pcs[2]]), "% variance")) +
      ggplot2::labs(title = plot_title, subtitle = plot_subtitle, caption = plot_caption) +
      ggplot2::xlim(xlim[1], xlim[2]) +
      ggplot2::ylim(ylim[1], ylim[2])

    pca_plot <-
      pca_plot +
      global_theme()

    if(label_samples){
      pca_plot <-
        pca_plot +
        ggrepel::geom_text_repel(
          data = plot_obj,
          ggplot2::aes(label = samplelabel),
          segment.size = label_segment_size,
          min.segment.length = 0,
          size = label_size,
          point.padding = label_point_padding, box.padding = label_box_padding,
          max.overlaps = max_overlaps,
          show.legend = FALSE,
          force = label_force_push, force_pull = label_force_pull,
          na.rm = TRUE
        )
    }

    if (draw_polygon) {
      pca_plot <-
        pca_plot +
        ggalt::geom_encircle(
          ggplot2::aes(group = Group, fill = Group, color = Group),
          alpha = circle_transparency,
          size = circle_border_width,
          show.legend = FALSE,
          na.rm = TRUE
        )
    }

    if (draw_ellipses){
      pca_plot <-
        pca_plot +
        ggplot2::stat_ellipse(
          ggplot2::aes(group = Group, fill = Group, color = Group),
          geom = "polygon",
          type = ellipse_type,
          size = circle_border_width,
          alpha = circle_transparency,
          show.legend = FALSE,
          na.rm = TRUE
        )
    }

    return(pca_plot)
  }
