#' Set Aesthetics
#'
#' @param label list of strings
#' @param shapes list of strings
#' @param fills list of strings
#' @param outlines list of strings
#'
#' @return data frame, aesthetics for use in plotting
#' @export
#'
#' @examples
#' set_aes(
#'   condition_labels,
#'   c("circle", "square", "triangle up", "x"),
#'   c("blue", "green", "red", "orange"),
#'   c("blue", "green", "red", "orange"))
#'
set_aes <- function(label, shapes, outlines, fills){
  # data frame, color choices
  colorpalette <-
    data.frame(
      colorID = c("#d50000", "#0000ca", "#2c902c", "#a300de", "#ffef0d", "#ed9a00", "#ff2fd2", "#79491c", "#000000", "#FFFFFF", "#9d9d9d"),
      colorName = c("red", "blue", "green", "purple", "yellow", "orange", "pink", "brown","black", "white", "grey"))

  # dataframe, shape choices for plots
  shape_palette <-
    data.frame(
      shapeID = c(21, 10, 23, 9, 22, 7,24, 25, 4, 3, 8),
      shapeDescription = c("circle", "crossed circle", "diamond", "crossed diamond",
                           "square", "crossed square", "triangle up", "triangle down", "x", "plus", "asterik"))

  plot_aes <- data.frame(labels = label, shapeDescription = shapes, fillName = fills, outlineName = outlines)
  plot_aes <- dplyr::left_join(plot_aes, shape_palette, by = "shapeDescription")
  plot_aes <- dplyr::left_join(plot_aes, colorpalette, by = c("fillName" = "colorName"))
  plot_aes <- dplyr::rename(plot_aes, fillID = colorID)
  plot_aes <- dplyr::left_join(plot_aes, colorpalette, by = c("outlineName" = "colorName"))
  plot_aes <- dplyr::rename(plot_aes, outlineID = colorID)
  plot_aes <- dplyr::select(plot_aes, labels, shapeDescription, shapeID, outlineName, outlineID, fillName, fillID)
  print(plot_aes)
  return(plot_aes)


}
