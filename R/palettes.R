#' Palettes
#'
#' just sets the color palettes and prints the options
#'
#' @return prints the options for colors
#' @export
#'
#' @examples
#' getPalettes()
getPalettes <- function(){
  colorpalette <-
    data.frame(
      colorID = c("#d50000", "#0000ca", "#2c902c", "#a300de", "#ffef0d", "#ed9a00", "#ff2fd2", "#79491c", "#000000", "#FFFFFF", "#9d9d9d"),
      colorName = c("red", "blue", "green", "purple", "yellow", "orange", "pink", "brown","black", "white", "grey"))

  # dataframe, shape choices for plots
  shape_palette <- data.frame(
    shapeID = c(21, 10, 23, 9, 22, 7,24, 25, 4, 3, 8),
    shapeDescription = c("circle", "crossed circle", "diamond", "crossed diamond", "square", "crossed square",
                         "triangle up", "triangle down", "x", "plus", "asterik"))

  print("Color Options: \n")
  print(colorpalette$colorName)
  print("Shape Options: \n")
  print(shape_palette$shapeDescription)
}
