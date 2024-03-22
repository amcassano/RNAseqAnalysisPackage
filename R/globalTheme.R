#' Global Theme
#'
#' Sets the global GGplot theme
#'
#' @return ggplot theme
#' @export
#'
#' @examples
#' global_theme()
global_theme <- function() {
  return(
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      axis.line = ggplot2::element_line("black", 0.75),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 13, face = "italic"),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11.75),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9.5)
    )
  )
}
