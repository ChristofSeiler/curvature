#' Plot curvature object.
#'
#' @import tibble
#' @import ggplot2
#' @export
#'
plot.curvature = function(res) {
  ggplot(tibble(sec = c(res)),aes(sec)) +
    geom_histogram(bins = 30)
}
