#' Summarize curvature object.
#'
#' @import tibble
#' @export
#'
summary.curvature = function(res) {
  cat("sectional curvature stats:\n")
  tibble(mean = mean(res),
         median = median(res),
         sd = sd(res),
         min = min(res),
         max = max(res))
}
