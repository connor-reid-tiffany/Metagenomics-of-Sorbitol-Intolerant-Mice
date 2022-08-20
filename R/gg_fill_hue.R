#' Generate color palettes
#' @param n the number of colors to generate
#' @return a color for each n
#' @importFrom grDevices hcl

gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
