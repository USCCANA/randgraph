
#' Density of a graph
#' @param x A graph.
#' @export
gdensity <- function(x) UseMethod("gdensity")

#' @export
#' @rdname gdensity
gdensity.igraph <- function(x) {
  n <- igraph::vcount(x)
  e <- igraph::ecount(x)

  e/n/n
}

