#' Get subset of adjacency matrix that is connected to one variable
#'
#' @param adjMat Full adjacency matrix
#' @param targetIndex index of variable of interest
#' @param nDeg degrees of separation to include
#' @return an adjacency matrix only including the variables within \code{nDeg}
#' of separation from the variable of interest

graph.minimal.n <- function(adjMat, targetIndex, nDeg = 1){

  edges.keep <- c(targetIndex)
  edges.deg <- list()
  for(i in 1:nDeg){
    edges.keep.old <- edges.keep
    edges.keep <- unique(c(edges.keep, which(rowSums(adjMat[, edges.keep, drop=FALSE]) > 0)))
    edges.deg[[i]] <- setdiff(edges.keep, edges.keep.old)
  }
  return(adjMat[edges.keep, edges.keep])
}
