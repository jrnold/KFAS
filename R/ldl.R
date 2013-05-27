#' LDL Decomposition of a Matrix
#'
#' Function \code{ldl} computes the LDL decomposition of a symmetric matrix.
#'
#'
#' @export
#' @param x Symmetrix matrix.
#' @param tol Tolerance parameter for LDL decomposition, determines which
#' diagonal values are counted as zero.
ldl <- function(x, tol = .Machine$double.eps^0.5) {
    if (!isSymmetric(x, tol = tol)) 
        stop("Matrix is not symmetric!")
    .Fortran("ldl", PACKAGE = "KFAS", x, as.integer(dim(x)[1]), tol = tol, as.integer(0))[[1]]
} 
