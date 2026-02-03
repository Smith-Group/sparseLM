#' Add values to an existing sparse column in place
#'
#' @param x dgCMatrix to be assigned
#' @param i rows to be assigned in strictly increasing order
#' @param j single integer giving column to assign
#' @param value numeric vector (same length as \code{i}) to increment rows i of column j of x
#' @param validate logical scalar; if TRUE, validates that i is strictly increasing and within bounds
#'
#' @return value passed to x
#' @description
#' Adds \code{value} to the existing nonzero entries of column \code{j} at rows \code{i}.
#' Only existing nonzeros are updated; no new nonzeros are created.
#' \code{add_assign_col_inplace_unsafe} skips the input validation performed by
#' \code{add_assign_col_inplace}.
#' When \code{validate = TRUE}, the function checks that \code{i} is strictly
#' increasing, in range, and contains no NA; that \code{value} contains no NA/NaN;
#' and that the column pointer \code{p} is nondecreasing.
#'
#' @examples
#' if (requireNamespace("Matrix", quietly = TRUE)) {
#'   x <- Matrix::Matrix(0, nrow = 4, ncol = 2, sparse = TRUE)
#'   x[1, 1] <- 1
#'   x[3, 1] <- 2
#'   # value is the same length as i
#'   add_assign_col_inplace(x, i = c(1, 3), j = 1, value = c(10, -5))
#' }
#'
#' @export
add_assign_col_inplace <- function(x, i, j, value, validate = TRUE) {
	.Call("add_assign_col_inplace", x, i, j, value, validate)
}

#' @rdname add_assign_col_inplace
#' @export
add_assign_col_inplace_unsafe <- function(x, i, j, value) {
	.Call("add_assign_col_inplace", x, i, j, value, FALSE)
}
