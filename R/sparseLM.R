#' @useDynLib sparseLM, .registration = TRUE
#' @import Matrix
NULL

#' Sparse Levenberg-Marquardt Algorithm Options
#'
#' @param mu scale factor for initial mu
#' @param epsilon1 stopping threshold for ||J^T e||_inf
#' @param epsilon2 stopping threshold for ||dp||_2
#' @param epsilon3 stopping threshold for ||e||_2
#' @param delta step used in difference approximation to the Jacobian
#'
#' @return numeric vector of length 6 with the above options and spsolver=1 (SuiteSparse CHOLMOD)
#'
#' @export
sparselm.opts <- function(mu=1e-03, epsilon1=1e-12, epsilon2=1e-12, epsilon3=1e-12, delta=1E-06) {

	c(mu=mu, epsilon1=epsilon1, epsilon2=epsilon2, epsilon3=epsilon3, delta=delta, spsolver=1)
}

#' Nonlinear Least Squares Fit with Sparse Levenberg-Marquardt Algorithm
#'
#' @param p initial parameter estimates length nvars
#' @param x measurement vector length nobs
#' @param func functional relation describing measurements given a parameter vector p, returning vector length nobs
#' @param fjac function to supply the nonzero pattern of the sparse Jacobian of func and optionally evaluate it at p, returning nobs by nvars dgCMatrix
#' @param Jnnz number of nonzeros for the Jacobian J
#' @param JtJnnz number of nonzeros for the product J^t*J, -1 if unknown
#' @param nconvars number of constrained variables (currently reserved for future use)
#' @param itmax maximum number of iterations
#' @param opts minim. options \code{mu, epsilon1, epsilon2, epsilon3, delta, spsolver}, not currently implemented!
#' @param dif logical indicating whether to use finite differences
#' @param ... additional arguments passed to func and fjac
#'
#' @return list with four elements: par, niter, info, and term
#'
#' @export
sparselm <- function(p, x, func, fjac, Jnnz, JtJnnz=-1, nconvars=0, itmax=100, opts=sparselm.opts(), dif=FALSE, ...) {

	func1 <- function(p) func(p, ...)
	fjac1 <- function(p) fjac(p, ...)
	
	out <- .Call("C_sparselm", func1, fjac1, p, x, as.integer(nconvars), as.integer(Jnnz), as.integer(JtJnnz), as.integer(itmax), as.numeric(opts), as.logical(dif), new.env(), PACKAGE = "sparseLM")

	names(out) <- c("par", "niter", "info")
	names(out[["info"]]) <- c("rssinit", "rss", "||J^T e||_inf", "||dp||_2", "mu/max[J^T J]_ii", "niter", "term", "nfunc", "nfjac", "nsys")
	out[["term"]] <- c(
		"stopped by small gradient J^T e",
		"stopped by small dp",
		"stopped by itmax",
		"singular matrix. Restart from current p with increased mu",
		"too many failed attempts to increase damping. Restart with increased mu",
		"stopped by small ||e||_2",
		"stopped by invalid (i.e. NaN or Inf) func values. User error"
	)[out[["info"]][["term"]]]
	
	out
}
