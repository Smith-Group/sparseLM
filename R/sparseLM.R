#' @useDynLib sparseLM, .registration = TRUE
#' @import Matrix
#' @importFrom methods as
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
#' @examples
#' set.seed(1)
#' t <- seq(0, 10, length.out = 80)
#' g <- function(A, mu, s, x) A * exp(-0.5 * ((x - mu) / s)^2)
#' y <- g(3, 3, 0.7, t) + g(2, 7, 1.0, t) + rnorm(length(t), sd = 0.2)
#' p0 <- c(2.5, 3.2, 0.8, 1.5, 6.8, 1.2)
#'
#' y_obs <- y
#' mask1 <- abs(t - p0[2]) <= 3 * p0[3]
#' mask2 <- abs(t - p0[5]) <= 3 * p0[6]
#' mask <- mask1 | mask2
#' y[!mask] <- 0
#'
#' idx1 <- which(mask1); idx2 <- which(mask2)
#' i <- c(idx1, idx1, idx1, idx2, idx2, idx2)
#' j <- c(rep(1, length(idx1)), rep(2, length(idx1)), rep(3, length(idx1)),
#'        rep(4, length(idx2)), rep(5, length(idx2)), rep(6, length(idx2)))
#' Jpat <- sparseMatrix(i = i, j = j, x = 1, dims = c(length(t), 6))
#'
#' func <- function(p, t, ...) {
#'   f <- g(p[1], p[2], p[3], t) + g(p[4], p[5], p[6], t)
#'   f[!mask] <- 0
#'   f
#' }
#'
#' fjac <- function(p, t, ...) {
#'   A1 <- p[1]; mu1 <- p[2]; s1 <- p[3]
#'   A2 <- p[4]; mu2 <- p[5]; s2 <- p[6]
#'   e1 <- exp(-0.5 * ((t[idx1] - mu1) / s1)^2)
#'   e2 <- exp(-0.5 * ((t[idx2] - mu2) / s2)^2)
#'   J <- Jpat
#'   J@x <- c(e1, A1 * e1 * ((t[idx1] - mu1) / s1^2),
#'            A1 * e1 * ((t[idx1] - mu1)^2 / s1^3),
#'            e2, A2 * e2 * ((t[idx2] - mu2) / s2^2),
#'            A2 * e2 * ((t[idx2] - mu2)^2 / s2^3))
#'   J
#' }
#'
#' fit <- sparselm(p0, y, func, fjac,
#'                 Jnnz = length(i), JtJnnz = -1,
#'                 nconvars = 0, t = t)
#' Jpat
#' fit$par
#'
#' f0 <- g(p0[1], p0[2], p0[3], t) + g(p0[4], p0[5], p0[6], t)
#' f1 <- g(fit$par[1], fit$par[2], fit$par[3], t) + g(fit$par[4], fit$par[5], fit$par[6], t)
#' plot(t, y_obs, pch = 16, cex = 0.6, col = "black", xlab = "t", ylab = "y")
#' lines(t, f0, col = "blue")
#' lines(t, f1, col = "red")
#'
#' @export
sparselm <- function(p, x, func, fjac, Jnnz, JtJnnz=-1, nconvars=0, itmax=100, opts=sparselm.opts(), dif=FALSE, ...) {

	p <- as.numeric(p)
	x <- as.numeric(x)
	nvars <- length(p)
	nobs <- length(x)
	opts_num <- as.numeric(opts)
	eps3 <- if (length(opts_num) >= 4L) abs(opts_num[4L]) else NA_real_
	ctx <- new.env(parent = parent.frame())
	ctx$callback_error <- NULL
	ctx$last_jac <- NULL
	ctx$use_prefetched_jac <- FALSE
	ctx$last_func <- NULL
	ctx$use_prefetched_func <- FALSE
	validate_func_result <- function(value) {
		value <- as.numeric(value)
		if (length(value) != nobs) {
			stop("evaluation of func function not expected length!", call. = FALSE)
		}
		value
	}
	normalize_jac <- function(jac) {
		if (inherits(jac, "dgCMatrix")) {
			jac
		} else {
			methods::as(methods::as(jac, "generalMatrix"), "CsparseMatrix")
		}
	}
	validate_jac <- function(jac) {
		jac <- normalize_jac(jac)
		if (!identical(dim(jac), c(nobs, nvars))) {
			stop("fjac did not return expected dgCMatrix object", call. = FALSE)
		}
		if (length(jac@x) != Jnnz || length(jac@i) != Jnnz || length(jac@p) != nvars + 1L) {
			stop("fjac did not return expected dgCMatrix object", call. = FALSE)
		}
		if (jac@p[1] != 0L || jac@p[nvars + 1L] != Jnnz || any(diff(jac@p) < 0L)) {
			stop("fjac did not return expected dgCMatrix object", call. = FALSE)
		}
		for (col in seq_len(nvars)) {
			start <- jac@p[col] + 1L
			end <- jac@p[col + 1L]
			if (start <= end) {
				rows <- jac@i[start:end]
				if (any(rows < 0L) || any(rows >= nobs) || any(diff(rows) <= 0L)) {
					stop("fjac did not return expected dgCMatrix object", call. = FALSE)
				}
			}
		}
		jac
	}
	prefetched_func <- tryCatch(
		validate_func_result(func(p, ...)),
		error = function(e) {
			ctx$callback_error <- conditionMessage(e)
			NULL
		}
	)
	if (!is.null(ctx$callback_error)) {
		stop(ctx$callback_error, call. = FALSE)
	}
	ctx$last_func <- prefetched_func
	ctx$use_prefetched_func <- TRUE
	need_initial_jac <- nobs >= nvars &&
		all(is.finite(prefetched_func)) &&
		!is.na(eps3) &&
		sum((x - prefetched_func)^2) > eps3

	if (need_initial_jac) {
		prefetched_jac <- tryCatch(
			validate_jac(fjac(p, ...)),
			error = function(e) {
				ctx$callback_error <- conditionMessage(e)
				NULL
			}
		)
		if (!is.null(ctx$callback_error)) {
			stop(ctx$callback_error, call. = FALSE)
		}
		ctx$last_jac <- prefetched_jac
		ctx$use_prefetched_jac <- TRUE
	}

	func1 <- function(p) {
		if (ctx$use_prefetched_func) {
			ctx$use_prefetched_func <- FALSE
			return(ctx$last_func)
		}
		if (!is.null(ctx$callback_error)) {
			return(rep(NaN, length(x)))
		}
		tryCatch(
			{
				value <- validate_func_result(func(p, ...))
				ctx$last_func <- value
				value
			},
			error = function(e) {
				ctx$callback_error <- conditionMessage(e)
				rep(NaN, length(x))
			}
		)
	}
	fjac1 <- function(p) {
		if (ctx$use_prefetched_jac) {
			ctx$use_prefetched_jac <- FALSE
			return(ctx$last_jac)
		}
		if (!is.null(ctx$callback_error)) {
			return(ctx$last_jac)
		}
		tryCatch(
			{
				jac <- validate_jac(fjac(p, ...))
				ctx$last_jac <- jac
				jac
			},
			error = function(e) {
				ctx$callback_error <- conditionMessage(e)
				ctx$last_jac
			}
		)
	}
	
	out <- .Call("C_sparselm", func1, fjac1, p, x, as.integer(nconvars), as.integer(Jnnz), as.integer(JtJnnz), as.integer(itmax), as.numeric(opts), as.logical(dif), ctx, PACKAGE = "sparseLM")
	if (!is.null(ctx$callback_error)) {
		stop(ctx$callback_error, call. = FALSE)
	}

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
