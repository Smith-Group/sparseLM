% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseLM.R
\name{sparselm.opts}
\alias{sparselm.opts}
\title{Sparse Levenberg-Marquardt Algorithm Options}
\usage{
sparselm.opts(
  mu = 0.001,
  epsilon1 = 1e-12,
  epsilon2 = 1e-12,
  epsilon3 = 1e-12,
  delta = 1e-06
)
}
\arguments{
\item{mu}{scale factor for initial mu}

\item{epsilon1}{stopping threshold for ||J^T e||_inf}

\item{epsilon2}{stopping threshold for ||dp||_2}

\item{epsilon3}{stopping threshold for ||e||_2}

\item{delta}{step used in difference approximation to the Jacobian}
}
\value{
numeric vector of length 6 with the above options and spsolver=1 (SuiteSparse CHOLMOD)
}
\description{
Sparse Levenberg-Marquardt Algorithm Options
}
