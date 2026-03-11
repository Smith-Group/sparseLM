# sparseLM R Package

`sparseLM` provides an R interface to the `sparseLM` library for solving
nonlinear least squares problems with sparse Jacobians. The package wraps the
underlying C solver and uses CHOLMOD through the `Matrix` package interface,
exposing the results through an R-friendly fitting interface.

The underlying `sparseLM` algorithms and C library were created by Manolis
Lourakis. The R package interface was created by Colin Smith.

## Documentation

Package documentation is available at:

<https://smith-group.github.io/sparseLM/>

## Upstream project

The original `sparseLM` project page is:

<https://users.ics.forth.gr/~lourakis/sparseLM/>

## Package scope

This package is focused on bringing `sparseLM` into R with:

- a wrapper for sparse Levenberg-Marquardt fitting
- support for user-supplied sparse Jacobians via `Matrix` classes

## Citation

If you use `sparseLM`, please cite the original sparseLM publication:

Lourakis, M. I. A. (2010). Sparse Non-linear Least Squares Optimization for
Geometric Vision. In *Computer Vision - ECCV 2010*, Lecture Notes in Computer
Science 6312, 43-56. Springer.
<https://doi.org/10.1007/978-3-642-15552-9_4>
