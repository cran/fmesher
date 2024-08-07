% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluator.R
\name{fm_basis}
\alias{fm_basis}
\alias{fm_basis.default}
\alias{fm_basis.fm_mesh_1d}
\alias{fm_basis.fm_mesh_2d}
\alias{fm_basis.inla.mesh.1d}
\alias{fm_basis.inla.mesh}
\alias{fm_basis.fm_evaluator}
\alias{fm_basis.fm_tensor}
\title{Compute mapping matrix between mesh function space and points}
\usage{
fm_basis(x, ...)

\method{fm_basis}{default}(x, loc, ...)

\method{fm_basis}{fm_mesh_1d}(x, loc, weights = NULL, derivatives = NULL, ...)

\method{fm_basis}{fm_mesh_2d}(x, loc, weights = NULL, derivatives = NULL, ...)

\method{fm_basis}{inla.mesh.1d}(x, loc, ...)

\method{fm_basis}{inla.mesh}(x, loc, ...)

\method{fm_basis}{fm_evaluator}(x, ...)

\method{fm_basis}{fm_tensor}(x, loc, weights = NULL, ...)
}
\arguments{
\item{x}{An object supported by the \code{\link[=fm_evaluator]{fm_evaluator()}} class}

\item{\dots}{Currently unused}

\item{loc}{A set of points of a class supported by \code{fm_evaluator(x, loc = loc)}}

\item{weights}{Optional weight vector to apply (from the left, one
weight for each row of the basis matrix)}

\item{derivatives}{If non-NULL and logical, return a list, optionally
including derivative matrices.}
}
\value{
A \code{sparseMatrix}

For \code{fm_mesh_1d}, a matrix, or if \code{derivatives} is \code{TRUE},
a list with elements
\item{A }{The projection matrix, \verb{u(loc_i)=sum_j A_ij w_i}}
\item{d1A, d2A }{Derivative weight matrices,
\verb{du/dx(loc_i)=sum_j dx_ij w_i}, etc.}

For \code{fm_mesh_2d}, a matrix, or if \code{derivatives} is \code{TRUE},
a list with elements
\item{A }{The projection matrix, \verb{u(loc_i)=sum_j A_ij w_i}}
\item{dx, dy, dz }{Derivative weight matrices, \verb{du/dx(loc_i)=sum_j dx_ij w_i}, etc.}

For \code{fm_tensor}, a matrix
}
\description{
Computes the basis mapping matrix between a function space on a mesh, and locations.
}
\examples{
# Compute basis mapping matrix
str(fm_basis(fmexample$mesh, fmexample$loc))
}
\seealso{
\code{\link[=fm_raw_basis]{fm_raw_basis()}}
}
