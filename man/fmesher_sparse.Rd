% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fm_as_dgCMatrix}
\alias{fm_as_dgCMatrix}
\alias{fm_as_dgTMatrix}
\alias{fm_as_unpackedMatrix}
\alias{fm_as_fmesher_sparse}
\alias{fm_as_dgCMatrix.default}
\alias{fm_as_dgCMatrix.fmesher_sparse}
\alias{fm_as_dgTMatrix.default}
\alias{fm_as_unpackedMatrix.default}
\alias{fm_as_unpackedMatrix.fmesher_sparse}
\alias{fm_as_dgTMatrix.fmesher_sparse}
\title{Conversion between sparse matrix types}
\usage{
fm_as_dgCMatrix(x)

fm_as_dgTMatrix(x, unique = TRUE, ...)

fm_as_unpackedMatrix(x)

fm_as_fmesher_sparse(x)

\method{fm_as_dgCMatrix}{default}(x)

\method{fm_as_dgCMatrix}{fmesher_sparse}(x)

\method{fm_as_dgTMatrix}{default}(x, unique = TRUE, ...)

\method{fm_as_unpackedMatrix}{default}(x)

\method{fm_as_unpackedMatrix}{fmesher_sparse}(x)

\method{fm_as_dgTMatrix}{fmesher_sparse}(x, unique = TRUE, ...)
}
\arguments{
\item{x}{Object to be converted}

\item{unique}{logical; if \code{TRUE}, ensures that the sparse triplet
representation has a single entry for each non-zero matrix element.}
}
\value{
\code{fm_as_dgCMatrix} returns a \link[Matrix:dgCMatrix-class]{Matrix::dgCMatrix} object.

\code{fm_as_dgTMatrix} returns a \link[Matrix:dgTMatrix-class]{Matrix::dgTMatrix} object.

\code{fm_as_unpackedMatrix} returns an object of virtual class
\link[Matrix:unpackedMatrix-class]{Matrix::unpackedMatrix}.

\code{fm_as_fmesher_sparse} returns an \code{fmesher_sparse} object.
}
\description{
Conversion between sparse matrix types
}
\examples{
library(Matrix)
str(A <- fm_as_dgCMatrix(matrix(c(1, 2, 0, 0, 0, 3, 4, 0, 5), 3, 3)))
str(fm_as_dgTMatrix(A))
str(fm_as_unpackedMatrix(A))
str(fm_as_fmesher_sparse(A))
}
\keyword{internal}
