% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbox.R
\name{fm_bbox}
\alias{fm_bbox}
\alias{fm_bbox.list}
\alias{fm_bbox.NULL}
\alias{fm_bbox.numeric}
\alias{fm_bbox.matrix}
\alias{fm_bbox.Matrix}
\alias{fm_bbox.fm_bbox}
\alias{fm_bbox.fm_mesh_1d}
\alias{fm_bbox.fm_mesh_2d}
\alias{fm_bbox.fm_mesh_3d}
\alias{fm_bbox.fm_segm}
\alias{fm_bbox.fm_lattice_2d}
\alias{fm_bbox.fm_lattice_Nd}
\alias{fm_bbox.fm_tensor}
\alias{fm_bbox.fm_collect}
\alias{fm_bbox.sf}
\alias{fm_bbox.sfg}
\alias{fm_bbox.sfc}
\alias{fm_bbox.bbox}
\alias{fm_as_bbox}
\alias{[.fm_bbox}
\alias{c.fm_bbox}
\alias{fm_as_bbox_list}
\title{Bounding box class}
\usage{
fm_bbox(...)

\method{fm_bbox}{list}(x, ...)

## S3 method for class 'NULL'
fm_bbox(...)

\method{fm_bbox}{numeric}(x, ...)

\method{fm_bbox}{matrix}(x, ...)

\method{fm_bbox}{Matrix}(x, ...)

\method{fm_bbox}{fm_bbox}(x, ...)

\method{fm_bbox}{fm_mesh_1d}(x, ...)

\method{fm_bbox}{fm_mesh_2d}(x, ...)

\method{fm_bbox}{fm_mesh_3d}(x, ...)

\method{fm_bbox}{fm_segm}(x, ...)

\method{fm_bbox}{fm_lattice_2d}(x, ...)

\method{fm_bbox}{fm_lattice_Nd}(x, ...)

\method{fm_bbox}{fm_tensor}(x, ...)

\method{fm_bbox}{fm_collect}(x, ...)

\method{fm_bbox}{sf}(x, ...)

\method{fm_bbox}{sfg}(x, ...)

\method{fm_bbox}{sfc}(x, ...)

\method{fm_bbox}{bbox}(x, ...)

fm_as_bbox(x, ...)

\method{[}{fm_bbox}(x, i)

\method{c}{fm_bbox}(..., .join = FALSE)

fm_as_bbox_list(x, ...)
}
\arguments{
\item{...}{Passed on to sub-methods}

\item{x}{\code{fm_bbox} object from which to extract element(s)}

\item{i}{indices specifying elements to extract}

\item{.join}{logical; if \code{TRUE}, concatenate the bounding boxes into a single
multi-dimensional bounding box. Default is \code{FALSE}.}
}
\value{
For \code{\link[=c.fm_bbox]{c.fm_bbox()}}, a \code{fm_bbox_list} object if \code{.join = FALSE} (the
default) or an \code{fm_bbox} object if \code{.join = TRUE}.
}
\description{
Simple class for handling bounding box information
}
\section{Methods (by class)}{
\itemize{
\item \code{fm_bbox(list)}: Construct a bounding box from
precomputed interval information, stored as a list of 2-vector ranges,
\code{list(xlim, ylim, ...)}.

}}
\section{Methods (by generic)}{
\itemize{
\item \code{[}: Extract sub-list

\item \code{c(fm_bbox)}: The \code{...} arguments should be \code{fm_bbox} objects, or
coercible with \code{fm_as_bbox(list(...))}.

}}
\section{Functions}{
\itemize{
\item \code{fm_as_bbox_list()}: Convert a list to a \code{fm_bbox_list} object, with
each element converted to an \code{fm_bbox} object.

}}
\examples{
fm_bbox(matrix(1:6, 3, 2))
m <- c(A = fm_bbox(cbind(1, 2)), B = fm_bbox(cbind(3, 4)))
str(m)
str(m[2])
m <- fm_as_bbox_list(list(
  A = fm_bbox(cbind(1, 2)),
  B = fm_bbox(cbind(3, 4))
))
str(fm_as_bbox_list(m))
}
