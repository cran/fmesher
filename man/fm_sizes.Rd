% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fem.R
\name{fm_sizes}
\alias{fm_sizes}
\alias{fm_sizes.fm_mesh_2d}
\alias{fm_sizes.fm_mesh_3d}
\title{fm_sizes}
\usage{
fm_sizes(...)

\method{fm_sizes}{fm_mesh_2d}(mesh, ...)

\method{fm_sizes}{fm_mesh_3d}(mesh, ...)
}
\arguments{
\item{...}{Passed on to submethods}

\item{mesh}{object of a supported mesh class}
}
\value{
A \code{list} with elements \code{face} and \code{vertex} for 2D meshes, or \code{cell}
and \code{vertex} for 3D meshes. The elements are vectors of effective sizes of
the faces/cells and vertices, respectively.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
Compute effective sizes of faces/cells and vertices in a mesh
}
\examples{
str(fm_sizes(fmexample$mesh))

}
