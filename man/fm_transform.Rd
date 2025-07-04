% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crs.R
\name{fm_transform}
\alias{fm_transform}
\alias{fm_transform.default}
\alias{fm_transform.NULL}
\alias{fm_transform.matrix}
\alias{fm_transform.sf}
\alias{fm_transform.sfc}
\alias{fm_transform.sfg}
\alias{fm_transform.Spatial}
\alias{fm_transform.fm_mesh_2d}
\alias{fm_transform.fm_collect}
\alias{fm_transform.fm_lattice_2d}
\alias{fm_transform.fm_segm}
\alias{fm_transform.fm_list}
\title{Object coordinate transformation}
\usage{
fm_transform(x, crs, ...)

\method{fm_transform}{default}(x, crs, ..., crs0 = NULL)

## S3 method for class 'NULL'
fm_transform(x, crs, ...)

\method{fm_transform}{matrix}(x, crs, ..., passthrough = FALSE, crs0 = NULL)

\method{fm_transform}{sf}(x, crs, ..., passthrough = FALSE)

\method{fm_transform}{sfc}(x, crs, ..., passthrough = FALSE)

\method{fm_transform}{sfg}(x, crs, ..., passthrough = FALSE)

\method{fm_transform}{Spatial}(x, crs, ..., passthrough = FALSE)

\method{fm_transform}{fm_mesh_2d}(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = fm_crs(x))

\method{fm_transform}{fm_collect}(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = NULL)

\method{fm_transform}{fm_lattice_2d}(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = fm_crs(x))

\method{fm_transform}{fm_segm}(x, crs = fm_crs(x), ..., passthrough = FALSE, crs0 = fm_crs(x))

\method{fm_transform}{fm_list}(x, crs, ...)
}
\arguments{
\item{x}{The object that should be transformed from it's current CRS to a new CRS}

\item{crs}{The target crs object}

\item{\dots}{Potential additional arguments}

\item{crs0}{The source crs object for spatial classes without crs information}

\item{passthrough}{Default is FALSE.
Setting to TRUE allows objects with no CRS information to be passed
through without transformation. Use with care!}
}
\value{
A transformed object, normally of the same class as
the input object.
}
\description{
Handle transformation of various inla objects according to coordinate
reference systems of \code{crs} (from \code{sf::st_crs()}), \code{fm_crs}, \code{sp::CRS},
\code{fm_CRS}, or \code{INLA::inla.CRS} class.
}
\examples{
fm_transform(
  rbind(c(0, 0), c(0, 90), c(0, 91)),
  crs = fm_crs("sphere"),
  crs0 = fm_crs("longlat_norm")
)

}
\seealso{
\code{\link[=fm_CRS]{fm_CRS()}}
}
