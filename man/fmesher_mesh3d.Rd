% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fmesher_mesh3d}
\alias{fmesher_mesh3d}
\title{3D tetrahedralisation storage}
\usage{
fmesher_mesh3d(options, loc, tv)
}
\arguments{
\item{options}{list of triangulation options}

\item{loc}{numeric matrix; initial points to include}

\item{tv}{4-column integer matrix with 0-based vertex indices for each triangle}
}
\value{
A list of information objects for a generated tetrahedralisation
}
\description{
(...)
}
\examples{
m <- fmesher_mesh3d(list(),
                    matrix(c(1,0,0,0,1,0,0,0,1,0,0,0), 4, 3, byrow=TRUE),
                    matrix(c(0,1,2,3), 1, 4, byrow=TRUE))
}
