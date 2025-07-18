% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{fm_int}
\alias{fm_int}
\alias{fm_int.list}
\alias{fm_int.numeric}
\alias{fm_int.character}
\alias{fm_int.factor}
\alias{fm_int.SpatRaster}
\alias{fm_int.fm_lattice_2d}
\alias{fm_int.fm_mesh_1d}
\alias{fm_int.fm_mesh_2d}
\title{Multi-domain integration}
\usage{
fm_int(domain, samplers = NULL, ...)

\method{fm_int}{list}(domain, samplers = NULL, ...)

\method{fm_int}{numeric}(domain, samplers = NULL, name = "x", ...)

\method{fm_int}{character}(domain, samplers = NULL, name = "x", ...)

\method{fm_int}{factor}(domain, samplers = NULL, name = "x", ...)

\method{fm_int}{SpatRaster}(domain, samplers = NULL, name = "x", ...)

\method{fm_int}{fm_lattice_2d}(domain, samplers = NULL, name = "x", ...)

\method{fm_int}{fm_mesh_1d}(
  domain,
  samplers = NULL,
  name = "x",
  int.args = NULL,
  format = NULL,
  ...
)

\method{fm_int}{fm_mesh_2d}(
  domain,
  samplers = NULL,
  name = NULL,
  int.args = NULL,
  format = NULL,
  ...
)
}
\arguments{
\item{domain}{Functional space specification; single domain or a named list
of domains}

\item{samplers}{For single domain \code{fm_int} methods, an object specifying one
or more subsets of the domain, and optional weighting in a \code{weight}
variable. For \code{fm_int.list}, a list of sampling definitions, where data
frame elements may contain information for multiple domains, in which case
each row represent a separate tensor product integration subspace.}

\item{\dots}{Additional arguments passed on to other methods}

\item{name}{For single-domain methods, the variable name to use for the
integration points. Default 'x'}

\item{int.args}{List of arguments passed to line and integration methods.
\itemize{
\item \code{method}: "stable" (to aggregate integration weights onto mesh nodes)
or "direct" (to construct a within triangle/segment integration scheme
without aggregating onto mesh nodes)
\item \code{nsub1}, \code{nsub2}: integers controlling the number of internal integration
points before aggregation. Points per triangle: \code{(nsub2+1)^2}.
Points per knot segment: \code{nsub1}
}}

\item{format}{character; determines the output format, as either "sf"
(default for \code{fm_mesh_2d} when the sampler is \code{NULL}),
"numeric" (default for \code{fm_mesh_1d}), "bary", or "sp".
When \code{NULL}, determined by the domain and sampler types.}
}
\value{
A \code{tibble}, \code{sf}, or \code{SpatialPointsDataFrame} of 1D
and 2D integration points, including a \code{weight} column, a\code{.block} column,
and a matrix column \code{.block_origin}.
The \code{.block} column is used to identify the integration
blocks defined by the samplers. The \code{.block_origin} collects the original
subdomain block information for tensor product blocks.
}
\description{
Construct integration points on tensor product spaces
}
\section{Methods (by class)}{
\itemize{
\item \code{fm_int(list)}: Multi-domain integration

\item \code{fm_int(numeric)}: Discrete double or integer space integration

\item \code{fm_int(character)}: Discrete character space integration

\item \code{fm_int(factor)}: Discrete factor space integration

\item \code{fm_int(SpatRaster)}: \code{SpatRaster} integration. Not yet implemented.

\item \code{fm_int(fm_lattice_2d)}: \code{fm_lattice_2d} integration. Not yet implemented.

\item \code{fm_int(fm_mesh_1d)}: \code{fm_mesh_1d} integration. Supported samplers:
\itemize{
\item \code{NULL} for integration over the entire domain;
\item A length 2 vector defining an interval;
\item A 2-column matrix with a single interval in each row;
\item A tibble with a named column containing a matrix, and optionally a
\code{weight} column.
}

\item \code{fm_int(fm_mesh_2d)}: \code{fm_mesh_2d} integration. Any sampler class with an
associated \code{\link[=fm_int_mesh_2d]{fm_int_mesh_2d()}} method is supported.

}}
\examples{
# Integration on the interval (2, 3.5) with Simpson's rule
ips <- fm_int(fm_mesh_1d(0:4), samplers = cbind(2, 3.5))
plot(ips$x, ips$weight)

# Create integration points for the two intervals [0,3] and [5,10]
ips <- fm_int(
  fm_mesh_1d(0:10),
  rbind(c(0, 3), c(5, 10))
)
plot(ips$x, ips$weight)

# Convert a 1D mesh into integration points
mesh <- fm_mesh_1d(seq(0, 10, by = 1))
ips <- fm_int(mesh, name = "time")
plot(ips$time, ips$weight)

if (require("ggplot2", quietly = TRUE)) {
  #' Integrate on a 2D mesh with polygon boundary subset
  ips <- fm_int(fmexample$mesh, fmexample$boundary_sf[[1]])
  ggplot() +
    geom_sf(data = fm_as_sfc(fmexample$mesh, multi = TRUE), alpha = 0.5) +
    geom_sf(data = fmexample$boundary_sf[[1]], fill = "red", alpha = 0.5) +
    geom_sf(data = ips, aes(size = weight)) +
    scale_size_area()
}

# Individual sampling points:
(ips <- fm_int(0:10, c(0, 3, 5, 6, 10)))
# Sampling blocks:
(ips <- fm_int(0:10, list(c(0, 3), c(5, 6, 10))))

# Continuous integration on intervals
ips <- fm_int(
  fm_mesh_1d(0:10, boundary = "cyclic"),
  rbind(c(0, 3), c(5, 10))
)
plot(ips$x, ips$weight)

}
