% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonconvex_hull.R
\name{fm_nonconvex_hull_inla}
\alias{fm_nonconvex_hull_inla}
\alias{fm_nonconvex_hull_inla_basic}
\title{Non-convex hull computation}
\usage{
fm_nonconvex_hull_inla(
  x,
  convex = -0.15,
  concave = convex,
  resolution = 40,
  eps = NULL,
  eps_rel = NULL,
  crs = NULL,
  ...
)

fm_nonconvex_hull_inla_basic(
  x,
  convex = -0.15,
  resolution = 40,
  eps = NULL,
  crs = NULL
)
}
\arguments{
\item{x}{A spatial object}

\item{convex}{numeric vector; How much to extend}

\item{concave}{numeric vector; The minimum allowed reentrant curvature. Default equal to \code{convex}}

\item{resolution}{The internal computation resolution.  A warning will be
issued when this needs to be increased for higher accuracy, with the
required resolution stated.}

\item{eps, eps_rel}{The polygonal curve simplification tolerances used for
simplifying the resulting boundary curve.  See \code{\link[=fm_simplify_helper]{fm_simplify_helper()}} for
details.}

\item{crs}{Options crs object for the resulting polygon}

\item{\dots}{Unused.}
}
\value{
\code{fm_nonconvex_hull_inla()} returns an \code{fm_segm}/\code{inla.mesh.segment}
object, for compatibility with \code{inla.nonconvex.hull()}.
}
\description{
Legacy method for \code{INLA::inla.nonconvex.hull()}
}
\details{
Requires \code{splancs::nndistF()}
}
\section{Functions}{
\itemize{
\item \code{fm_nonconvex_hull_inla_basic()}: Special method for \code{convex = 0}.

}}
\section{INLA compatibility}{

For mesh and curve creation, the \code{\link[=fm_rcdt_2d_inla]{fm_rcdt_2d_inla()}}, \code{\link[=fm_mesh_2d_inla]{fm_mesh_2d_inla()}},
and \code{\link[=fm_nonconvex_hull_inla]{fm_nonconvex_hull_inla()}} methods will keep the interface syntax used by
\code{INLA::inla.mesh.create()}, \code{INLA::inla.mesh.2d()}, and
\code{INLA::inla.nonconvex.hull()} functions, respectively, whereas the
\code{\link[=fm_rcdt_2d]{fm_rcdt_2d()}}, \code{\link[=fm_mesh_2d]{fm_mesh_2d()}}, and \code{\link[=fm_nonconvex_hull]{fm_nonconvex_hull()}} interfaces may be
different, and potentially change in the future.
}

\examples{
\dontshow{if (require("splancs")) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
fm_nonconvex_hull_inla(cbind(0, 0), convex = 1)
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=fm_nonconvex_hull]{fm_nonconvex_hull()}}

Other nonconvex inla legacy support: 
\code{\link{fm_segm_contour_helper}()},
\code{\link{fm_simplify_helper}()}
}
\concept{nonconvex inla legacy support}
