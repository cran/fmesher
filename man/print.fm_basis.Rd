% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.fm_basis}
\alias{print.fm_basis}
\title{Print method for \code{fm_basis}}
\usage{
\method{print}{fm_basis}(x, ..., prefix = "")
}
\arguments{
\item{x}{\code{\link[=fm_basis]{fm_basis()}} object}

\item{\dots}{Unused}

\item{prefix}{a prefix to be used for each line. Default is an empty string.}
}
\value{
\code{invisible(x)}
}
\description{
Prints information for an \link{fm_basis} object.
}
\examples{
print(fm_basis(fmexample$mesh, fmexample$loc, full = TRUE))
}
\seealso{
\code{\link[=fm_basis]{fm_basis()}}
}
