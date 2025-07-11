% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis.R
\name{fm_block}
\alias{fm_block}
\alias{fm_block_eval}
\alias{fm_block_logsumexp_eval}
\alias{fm_block_weights}
\alias{fm_block_log_weights}
\alias{fm_block_log_shift}
\alias{fm_block_prep}
\title{Blockwise aggregation matrices}
\usage{
fm_block(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL
)

fm_block_eval(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL,
  values = NULL
)

fm_block_logsumexp_eval(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL,
  values = NULL,
  log = TRUE
)

fm_block_weights(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL
)

fm_block_log_weights(
  block = NULL,
  weights = NULL,
  log_weights = NULL,
  rescale = FALSE,
  n_block = NULL
)

fm_block_log_shift(block = NULL, log_weights = NULL, n_block = NULL)

fm_block_prep(
  block = NULL,
  log_weights = NULL,
  weights = NULL,
  n_block = NULL,
  values = NULL,
  n_values = NULL,
  force_log = FALSE
)
}
\arguments{
\item{block}{integer vector; block information. If \code{NULL},
\code{rep(1L, block_len)} is used, where \code{block_len} is determined by
\verb{length(log_weights)))} or \verb{length(weights)))}. A single scalar is also
repeated to a vector of corresponding length to the weights.

Note: from version \verb{0.2.0.9017} to \verb{0.4.0.9005}, 'character'
input was converted to integer with \code{as.integer(factor(block))}. As this
could lead to unintended ordering of the output, this is no longer allowed.}

\item{weights}{Optional weight vector}

\item{log_weights}{Optional \code{log(weights)} vector. Overrides \code{weights} when
non-NULL.}

\item{rescale}{logical; If \code{TRUE}, normalise the weights by \code{sum(weights)}
or \code{sum(exp(log_weights))} within each block.
Default: \code{FALSE}}

\item{n_block}{integer; The number of conceptual blocks. Only needs to be
specified if it's larger than \code{max(block)}, or to keep the output of
consistent size for different inputs.}

\item{values}{Vector to be blockwise aggregated}

\item{log}{If \code{TRUE} (default), return log-sum-exp. If \code{FALSE},
return sum-exp.}

\item{n_values}{When supplied, used instead of \code{length(values)} to determine
the value vector input length.}

\item{force_log}{When \code{FALSE} (default),
passes either \code{weights} and \code{log_weights} on, if provided, with \code{log_weights}
taking precedence. If \code{TRUE}, forces the computation of \code{log_weights},
whether given in the input or not.}
}
\value{
A (sparse) matrix
}
\description{
Creates an aggregation matrix for blockwise aggregation, with optional
weighting.
}
\section{Functions}{
\itemize{
\item \code{fm_block()}: A (sparse) matrix of size \code{n_block} times
\code{length(block)}.

\item \code{fm_block_eval()}: Evaluate aggregation. More efficient alternative to to
\code{as.vector(fm_block(...) \%*\% values)}.

\item \code{fm_block_logsumexp_eval()}: Evaluate log-sum-exp aggregation.
More efficient and numerically stable alternative to to
\code{log(as.vector(fm_block(...) \%*\% exp(values)))}.

\item \code{fm_block_weights()}: Computes (optionally) blockwise renormalised weights

\item \code{fm_block_log_weights()}: Computes (optionally) blockwise renormalised log-weights

\item \code{fm_block_log_shift()}: Computes shifts for stable blocked log-sum-exp.
To compute \eqn{\log(\sum_{i; \textrm{block}_i=k} \exp(v_i) w_i)}{
log(sum_(i;block_i=k) exp(v_i) w_i)
} for
each block \code{k}, first compute combined values and weights, and a shift:

\if{html}{\out{<div class="sourceCode">}}\preformatted{w_values <- values + fm_block_log_weights(block, log_weights = log_weights)
shift <- fm_block_log_shift(block, log_weights = w_values)
}\if{html}{\out{</div>}}

Then aggregate the values within each block:

\if{html}{\out{<div class="sourceCode">}}\preformatted{agg <- aggregate(exp(w_values - shift[block]),
                 by = list(block = block),
                 \\(x) log(sum(x)))
agg$x <- agg$x + shift[agg$block]
}\if{html}{\out{</div>}}

The implementation uses a faster method:

\if{html}{\out{<div class="sourceCode">}}\preformatted{as.vector(
  Matrix::sparseMatrix(
    i = block,
    j = rep(1L, length(block)),
    x = exp(w_values - shift[block]),
    dims = c(n_block, 1))
) + shift
}\if{html}{\out{</div>}}

\item \code{fm_block_prep()}: Helper function for preparing \code{block}, \code{weights}, and
\code{log_weights}, \code{n_block} inputs.

}}
\examples{
block <- rep(1:2, 3:2)
fm_block(block)
fm_block(block, rescale = TRUE)
fm_block(block, log_weights = -2:2, rescale = TRUE)
fm_block_eval(
  block,
  weights = 1:5,
  rescale = TRUE,
  values = 11:15
)
fm_block_logsumexp_eval(
  block,
  weights = 1:5,
  rescale = TRUE,
  values = log(11:15),
  log = FALSE
)
}
