## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(12345L)

## ----setup--------------------------------------------------------------------
library(fmesher)

## -----------------------------------------------------------------------------
# Custom class for harmonic functions up to order `n`
create_custom <- function(n) {
  stopifnot(n >= 0)
  structure(
    list(n = n),
    class = "custom"
  )
}
fm_dof.custom <- function(x) {
  # Return the number of degrees of freedom
  1L + 2L * x[["n"]]
}
fm_basis.custom <- function(x, loc, ..., full = FALSE) {
  # Return the evaluated basis functions
  A <- Matrix::Matrix(0.0, NROW(loc), fm_dof(x))
  ok <- !is.na(loc)
  A[ok, 1L] <- 1.0
  for (k in seq_len(x[["n"]])) {
    A[ok, 2 * k] <- cos(2 * pi * k * loc[ok])
    A[ok, 2 * k + 1L] <- sin(2 * pi * k * loc[ok])
  }
  result <- structure(
    list(
      A = A,
      ok = ok, # Required prior to version 0.2.0.9003
      loc = loc
    ),
    class = "fm_basis"
  )

  # Use the fm_basis method to extract the A matrix if full is FALSE:
  fm_basis(result, full = full)
}

## ----eval=FALSE---------------------------------------------------------------
#  # 'matrix' and 'Matrix' methods:
#  fm_basis(
#    A = A,
#    ok = ok, # If missing or NULL, inferred to be all TRUE
#    loc = loc, # Optional additional content
#    full = full
#  )
#  # 'list' method:
#  fm_basis(
#    list(
#      A = A,
#      ok = ok, # If missing or NULL, inferred to be all TRUE
#      loc = loc
#    ),
#    full = full
#  )

## -----------------------------------------------------------------------------
.S3method("fm_dof", "custom", "fm_dof.custom")
.S3method("fm_basis", "custom", "fm_basis.custom")

## -----------------------------------------------------------------------------
#' @rawNamespace S3method(fmesher::fm_dof, custom)
#' @rawNamespace S3method(fmesher::fm_basis, custom)

## ----eval = FALSE-------------------------------------------------------------
#  #' @title Degrees of freedom for custom mesh
#  #' @description the number of degrees of freedom
#  #' # The rest of the documentation goes here
#  #' @exportS3method fmesher::fm_dof
#  fm_dof.custom <- function(x) {
#    1L + 2L * x[["n"]]
#  }

## -----------------------------------------------------------------------------
m <- create_custom(2)

# How many latent variables are needed?
fm_dof(m)

# Evaluate the basis functions at some locations:
fm_basis(m, seq(0, 1, length.out = 6))
fm_basis(m, seq(0, 1, length.out = 6), full = TRUE)

# Check if missing values are handled correctly:
fm_basis(m, c(0.1, NA, 0.2))
fm_basis(m, c(0.1, NA, 0.2), full = TRUE)

## ----echo=FALSE---------------------------------------------------------------
out <- data.frame(
  x = c(
    seq(0, 0.3, by = 0.1),
    seq(0.3, 0.5, by = 0.1),
    seq(0.5, 1, by = 0.1)
  ),
  weight = rep(
    c(0.05, 0.1, 0.05, 0.05, 0.1, 0.05, 0.05, 0.1, 0.05),
    c(1, 2, 1, 1, 1, 1, 1, 4, 1)
  ),
  .block = rep(c(1, 2, 3), c(4, 3, 6))
)
out

## ----eval=TRUE----------------------------------------------------------------
values <- fm_evaluate(
  m,
  field = c(1, 1, 0, 0, 0),
  loc = out[["x"]]
)

# Blockwise aggregation:
fm_block_eval(
  block = out$.block,
  weights = out$weight,
  values = values
)

# Exact integrals:
c(0.3, 0.2, 0.5) +
  c(
    sin(2 * pi * 0.3),
    sin(2 * pi * 0.5) - sin(2 * pi * 0.3),
    sin(2 * pi) - sin(2 * pi * 0.5)
  ) / (2 * pi)

