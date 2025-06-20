#' @include deprecated.R

# fm_mesh_1d ####

`match.arg.vector` <- function(arg = NULL,
                               choices,
                               length = NULL) {
  ## Like match.arg, but for a vector of options 'arg'
  if (is.null(length)) {
    length <- ifelse(is.null(arg), 1L, length(arg))
  }
  if (is.null(arg)) {
    arg <- match.arg(arg, choices)
  } else {
    for (k in seq_along(arg)) {
      arg[k] <- match.arg(arg[k], choices)
    }
  }
  if (length(arg) < length) {
    arg <- c(arg, rep(arg, length - length(arg)))
  } else if (length(arg) > length) {
    stop("Option list too long.")
  }
  return(arg)
}

#' @title Make a 1D mesh object
#' @description
#' Create a `fm_mesh_1d` object.
#'
#' @param loc B-spline knot locations.
#' @param interval Interval domain endpoints.
#' @param boundary Boundary condition specification.  Valid conditions are
#' `c('neumann', 'dirichlet', 'free', 'cyclic')`.  Two separate values can
#' be specified, one applied to each endpoint.
#' @param degree The B-spline basis degree.  Supported values are 0, 1, and 2.
#' @param free.clamped If `TRUE`, for `'free'` boundaries, clamp the
#' basis functions to the interval endpoints.
#' @param \dots Additional options, currently unused.
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#' @returns An `fm_mesh_1d` object
#' @export
#' @family object creation and conversion
#' @examples
#' if (require("ggplot2")) {
#'   m1 <- fm_mesh_1d(c(1, 2, 3, 5, 8, 10),
#'     boundary = c("neumann", "free")
#'   )
#'   weights <- c(2, 3, 6, 3, 4, 7)
#'   ggplot() +
#'     geom_fm(data = m1, xlim = c(0.5, 11), weights = weights)
#'
#'   m2 <- fm_mesh_1d(c(1, 2, 3, 5, 8, 10),
#'     boundary = c("neumann", "free"),
#'     degree = 2
#'   )
#'   ggplot() +
#'     geom_fm(data = m2, xlim = c(0.5, 11), weights = weights)
#'
#'   # The knot interpretation is different for degree=2 and degree=1 meshes:
#'   ggplot() +
#'     geom_fm(data = m1, xlim = c(0.5, 11), weights = weights) +
#'     geom_fm(data = m2, xlim = c(0.5, 11), weights = weights)
#'
#'   # The `mid` values are the representative basis function midpoints,
#'   # and can be used to connect degree=2 and degree=1 mesh interpretations:
#'   m1b <- fm_mesh_1d(m2$mid,
#'     boundary = c("neumann", "free"),
#'     degree = 1
#'   )
#'   ggplot() +
#'     geom_fm(data = m2, xlim = c(0.5, 11), weights = weights) +
#'     geom_fm(data = m1b, xlim = c(0.5, 11), weights = weights)
#' }
#'
fm_mesh_1d <- function(loc,
                       interval = range(loc),
                       boundary = NULL,
                       degree = 1,
                       free.clamped = FALSE,
                       ...) {
  ## Note: do not change the order of these options without also
  ## changing 'basis.reduction' below.
  boundary.options <- c("neumann", "dirichlet", "free", "cyclic")

  boundary <- match.arg.vector(boundary, boundary.options, length = 2)
  cyclic <- !is.na(pmatch(boundary[1], "cyclic"))
  if (cyclic && is.na(pmatch(boundary[2], "cyclic"))) {
    stop("Inconsistent boundary specification 'boundary=c(",
      paste(boundary, collapse = ","), ")'.",
      sep = ""
    )
  }

  loc.orig <- loc
  if (cyclic) {
    if (diff(interval) < diff(range(loc))) {
      warning(
        "Given cyclic interval is narrower than the range of knot locations."
      )
    }
    loc_1 <- min(loc)
    if (loc_1 < interval[1]) {
      # Keep the point to the left of the interval, but adjacent
      loc_1 <- (loc_1 - interval[1]) %% diff(interval) - diff(interval) +
        interval[1]
    }
    if (loc_1 > interval[2]) {
      # Move the point into the interval
      loc_1 <- (loc_1 - interval[1]) %% diff(interval) + interval[1]
    }
    loc <- sort(unique((loc - loc_1) %% diff(interval))) + loc_1
  } else {
    if (loc[1] < interval[1]) {
      warning(
        "fm_mesh_1d: All 'loc' should be >= interval[1].",
        " Moving to interval edge."
      )
    }
    if (loc[2] > interval[2]) {
      warning(
        "fm_mesh_1d: All 'loc' should be <= interval[2].",
        " Moving to interval edge."
      )
    }
    if (min(loc) > interval[1]) {
      warning(
        "fm_mesh_1d: 'min(loc)' should be == interval[1].",
        " Adding knot at interval edge."
      )
    }
    if (max(loc) < interval[2]) {
      warning(
        "fm_mesh_1d: 'max(loc)' should be == interval[2].",
        " Adding knot at interval edge."
      )
    }
    loc <-
      sort(unique(c(
        interval,
        pmax(
          interval[1], pmin(interval[2], loc)
        )
      )))
  }

  n <- length(loc)

  if ((degree < 0) || (degree > 3)) {
    stop(paste("'degree' must be 0, 1, 2, 3.  'degree=",
      degree,
      "' is not supported.",
      sep = ""
    ))
  }

  if (length(free.clamped) == 1L) {
    free.clamped <- rep(free.clamped, 2)
  }


  ## Number of basis functions
  stopifnot(degree >= 0)
  stopifnot(degree <= 3)
  m_free <- n + c(0, 0, 1, 2)
  # How many more basis functions at each end, compared with the number of
  # knots; for cyclic, is only applied once:
  m_adjust <- list(
    c(0, -1, 0, 0), ## neu, dir, free, cyclic
    c(0, -1, 0, 0), ## neu, dir, free, cyclic
    c(-1, -1, 0, -1), ## neu, dir, free, cyclic
    c(-1, -1, 0, -2) ## neu, dir, free, cyclic
  )
  if (cyclic) {
    m <- m_free[degree + 1] + m_adjust[[degree + 1]][4]
  } else {
    i1 <- pmatch(boundary[1], boundary.options)
    i2 <- pmatch(boundary[2], boundary.options)
    m <- m_free[degree + 1] +
      m_adjust[[degree + 1]][i1] +
      m_adjust[[degree + 1]][i2]
  }
  if (m < 1L) {
    stop("Degree ", degree,
      " meshes must have at least ", 1L,
      " basis functions, not 'm=", m, "'.",
      sep = ""
    )
  }

  ## Compute representative basis midpoints.
  if ((degree == 0) || (degree == 1)) {
    mid <- loc
    if (boundary[1] == "dirichlet") {
      mid <- mid[-1]
    }
    if (boundary[2] == "dirichlet") {
      mid <- mid[-length(mid)]
    }
  } else if (degree == 3) {
    mid <- loc
    if (boundary[1] == "dirichlet") {
      mid <- mid[-1]
    }
    if (boundary[2] == "dirichlet") {
      mid <- mid[-length(mid)]
    }
  } else { ## degree==2
    if (cyclic) {
      mid <- (loc + c(loc[-1], interval[2])) / 2
    } else {
      mid <- (loc[-n] + loc[-1]) / 2
      mid <-
        switch(boundary[1],
          neumann = mid,
          dirichlet = mid,
          free = if (free.clamped[1]) {
            c(loc[1], mid)
          } else {
            c(loc[1] - diff(loc[1:2]) / 2, mid)
          }
        )
      mid <-
        switch(boundary[2],
          neumann = mid,
          dirichlet = mid,
          free = if (free.clamped[2]) {
            c(mid, loc[n])
          } else {
            c(mid, loc[n] + diff(loc[(n - 1):n]) / 2)
          }
        )
    }
  }

  mesh <-
    structure(
      list(
        n = n,
        m = m,
        loc = loc,
        mid = mid,
        interval = interval,
        boundary = boundary,
        cyclic = cyclic,
        manifold = ifelse(cyclic, "S1", "R1"),
        degree = degree,
        free.clamped = free.clamped,
        idx = list(loc = NULL)
      ),
      class = c("fm_mesh_1d", "inla.mesh.1d")
    )

  if (degree != 2) {
    mesh$idx$loc <-
      fm_bary(mesh, loc.orig, method = "nearest")$index
  } else {
    if (length(mid) >= 2) {
      mesh$idx$loc <-
        fm_bary(fm_mesh_1d(mid, degree = 0),
          loc.orig,
          method = "nearest"
        )$index
    } else {
      mesh$idx$loc <- rep(1, length(loc.orig))
    }
  }

  return(mesh)
}



#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_mesh_1d Convert an object to `fm_mesh_1d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @returns An `fm_mesh_1d` or `fm_mesh_1d_list` object
#' @export
#' @family object creation and conversion
#' @export
#' @examples
#' fm_as_mesh_1d_list(list(fm_mesh_1d(1:4)))
fm_as_mesh_1d <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("fm_as_mesh_1d")
}
#' @describeIn fm_as_mesh_1d Convert each element of a list
#' @export
fm_as_mesh_1d_list <- function(x, ...) {
  fm_as_list(x, ..., .class_stub = "mesh_1d")
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
fm_as_mesh_1d.fm_mesh_1d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_1d", setdiff(class(x), "fm_mesh_1d"))
  x
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
#' @method fm_as_mesh_1d inla.mesh.1d
fm_as_mesh_1d.inla.mesh.1d <- function(x, ...) {
  class(x) <- c("fm_mesh_1d", class(x))
  x
}
