#' Print objects
#'
#' Print objects
#'
#' @param x an object used to select a method.
#' @param \dots further arguments passed to or from other methods.
#' @param verbose logical
#' @param digits a positive integer indicating how many significant digits are
#'   to be used for numeric and complex x. The default, NULL, uses
#'   `getOption("digits")`.
#' @param newline logical; if `TRUE` (default), end the printing with `\n`
#'
#' @returns The input object `x`
#' @name fmesher-print
#' @examples
#' fm_bbox(matrix(1:6, 3, 2))
#' print(fm_bbox(matrix(1:6, 3, 2)), verbose = FALSE)
#'
#' print(fmexample$mesh)
#' print(fmexample$boundary_fm)
#'
#' print(fm_mesh_1d(c(1, 2, 3, 5, 7), degree = 2))
NULL

#' @export
#' @rdname fmesher-print
print.fm_segm <- function(x,
                          ...,
                          digits = NULL,
                          verbose = TRUE,
                          newline = TRUE) {
  my.segm <- function(x) {
    if (is.null(x)) {
      return(list(n = 0, grps = NULL))
    }
    n <- max(0, nrow(x$idx))
    if (max(0, length(unique(x$grp))) > 0) {
      grps <- unique(x$grp)
    } else {
      grps <- NULL
    }
    list(n = n, grps = grps)
  }

  ret <- my.segm(x)

  if (verbose) {
    cat("fm_segm object:\n", sep = "")
  }
  if (verbose && !fm_crs_is_null(fm_crs(x))) {
    crs <- fm_wkt(x$crs)
    ret <- c(ret, list(crs = as.character(fm_wkt(crs))))
    ret <- c(ret, list(crs_proj4 = as.character(fm_proj4string(crs))))

    cat("  CRS:\n    LegacyPROJ4:\t", ret$crs_proj4, "\n", sep = "")
    if (verbose) {
      if (is.na(ret$crs)) {
        cat("    WKT:\tNA\n", sep = "")
      } else {
        cat("    WKT:\n", ret$crs, "\n", sep = "")
      }
    } else {
      cat("    WKT: (only shown with verbose = TRUE)", "\n", sep = "")
    }
  }
  if (verbose) {
    cat("  ", sep = "")
  }
  if (ret$n > 0) {
    extra <- ""
    if (sum(x$is.bnd) > 0) {
      cat(sum(x$is.bnd),
        " boundary edge",
        if (sum(x$is.bnd) != 1) {
          "s"
        },
        sep = ""
      )
      extra <- ", "
    }
    if (sum(!x$is.bnd) > 0) {
      cat(extra, sum(!x$is.bnd),
        " interior edge",
        if (sum(!x$is.bnd) != 1) {
          "s"
        },
        sep = ""
      )
    }
  } else {
    cat("0 edges", sep = "")
  }
  if (!is.null(ret$grps)) {
    n <- length(ret$grps)
    cat(" (", n, " group", if (n != 1) "s" else "", sep = "")
    if (n <= 10) {
      cat(":", ret$grps, sep = " ")
    } else {
      cat(":", ret$grps[1:10], "...", sep = " ")
    }
    cat(")")
  }
  if (verbose) {
    cat("\n  ", sep = "")
    cat("Bounding box = ", sep = "")
    print(
      fm_bbox(x),
      verbose = FALSE,
      newline = FALSE,
      digits = digits
    )
  }
  if (newline) {
    cat("\n", sep = "")
  }
  return(invisible(x))
}


#' @export
#' @rdname fmesher-print
print.fm_segm_list <- function(x,
                               ...,
                               digits = NULL,
                               verbose = FALSE,
                               newline = TRUE) {
  if (verbose) {
    cat("list of ", length(x), " fm_segm objects:\n", sep = "")
    if (!is.null(names(x)[k])) {
      cat(names(x)[k], ": ", sep = "")
    }
    lapply(x, function(xx) {
      print(
        xx,
        digits = digits,
        verbose = TRUE,
        newline = TRUE
      )
    })
  } else {
    for (k in seq_along(x)) {
      if (!is.null(names(x)[k])) {
        cat(names(x)[k], ": ", sep = "")
      }
      print(
        x[[k]],
        digits = digits,
        verbose = FALSE,
        newline = newline
      )
      if (!newline && (k < length(x))) {
        cat(", ", sep = "")
      }
    }
  }
  return(invisible(x))
}


#' @export
#' @rdname fmesher-print
print.fm_list <- function(x,
                          ...,
                          digits = NULL,
                          verbose = FALSE,
                          newline = TRUE) {
  if (verbose) {
    cat("list of ", length(x), " fmesher objects:\n", sep = "")
    lapply(x, function(xx) {
      print(
        xx,
        digits = digits,
        verbose = TRUE,
        newline = TRUE
      )
    })
  } else {
    for (k in seq_along(x)) {
      print(
        x[[k]],
        digits = digits,
        verbose = FALSE,
        newline = newline
      )
      if (!newline && (k < length(x))) {
        cat(", ", sep = "")
      }
    }
  }
  return(invisible(x))
}


#' @param verbose logical
#' @param digits a positive integer indicating how many significant digits are
#'   to be used for numeric and complex x. The default, NULL, uses
#'   `getOption("digits")`.
#'
#' @export
#' @rdname fmesher-print
print.fm_mesh_2d <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  if (verbose) {
    ret <- c(ret, list(time = x$meta$time, is.refined = x$meta$is.refined))
  }
  ret <-
    c(
      ret,
      list(
        manifold = x$manifold,
        nV = nrow(x$loc),
        nT = nrow(x$graph$tv)
      )
    )
  crs <- fm_wkt(x$crs)
  ret <- c(ret, list(crs = as.character(fm_wkt(crs))))
  ret <- c(ret, list(crs_proj4 = as.character(fm_proj4string(crs))))

  if (!is.null(x$segm)) {
    ret$segm <- fm_as_segm_list(list(
      "Boundary" = x$segm$bnd,
      "Interior" = x$segm$int
    ))
  } else {
    ret$segm <- fm_as_segm_list(list(
      Boundary = fm_segm(is.bnd = TRUE),
      Interior = fm_segm(is.bnd = FALSE)
    ))
  }

  my.print.proc_time <- function(x, ...) {
    if (is.null(x)) {
      y <- matrix(0, 1, 5)
    } else if (!is.matrix(x)) {
      y <- matrix(x, 1, 5)
    } else {
      y <- x
    }
    for (k in seq_len(nrow(y))) {
      if (!is.na(y[k, 4L])) {
        y[k, 1L] <- y[k, 1L] + y[k, 4L]
      }
      if (!is.na(y[k, 5L])) {
        y[k, 2L] <- y[k, 2L] + y[k, 5L]
      }
    }
    y <- y[, 1L:3L, drop = FALSE]
    colnames(y) <- c(
      gettext("user"), gettext("system"),
      gettext("elapsed")
    )
    print(y, ..., digits = digits)
    invisible(x)
  }

  cat("fm_mesh_2d object:\n", sep = "")
  if (!fm_crs_is_null(x$crs)) {
    cat("  CRS:\n    LegacyPROJ4:\t", ret$crs_proj4, "\n", sep = "")
    if (ret$verbose) {
      if (is.na(ret$crs)) {
        cat("    WKT:\tNA\n", sep = "")
      } else {
        cat("    WKT:\n", ret$crs, "\n", sep = "")
      }
    } else {
      cat("    WKT: (only shown with verbose = TRUE)", "\n", sep = "")
    }
  }
  if (ret$verbose) {
    cat("  Timings:\n")
    my.print.proc_time(ret$time)
  }
  cat("  Manifold:\t", ret$manifold, "\n", sep = "")
  if (ret$verbose) {
    cat("  Refined:\t", ret$is.refined, "\n", sep = "")
  }
  nV <- ret$nV
  nE <- as.integer(sum(x$graph$vv) / 2L)
  nF <- ret$nT
  cat("  V / E / T:\t", as.character(ret$nV), " / ", sep = "")
  cat(as.character(nE), " / ", sep = "")
  cat(as.character(ret$nT), "\n", sep = "")
  cat("  Euler char.:\t", as.character(nV - nE + nF), "\n", sep = "")

  cat("  Constraints:\t")
  print(ret$segm, newline = FALSE, digits = digits, verbose = FALSE)
  cat("\n  ", sep = "")
  print(fm_bbox(x), digits = digits)
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")
  invisible(x)
}



#' @export
#' @rdname fmesher-print
print.fm_mesh_3d <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  if (verbose) {
    ret <- c(ret, list())
  }
  ret <-
    c(
      ret,
      list(
        manifold = fm_manifold(x),
        nV = nrow(x$loc),
        nT = nrow(x$graph$tv)
      )
    )

  cat("fm_mesh_3d object:\n", sep = "")
  cat("  Manifold:\t", ret$manifold, "\n", sep = "")
  nV <- ret$nV
  nE <- as.integer(sum(x$graph$vv) / 2L)
  nF <- nrow(x$graph$tv) * 4L - sum(!is.na(x$graph$tt)) / 2L
  nC <- ret$nT
  cat(
    "  V / E / T / Tet:\t",
    paste0(c(nV, nE, nF, nC), collapse = " / ", sep = ""),
    "\n",
    sep = ""
  )
  cat("  Euler char.:\t", as.character(nV - nE + nF - nC), sep = "")
  cat("\n  ", sep = "")
  print(fm_bbox(x), digits = digits)
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")
  invisible(x)
}



#' @param verbose logical
#'
#' @export
#' @rdname fmesher-print
print.fm_mesh_1d <- function(x, ..., digits = NULL, verbose = FALSE) {
  cat("fm_mesh_1d object:\n", sep = "")

  cat("  Manifold:\t", x$manifold, "\n", sep = "")
  cat("  #{knots}:\t", length(x$loc), "\n", sep = "")
  cat("  Interval:\t(", paste0(format(x$interval, digits = digits),
    collapse = ", "
  ), ")\n", sep = "")
  clamped <- x$free.clamped & (x$boundary == "free")
  clamped <- c("", " and clamped")[clamped + 1]
  cat("  Boundary:\t(",
    paste0(x$boundary, clamped, collapse = ", "),
    ")\n",
    sep = ""
  )
  cat("  B-spline degree:\t", x$degree, "\n", sep = "")
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")

  invisible(x)
}



#' @export
#' @rdname fmesher-print
print.fm_bbox <- function(x,
                          ...,
                          digits = NULL,
                          verbose = TRUE,
                          newline = TRUE) {
  if (verbose) {
    cat("Bounding box: ", sep = "")
  }
  if (length(x) == 0) {
    cat("NULL")
  } else {
    for (k in seq_along(x)) {
      if (k > 1) {
        cat(" x ", sep = "")
      }
      cat("(",
        paste0(format(x[[k]], digits = digits),
          collapse = ","
        ),
        ")",
        sep = ""
      )
    }
  }
  if (newline) {
    cat("\n", sep = "")
  }
  return(invisible(x))
}


#' @export
#' @rdname fmesher-print
print.fm_tensor <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  ret <-
    c(
      ret,
      list(
        manifold = fm_manifold(x),
        sub_manifolds = vapply(x[["fun_spaces"]], fm_manifold, character(1)),
        sub_dof = vapply(x[["fun_spaces"]], fm_dof, 0L)
      )
    )


  cat("fm_tensor object:\n", sep = "")
  cat("  Manifold:\t", ret$manifold, " = ",
    paste0(ret$sub_manifolds, collapse = " x "),
    "\n",
    sep = ""
  )
  cat("  ", sep = "")
  print(fm_bbox(x), digits = digits)
  cat("  Basis d.o.f.:\t", fm_dof(x), " = ",
    paste0(ret$sub_dof, collapse = " x "),
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @export
#' @rdname fmesher-print
print.fm_collect <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  ret <-
    c(
      ret,
      list(
        manifold = fm_manifold(x),
        size = length(x[["fun_spaces"]]),
        sub_dof = vapply(x[["fun_spaces"]], fm_dof, 0L)
      )
    )


  cat("fm_collect object:\n", sep = "")
  cat("  Manifold:\t", ret$manifold, " x ", ret$size, "\n", sep = "")
  cat("  Bounding boxes:\n", sep = "")
  print(fm_bbox(x), digits = digits, verbose = FALSE)
  cat("  Basis d.o.f.:\t", fm_dof(x), " = ",
    paste0(ret$sub_dof, collapse = " + "),
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @export
#' @rdname fmesher-print
print.fm_lattice_2d <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  ret <-
    c(
      ret,
      list(
        manifold = fm_manifold(x),
        dim = x$dims,
        bbox = fm_bbox(x)
      )
    )

  cat("fm_lattice_2d object:\n", sep = "")
  cat("  Manifold:\t", ret$manifold, "\n", sep = "")
  cat("  Dimensions:\t", paste0(ret$dim, collapse = " x "), "\n", sep = "")
  cat("  ", sep = "")
  print(fm_bbox(x), digits = digits)
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")
  invisible(x)
}

#' @export
#' @rdname fmesher-print
print.fm_lattice_Nd <- function(x, ..., digits = NULL, verbose = FALSE) {
  ret <- list(verbose = verbose)
  ret <-
    c(
      ret,
      list(
        manifold = fm_manifold(x),
        dim = x$dims,
        bbox = fm_bbox(x)
      )
    )

  cat("fm_lattice_Nd object:\n", sep = "")
  cat("  Manifold:\t", ret$manifold, "\n", sep = "")
  cat("  Dimensions:\t", paste0(ret$dim, collapse = " x "), "\n", sep = "")
  cat("  ", sep = "")
  print(fm_bbox(x), digits = digits)
  cat("  Basis d.o.f.:\t", fm_dof(x), "\n", sep = "")
  invisible(x)
}



#' Print method for `fm_basis`
#'
#' Prints information for an [fm_basis] object.
#'
#' @param x [fm_basis()] object
#' @param prefix a prefix to be used for each line. Default is an empty string.
#' @param \dots Unused
#' @returns `invisible(x)`
#' @seealso [fm_basis()]
#' @export
#' @examples
#' print(fm_basis(fmexample$mesh, fmexample$loc, full = TRUE))
print.fm_basis <- function(x, ..., prefix = "") {
  cat(prefix, "fm_basis object\n", sep = "")
  cat(prefix,
    "  Projection matrix (A): ",
    paste0(dim(x$A), collapse = "-by-"),
    "\n",
    sep = ""
  )
  cat(prefix,
    "  Valid evaluations (ok): ",
    sum(x$ok),
    " out of ",
    length(x$ok),
    "\n",
    sep = ""
  )
  cat(
    prefix, "  Additional information: ",
    paste(names(x)[!names(x) %in% c("A", "ok")], collapse = ", "),
    "\n",
    sep = ""
  )

  invisible(x)
}

#' Print method for `fm_evaluator()`
#'
#' Prints information for an [fm_evaluator] object.
#'
#' @param x [fm_evaluator()] object
#' @param \dots Unused
#' @returns `invisible(x)`
#' @seealso [fm_evaluator()]
#' @export
#' @examples
#' print(fm_evaluator(fmexample$mesh, fmexample$loc))
print.fm_evaluator <- function(x, ...) {
  cat("fm_evaluator object\n  proj:\n", sep = "")
  print(fm_basis(x, full = TRUE), prefix = "    ")
  cat(
    "  Additional evaluator information: ",
    paste(names(x)[!names(x) %in% c("proj")], collapse = ", "),
    "\n",
    sep = ""
  )

  invisible(x)
}
