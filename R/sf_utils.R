# helper functions for working with sf objects

# Calculate signed area for polygon
#
# @aliases st_signed_area
# @export
# @param sfg A POLYGON sfg object
# @returns Returns the signed area.  Negative values indicate
# anti-clockwise winding direction.
# @author Andrew Seaton \email{Andrew.Seaton.2@@glasgow.ac.uk}
# @author Finn Lindgren <Finn.Lindgren@@gmail.com>
# @keywords internal

st_signed_area <- function(sfg) {
  warning("st_signed_area is not fully implemented,",
    " and should be avoided until it is.",
    immediate. = TRUE
  )
  if (!inherits(sfg, c("POLYGON", "sfg"))) {
    stop("Signed area only implemented for POLYGON sfg objects")
  }

  coords <- as.matrix(sfg)
  i <- seq_len(nrow(coords) - 1)
  edges <- cbind(coords[i, , drop = FALSE], coords[i + 1, , drop = FALSE])
  area <- sum((edges[, 3] - edges[, 1]) * (edges[, 2] + edges[, 4]) / 2)
  area
}


# Check sfg polygon satisfies standards for POLYGON simple features
#
# @description `r lifecycle::badge("experimental")`
# `sf::st_polygon()` doesn't fully check polygon construction validity.
# For now only implements a basic check for disjoint regions using `st_within()`
#
# @aliases st_check_polygon
# @export
# @param sfg A POLYGON sfg object
# @returns logical; `TRUE` if the `sfg` holes are entirely inside the outer
#   ring, and are disjoint, otherwise `FALSE`. When `FALSE`, the attribute
#   `Message` is set to a character vector describing the detected reasons.
# @keywords internal

st_check_polygon <- function(sfg) {
  if (!inherits(sfg, c("POLYGON", "sfg"))) {
    stop(
      "Requires POLYGON sfg object"
    )
  }

  np <- length(sfg)

  # 1st is outer boundary
  main <- sf::st_polygon(list(as.matrix(sfg[[1]])))
  ok <- TRUE
  msg <- NULL

  # Rest should be holes
  if (length(sfg) > 1) {
    holes <- lapply(
      sfg[-1],
      FUN = as.matrix
    )
    holes <- lapply(
      holes,
      FUN = function(x) sf::st_geometry(sf::st_polygon(list(x)))
    )
    holes <- do.call(
      c, holes
    )

    check_within <- sf::st_within(
      holes, main,
      sparse = FALSE
    )
    if (!all(check_within)) {
      ok <- FALSE
      msg <- c(msg, "Holes not entirely within the outer ring")
    }

    if (length(holes) > 1) {
      check_overlap <- vapply(
        seq_along(holes),
        function(k) {
          !any(sf::st_intersects(
            holes[-k], holes[k],
            sparse = FALSE
          ))
        },
        TRUE
      )
      if (!all(check_overlap)) {
        ok <- FALSE
        msg <- c(msg, "Some holes overlap")
      }
    }
  }

  attr(ok, "Message") <- msg
  ok
}
