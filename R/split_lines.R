#' @include deprecated.R

# fm_split_lines ####

#' @title Split lines at triangle edges
#'
#' @description Compute intersections between line segments and triangle edges,
#' and filter out segment of length zero.
#'
#' @param mesh An [fm_mesh_2d] object
#' @param segm An [fm_segm()] object with segments to be split
#' @param ... Unused.
#' @returns An [fm_segm()] object with the same crs as the mesh,
#' with an added field `origin`, that for each new segment gives the
#' originator index into to original `segm` object for each new line segment.
#'
#' @author Finn Lindgren <Finn.Lindgren@@gmail.com>
#'
#' @export
#' @examples
#' mesh <- fm_mesh_2d(
#'   boundary = fm_segm(
#'     rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
#'     is.bnd = TRUE
#'   )
#' )
#' splitter <- fm_segm(rbind(c(0.8, 0.2), c(0.2, 0.8)))
#' segm_split <- fm_split_lines(mesh, splitter)
#'
#' plot(mesh)
#' lines(splitter)
#' points(segm_split$loc)
fm_split_lines <- function(mesh, ...) {
  UseMethod("fm_split_lines")
}


#' @rdname fm_split_lines
#' @export
fm_split_lines.fm_mesh_2d <- function(mesh, segm, ...) {
  segm <-
    fm_transform(
      fm_as_segm(segm),
      crs = fm_crs(mesh),
      crs0 = fm_crs(segm),
      passthrough = TRUE
    )
  origin <- seq_len(NROW(segm$idx))
  if (NROW(segm$loc) > 0) {
    # Filter out segments not on the mesh
    t1 <- fm_bary(mesh, loc = segm$loc, crs = fm_crs(segm))$index
    keep <- !(is.na(t1[segm$idx[, 1]]) | is.na(t1[segm$idx[, 2]]))
    # if (any(!keep)) { warning("points outside boundary! filtering...")}
    segm <- fm_segm(
      loc = segm$loc,
      idx = segm$idx[keep, , drop = FALSE],
      grp = segm$grp[keep],
      is.bnd = all(segm$is.bnd),
      crs = fm_crs(segm)
    )
    origin <- origin[keep]
  }

  if (NROW(segm$idx) == 0) {
    segm$origin <- origin
    return(segm)
  }

  # Split the segments into parts
  if (NCOL(segm$loc) == 2) {
    segm$loc <- cbind(segm$loc, rep(0, NROW(segm$loc)))
  }
  splt <- fmesher_split_lines(
    mesh_loc = mesh$loc,
    mesh_tv = mesh$graph$tv - 1L,
    loc = segm$loc,
    idx = segm$idx - 1L,
    options = list()
  )

  stopifnot(!any(is.na(splt$split.loc)))

  indexoutput <- list("split.idx", "split.t", "split.origin")
  for (name in intersect(names(splt), indexoutput)) {
    splt[[name]] <- splt[[name]] + 1L
  }

  segm.split <- fm_segm(
    loc = splt$split.loc,
    idx = splt$split.idx,
    grp = segm$grp[splt$split.origin],
    is.bnd = all(segm$is.bnd),
    crs = fm_crs(segm)
  )
  origin <- origin[splt$split.origin]

  #  plot(mesh)
  #  lines(segm)
  #  points(segm.split$loc[, 1:2], col="red", pch = 20)
  #  points(segm$loc[, 1:2], col="blue", pch = 20)

  # Filter out zero length segments
  keep <- rowSums((segm.split$loc[segm.split$idx[, 2], , drop = FALSE] -
    segm.split$loc[segm.split$idx[, 1], , drop = FALSE])^2) > 0
  segm.split <- fm_segm(
    loc = segm.split$loc,
    idx = segm.split$idx[keep, , drop = FALSE],
    grp = segm.split$grp[keep],
    is.bnd = all(segm.split$is.bnd),
    crs = fm_crs(segm)
  )
  segm.split$origin <- origin[which(keep)]

  segm.split
}
