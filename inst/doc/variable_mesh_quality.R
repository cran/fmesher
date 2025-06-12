## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  collapse = TRUE,
  comment = "#>",
  dev = "png",
  dev.args = list(type = "cairo-png"),
  fig.width = 5,
  fig.height = 5 * 5 / 7
)
options(bitmapType = "cairo")

## ----label1, echo=TRUE, message=FALSE-----------------------------------------
library(ggplot2)
library(fmesher)

loc <- as.matrix(expand.grid(1:10, 1:10))
bnd_sf <- fm_nonconvex_hull(loc, convex = 1, concave = 10)
bnd <- fm_as_segm(bnd_sf)
loc <- fm_hexagon_lattice(bnd_sf, edge_len = 1)
mesh1 <- fm_rcdt_2d_inla(
  loc = loc,
  boundary = bnd,
  refine = list(max.edge = Inf)
)
ggplot() +
  geom_fm(data = mesh1)

## ----echo=TRUE, warning=FALSE,message=FALSE-----------------------------------
mesh2 <- fm_rcdt_2d_inla(
  loc = loc,
  boundary = bnd,
  refine = list(max.edge = 0.5)
)
ggplot() +
  geom_fm(data = mesh2) +
  coord_sf(default = TRUE)

## ----echo=TRUE, warning=FALSE,message=FALSE-----------------------------------
qual_loc <- function(loc) {
  if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
  }
  pmax(0.05, (loc[, 1] * 2 + loc[, 2]) / 16)
}
mesh3 <- fm_rcdt_2d_inla(
  loc = loc,
  boundary = bnd,
  refine = list(max.edge = Inf),
  quality.spec = list(
    loc = qual_loc(loc),
    segm = qual_loc(bnd$loc)
  )
)
ggplot() +
  geom_fm(data = mesh3) +
  coord_sf(default = TRUE)

## ----echo=TRUE, message=FALSE, size="small",warning=FALSE---------------------
qual_bnd <- function(loc) {
  rep(Inf, nrow(loc))
}
mesh5 <- fm_rcdt_2d_inla(
  loc = loc,
  boundary = bnd,
  refine = list(max.edge = Inf),
  quality.spec = list(
    loc = qual_loc(loc),
    segm = qual_bnd(bnd$loc)
  )
)
ggplot() +
  geom_fm(data = mesh5) +
  coord_sf(default = TRUE)

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
qual_bnd <- function(loc) {
  if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
  }
  pmax(0.1, 1 - abs(loc[, 2] / 10)^2)
}
mesh5 <- fm_rcdt_2d_inla(
  loc = loc,
  boundary = bnd,
  refine = list(
    max.edge = Inf,
    max.n.strict = 5000
  ),
  quality.spec = list(
    loc = qual_loc(loc),
    segm = qual_bnd(bnd$loc)
  )
)
ggplot() +
  geom_fm(data = mesh5) +
  coord_sf(default = TRUE)

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
out <- fm_assess(mesh5,
  spatial.range = 5,
  alpha = 2,
  dims = c(200, 200)
)
print(names(out))

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
ggplot() +
  geom_tile(
    data = out,
    aes(geometry = geometry, fill = edge.len),
    stat = "sf_coordinates"
  ) +
  coord_sf(default = TRUE) +
  scale_fill_distiller(
    type = "seq",
    direction = 1,
    palette = 1,
    na.value = "transparent"
  )

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
sd.dev.limits <- 1 + c(-1, 1) * max(abs(range(out$sd.dev, na.rm = TRUE) - 1))
col.values <- 2 * seq(0, 1, length.out = 100) - 1
col.values <- (sign(col.values) * abs(col.values)^1.5 + 1) / 2
ggplot() +
  geom_tile(
    data = out,
    aes(geometry = geometry, fill = sd.dev),
    stat = "sf_coordinates"
  ) +
  coord_sf(default = TRUE) +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    na.value = "transparent",
    limits = sd.dev.limits,
    values = col.values
  )

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
mesh6 <- fm_rcdt_2d_inla(
  loc = loc,
  boundary = bnd,
  refine = list(
    min.angle = 30,
    max.edge = Inf,
    max.n.strict = 5000
  ),
  quality.spec = list(
    loc = qual_loc(loc),
    segm = qual_bnd(bnd$loc)
  )
)
out6 <- fm_assess(mesh6,
  spatial.range = 5,
  alpha = 2,
  dims = c(200, 200)
)

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
ggplot() +
  geom_tile(
    data = out6,
    aes(geometry = geometry, fill = sd.dev),
    stat = "sf_coordinates"
  ) +
  coord_sf(default = TRUE) +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    na.value = "transparent",
    limits = sd.dev.limits,
    values = col.values
  )

## ----echo=TRUE, message=FALSE, size="small"-----------------------------------
mesh5
mesh6

