## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
suppressPackageStartupMessages(library(fmesher))
suppressPackageStartupMessages(library(tibble))

convert_fun_link <- function(x, webref) {
  pkg_name <- sub(
    pattern = "^(([^:]*)::)?([^(]+)(\\(.*\\))(<-)?$",
    replacement = "\\2",
    x = x
  )
  if (identical(pkg_name, "")) {
    pkg_name <- "fmesher"
  }
  nm <- sub(
    pattern = "^(([^:]*)::)?([^(]+)(\\(.*\\))(<-)?$",
    replacement = "\\3",
    x = x
  )
  args <- sub(
    pattern = "^(([^:]*)::)?([^(]+)(\\(.*\\))(<-)?$",
    replacement = "\\4\\5",
    x = x
  )
  help_name <- tryCatch(
    {
      basename(utils::help((nm), package = (pkg_name)))
    },
    error = function(e) {
      nm
    }
  )
  paste0(
    '<a href="', webref[[pkg_name]], help_name, '.html">`',
    if (!identical(pkg_name, "fmesher")) {
      paste0(pkg_name, "::")
    } else {
      ""
    },
    nm,
    args,
    "`</a>"
  )
}

convert_fun_links <- function(df, packages = colnames(df)) {
  webref <- list(
    fmesher = "https://inlabru-org.github.io/fmesher/reference/",
    inlabru = "https://inlabru-org.github.io/inlabru/reference/"
  )
  packages <- setdiff(packages, "Comments")
  for (k in seq_len(nrow(df))) {
    for (pkg in packages) {
      if (identical(pkg, "INLA")) {
        df[["INLA"]][[k]] <- paste0("`", df$INLA[[k]], "`")
      } else {
        df[[pkg]][[k]] <-
          vapply(
            df[[pkg]][[k]],
            function(x) convert_fun_link(x, webref),
            ""
          )
      }
    }
  }
  df
}

## ----echo = FALSE-------------------------------------------------------------
df <- tribble(
  ~INLA, ~fmesher,
  "inla.mesh.create()", c("fm_rcdt_2d()", "fm_rcdt_2d_inla()"),
  "inla.mesh.2d()", c("fm_mesh_2d()", "fm_mesh_2d_inla()"),
  "inla.delaunay()", "fm_delaunay_2d()",
  "inla.mesh.1d()", "fm_mesh_1d()",
  "inla.mesh.lattice()", "fm_lattice_2d()",
  "inla.mesh.segment()", "fm_segm()",
  "inla.nonconvex.hull()", c(
    "fm_nonconvex_hull()",
    "fm_extensions()",
    "fm_simplify()"
  ),
  c(
    "inla.nonconvex.hull()",
    "inla.contour.segment()",
    "inla.simplify.curve()"
  ),
  c(
    "fm_nonconvex_hull_inla()",
    "fm_simplify_helper()",
    "fm_segm_contour_helper()"
  )
)
df <- convert_fun_links(df)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(df)

## ----echo = FALSE-------------------------------------------------------------
df <- tribble(
  ~INLA, ~fmesher, ~inlabru,
  "inla.mesh.projector()", "fm_evaluator()", character(0),
  "inla.mesh.project()", "fm_evaluate()", character(0),
  "inla.spde.make.A()", c(
    "fm_basis()",
    "fm_row_kron()",
    "fm_block()",
    "fm_block_eval()"
  ), c(
    "inlabru::bru_mapper_multi()",
    "inlabru::ibm_jacobian()",
    "inlabru::bru_mapper_aggregate()"
  ),
  "inla.mesh.deriv()", "fm_basis()", character(0)
)
df <- convert_fun_links(df)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(df)

## ----echo = FALSE-------------------------------------------------------------
df <- tribble(
  ~INLA, ~fmesher, ~Comments,
  c("inla.mesh.fem()", "inla.mesh.1d.fem()"), "fm_fem()", " ",
  " ", "fm_matern_precision()", " ",
  " ", "fm_matern_sample()", paste0(
    "Convenience function that combines",
    " `fm_matern_precision()`",
    " and `fm_sample()`."
  ),
  " ", "fm_covariance()",
  paste0(
    "Basic helper function for computing",
    " covariances between different locations.",
    " Can produce sparse inverses like",
    " `inla.qinv()`, but currently (version 0.1.1) only",
    " by a 'brute force' method."
  ),
  "inla.qinv", "fm_qinv()",
  paste0(
    "Produce sparse inverses like",
    " `inla.qinv()`, but currently (version 0.2.0.9010)",
    " by an R implementation."
  ),
  " ", "fm_sample()", "Basic sampling method, like `inla.qsample()`"
)
df <- convert_fun_links(df)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(df)

## ----echo = FALSE-------------------------------------------------------------
df <- tribble(
  ~INLA,
  ~fmesher,
  ~Comments,
  "summary.inla.mesh()",
  c(
    "print.fm_mesh_2d()",
    "print.fm_segm()",
    "print.fm_mesh_1d()"
  ),
  "Use `print(mesh)` etc."
)
df <- convert_fun_links(df)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(df)

## ----echo = FALSE-------------------------------------------------------------
df <-
  tibble::tribble(
    ~INLA, ~fmesher, ~Comments,
    "inla.spTransform()", "fm_transform()", "",
    "mesh$crs", c("fm_crs(mesh)", "fm_CRS(mesh)"),
    paste0(
      "The crs may now be stored in different formats;",
      " use `fm_crs()` for `sf` format, and `fm_CRS()` for `sp` format.",
      " `fmesher` will attempt to convert when needed."
    ),
    "mesh$crs<-", "fm_crs(mesh)<-",
    paste0(
      "Direct assignment of crs information should be avoided,",
      " but is allowed as long as its compatible with the actual",
      " mesh coordinates."
    )
  )
df <- convert_fun_links(df)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(df)

## ----echo = FALSE-------------------------------------------------------------
df <- tribble(
  ~INLA, ~fmesher, ~inlabru, ~Comments,
  "No ggplot support",
  c("geom_fm(data = mesh)", "geom_fm(data = segm)"),
  "inlabru::gg(mesh)",
  "Use `ggplot() + geom_fm(data = mesh)` and `inlabru::gg()` methods",
  "plot.inla.mesh(rgl = FALSE)",
  c("plot.fm_mesh_2d()", "lines.fm_mesh_2d()"),
  character(0),
  "Use `plot()` or `lines()`",
  "lines.inla.mesh.segment(rgl = FALSE)",
  c(
    "plot.fm_segm()",
    "lines.fm_segm()"
  ),
  character(0),
  "Use `plot()` or `lines()`",
  "plot.inla.mesh(rgl = TRUE)",
  c("plot_rgl()", "lines_rgl()"),
  character(0),
  "",
  "lines.inla.mesh.segment(rgl = TRUE)",
  c(
    "plot_rgl()",
    "lines_rgl()"
  ),
  character(0),
  ""
)
df <- convert_fun_links(df)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(df)

