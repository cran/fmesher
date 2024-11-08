# fmesher 0.2.0

## New methods

* Add print methods for `fm_basis` and `fm_evaluator` objects (version `0.1.7.9003`)
* Add `fm_manifold_get()` generic method to extract manifold information from general
  objects, so that external objects can implement their own manifold information
  storage (version `0.1.7.9005`)
* Add `fm_crs_plot()` method for plotting `fm_crs` objects with optional
  graticules and Tissot indicatrices (version `0.1.7.9010`)

## New method options

* Add `full` argument to `fm_basis()` to toggle between matrix and full `fm_basis`
  object output (version `0.1.7.9002`)
* Add `loc` plotting option to `geom_fm.fm_mesh_2d` and modify the `ggplot`
  mapping interface for interior and boundary segments in the same method
  (version `0.1.7.9009`)
* Add `format="loc"` argument to `fm_as_sfc.fm_mesh_2d()` for converting mesh
  node coordinates to `sfc_POINT` format (version `0.1.7.9009`)

## Minor updates

* Update documentation and vector coordinate inputs to `fm_lattice_2d()`
  to clarify input interpretation and ensure correct boundary orientation (version `0.1.7.9001`)
* Add some length unit handling to `fm_crs_bounds()` (version `0.1.7.9008`)
* Add control argument `max_batch_size` to `fm_bary.fm_mesh_2d()`, that can be
  supplied via `fm_basis()`, for optional override of the default maximal batch
  calculation size, see #14 (version `0.1.7.9011`)

## Deprecation updates

* Remove `sp` objects from `fmexample` data. Use `fmexample_sp()` to access them if needed
  (version `0.1.7.9004`)
* Move `sp` dependency to Suggests, and remove `inlabru` dependency (version `0.1.7.9006`)
* Further `sp` use protection (version `0.1.7.9007`)
* Increased deprecation warning and error messages for old unsupported methods

# fmesher 0.1.7

* Fix regression bug in `fm_manifold()` that made it ignore all but the first given type options. See #16 (version `0.1.6.9001)
* Fix `plot.fm_mesh_2d` vectorisation bug (version `0.1.6.9002`)
* Add new `fm_subdivide()` method for `fm_mesh_2d` meshes (version `0.1.6.9003`)

# fmesher 0.1.6

* Fix for hiding away-facing triangle edges in `plot.fm_mesh_2d()` and `lines.fm_segm()`.
* Fix duplicated `fm_unify_coords.sf()` method. Thanks to Pablo Paccioretti for debugging and reporting, see #13 (version `0.1.5.9001`)
* Use batched `fm_bary.fm_mesh_2d()` computations for data sizes of `2e5` and upwards.
  This improves performance of `fm_basis()`/`inla.spde.make.A()` for large data sets, see #14 (version `0.1.5.9002`)
* Automatically convert raw `sfc` samplers objects to `sf` objects in `fm_int.list()` (version `0.1.5.9003`)
* Detect and warn about unsupported use of `fm_segm` objects in `fm_int.list()` (version `0.1.5.9004`)
* Add `fm_basis()` and `fm_fem()` support for `fm_tensor` function spaces (version `0.1.5.9005`)
* Add `fm_CRS()` support for `terra` objects (version `0.1.5.9006`)

# fmesher 0.1.5

* Fix bug in `fm_matern_sample()` that prevented the use of a non-NULL `loc` argument.
  For earlier versions, the workaround is to make a separate call to `fm_evaluate()`
  after calling `fm_matern_sample()`.
* Improve speed of `fm_block_log_shift()` by an order of magnitude for multi-block
  cases.
* Fix bug in `plot.fm_mesh_2d` for meshes with 2D coordinate storage (version `0.1.4.9002`)

# fmesher 0.1.4

* Work around `std::get<variant>()` lack of support for MacOS `< 10.14`.

# fmesher 0.1.3

* Fix bug in `fm_rcdt_2d_inla()` that gave different defaults for the `extend` and
  `refine` arguments when equal to `TRUE`, that should be treated the same as `list()`
  (version `0.1.2.9001`)
* Fix bug in `fm_rcdt_2d_inla()` to properly update `n` basis counter element
  when removing unused vertices. (version `0.1.2.9002`)
* Fix bug giving a spurious warning about inconsistent `is.bnd` for empty
  `fm_segm` objects, and inconsistent `grp` vector lengths. (version `0.1.2.9003`)
* Convert some of the old potentially unsafe C++ pointer methods to type safe
  C++17 features

# fmesher 0.1.2

* Various documentation improvements, in particular for INLA compatibility
* Modify mesh refinement tests to directly check the refinement criteria
  instead of the specific mesh result, to account for differences in
  floating point behaviour on M1 processors.
* Modify tests of non-mesh-generation features to use precomputed meshes
  or meshes with stable properties
* Protect against invalid `tv` inputs
* Revert from `\text{}` to `\textrm{}`, as AMS extensions are only supported
  from R 4.2.2 (https://www.stats.bris.ac.uk/R/doc/manuals/r-devel/R-exts.pdf
  2023-08-24, page 90), and CRAN oldrel for macOS is 4.2.0, not 4.2.3

# fmesher 0.1.1

* Simplify LICENSE information in the R package DESCRIPTION
* Fix example for legacy `plot_PolySet` method
* Fix C++ unused-variable warning detected by clang-tidy version 16.0.6
* Add workaround for `sf::st_buffer` not supporting negative distance
  for `s2` with longlat coordinates. Fixes #5

# fmesher 0.1.0

* Full conversion of the fmesher interface from the standalone binary in R-INLA
  (https://www.r-inla.org) to an interface powered by `Rcpp`
  (https://cran.r-project.org/package=Rcpp)
* Tools for 2D and 1D function spaces, see `fm_mesh_2d()` and `fm_mesh_1d()`
* Tools for CRS handling see `fm_crs()` and `fm_transform()`
* Plotting support for base graphics, `rgl`, and `ggplot2`
* Added basic GMRF tools, see `fm_matern_precision()`

# fmesher 0.0.9

* Basic fmesher library I/O interface
* Added a `NEWS.md` file to track changes to the package.
