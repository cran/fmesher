# sp sf test
test_that("fm_cprod(..., na.rm = TRUE) sf output can be generated", {
  sf_obj1 <- sf::st_as_sf(data.frame(x = 1:3, y = 3:5), coords = c("x", "y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x = 3:6, y = 5:8), coords = c("x", "y"))
  ips <- fm_cprod(sf_obj1, sf_obj2, na.rm = TRUE)

  expect_s3_class(ips, "sf")
  expect_equal(nrow(ips), 1)
  expect_equal(names(ips), c("geometry", "weight", ".block", ".block_origin"))
  expect_equal(as.numeric(unlist(sf::st_geometry(ips))), c(3, 5))
  expect_equal(as.numeric(ips[["weight"]]), 1)
})

test_that("fm_cprod(..., na.rm = FALSE) sf output can be generated", {
  sf_obj1 <- sf::st_as_sf(data.frame(x = 1:3, y = 3:5), coords = c("x", "y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x = 3:6, y = 5:8), coords = c("x", "y"))
  ips <- fm_cprod(sf_obj1, sf_obj2, na.rm = FALSE)

  expect_s3_class(ips, "sf")
  expect_equal(nrow(ips), 6)
  expect_equal(names(ips), c("geometry", "weight", ".block", ".block_origin"))
  expect_equal(
    as.numeric(unlist(sf::st_geometry(ips))),
    c(3, 5, 4, 6, 5, 7, 6, 8, 1, 3, 2, 4)
  )
  expect_equal(as.numeric(ips[["weight"]]), rep(c(1, NA), c(1, 5)))
})

test_that("fm_cprod(..., na.rm = FALSE) sf output with different geometry", {
  sf_obj1 <- sf::st_as_sf(data.frame(x = 1:3, y = 3:5),
    coords = c("x", "y")
  )
  sf_obj2 <- sf::st_as_sf(data.frame(x = 3:6, y = 5:8),
    coords = c("x", "y")
  )
  sf_obj1 <- sf::st_as_sf(tibble::tibble(geometry1 = sf_obj1$geometry))
  sf_obj2 <- sf::st_as_sf(tibble::tibble(geometry2 = sf_obj2$geometry))
  ips <- fm_cprod(sf_obj1, sf_obj2, na.rm = FALSE)

  expect_s3_class(ips, "sf")
  expect_equal(nrow(ips), 12)
  expect_equal(
    sort(names(ips)),
    sort(c(
      "geometry1", "geometry2", "weight",
      ".block", ".block_origin"
    ))
  )
  expect_equal(
    as.numeric(unlist(sf::st_geometry(ips))),
    unlist(
      rep(list(c(1, 3), c(2, 4), c(3, 5)), times = 4)
    )
  )
  expect_equal(
    as.numeric(unlist(ips$geometry2)),
    c(
      rep(c(3, 5), 3),
      rep(c(4, 6), 3),
      rep(c(5, 7), 3),
      rep(c(6, 8), 3)
    )
  )
  expect_equal(as.numeric(ips[["weight"]]), rep(1, 12))
})

test_that("fm_cprod(na.rm = TRUE) sp output can be generated", {
  skip_if_not(fm_safe_sp())

  sf_obj1 <- sf::st_as_sf(data.frame(x = 1:3, y = 3:5), coords = c("x", "y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x = 3:6, y = 5:8), coords = c("x", "y"))
  if (require(sp, quietly = TRUE)) {
    sp_obj1 <- as(sf_obj1, "Spatial")
    sp_obj2 <- as(sf_obj2, "Spatial")
  }
  ips <- fm_cprod(sp_obj1, sp_obj2, na.rm = TRUE)

  expect_s4_class(ips, "Spatial") # or precisely SpatialPointsDataFrame
  expect_equal(nrow(ips), 1)
  expect_equal(names(ips), c("weight", ".block", ".block_origin"))
  expect_equal(as.numeric(sp::coordinates(ips)), c(3, 5))
  expect_equal(as.numeric(ips[["weight"]]), 1)
})

test_that("fm_cprod(na.rm = FALSE) sp output can be generated", {
  skip_if_not(fm_safe_sp())
  sf_obj1 <- sf::st_as_sf(data.frame(x = 1:3, y = 3:5), coords = c("x", "y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x = 3:6, y = 5:8), coords = c("x", "y"))
  if (require(sp, quietly = TRUE)) {
    sp_obj1 <- as(sf_obj1, "Spatial")
    sp_obj2 <- as(sf_obj2, "Spatial")
  }
  ips <- fm_cprod(sp_obj1, sp_obj2, na.rm = FALSE)

  expect_s4_class(ips, "Spatial") # or precisely SpatialPointsDataFrame
  expect_equal(nrow(ips), 6)
  expect_equal(names(ips), c("weight", ".block", ".block_origin"))
  expect_equal(as.numeric(sp::coordinates(ips)), c(3:6, 1:2, 5:8, 3:4))
  expect_equal(as.numeric(ips[["weight"]]), rep(c(1, NA), c(1, 5)))
})
