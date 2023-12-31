test_that("fm_pixels sp vs sf", {
  skip_on_cran()
  mesh <- fm_mesh_2d_inla(
    cbind(0, 0),
    offset = 10,
    max.edge = 1,
    crs = fm_CRS("longlat_globe")
  )

  mydata <- sp::SpatialPointsDataFrame(
    mesh$loc,
    data = data.frame(y = rnorm(mesh$n) + 10),
    proj4string = fm_CRS("longlat_globe")
  )

  system.time({
    surface1 <- fm_pixels(mesh, dims = c(5, 5), mask = TRUE, format = "sp")
  })

  system.time({
    surface2 <- fm_pixels(mesh, dims = c(5, 5), mask = TRUE, format = "sf")
  })

  expect_equal(
    sf::st_coordinates(sf::st_as_sf(surface1)),
    sf::st_coordinates(surface2)
  )
})
