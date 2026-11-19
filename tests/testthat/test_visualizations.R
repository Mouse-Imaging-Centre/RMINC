library(testthat)
library(dplyr)

context("visualizations")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")
likeVol <- file.path(dataPath, "testsummaryminc1.mnc")

uniformData <-
  runif(3375) %>%
  `dim<-`(c(15, 15, 15)) %>%
  mincArray

discreteData <-
  sample(1:10, 3375, replace = TRUE) %>%
  `dim<-`(c(15, 15, 15)) %>%
  mincArray


test_that("Slice Extraction Works", {
  expect_equal(RMINC:::getSlice(uniformData, 5, 1)$s, uniformData[5, , ])
  expect_equal(RMINC:::getSlice(uniformData, 2, 2)$s, uniformData[, 2, ])
  expect_equal(RMINC:::getSlice(uniformData, 7, 3)$s, uniformData[,, 7])
})

uniformSlice <- RMINC:::getSlice(uniformData, 13, 1)$s
discreteSlice <- RMINC:::getSlice(discreteData, 10, 2)$s

test_that("Slice Shifting/Clamping Works", {
  unifClamp <- RMINC:::scaleSlice(
    uniformSlice,
    low = .2,
    high = .7,
    underTransparent = FALSE
  )

  expect_true(all(unifClamp <= abs(.7 - .2)))
  expect_true(all(unifClamp >= 0))

  discreteClamp <- RMINC:::scaleSlice(
    discreteSlice,
    low = 3,
    high = 4,
    underTransparent = TRUE
  )
  expect_true(all(discreteClamp > 0, na.rm = TRUE))
  expect_true(all(discreteClamp <= 1, na.rm = TRUE))
  expect_true(any(discreteClamp == 1, na.rm = TRUE))
})

test_that("Slice Colour Scaling Works", {
  scaledUnif <-
    RMINC:::scaleSliceToPalette(
      RMINC:::scaleSlice(
        uniformSlice,
        low = .5,
        high = .7,
        underTransparent = FALSE
      ),
      low = .5,
      high = .7,
      palette = heat.colors(300)
    )

  expect_true(max(scaledUnif) <= 300)
  expect_true(min(scaledUnif) >= 0)

  scaledDiscrete <-
    RMINC:::scaleSliceToPalette(
      RMINC:::scaleSlice(
        discreteSlice,
        low = 3,
        high = 4,
        underTransparent = TRUE
      ),
      low = 3,
      high = 4,
      palette = heat.colors(300)
    )

  expect_true(all(scaledDiscrete <= 300, na.rm = TRUE))
  expect_true(any(scaledDiscrete == 300, na.rm = TRUE))
  expect_true(all(scaledDiscrete > 0, na.rm = TRUE))
})

signedData <-
  runif(3375, -5, 5) %>%
  `dim<-`(c(15, 15, 15)) %>%
  mincArray

# Run mincPlotSliceSeries on a null device and record each call to
# plotrix::color.legend: its labels and the figure region it was drawn in
recordSliceSeriesLegends <- function(...) {
  calls <- list()
  testthat::local_mocked_bindings(
    color.legend = function(xl, yb, xr, yt, legend, rect.col, ...) {
      calls[[length(calls) + 1]] <<- list(legend = legend, fig = par("fig"))
      invisible(NULL)
    },
    .package = "plotrix"
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  mincPlotSliceSeries(...)
  calls
}

test_that("the symmetric legend shows the computed range (#130)", {
  skip_if_not_installed("plotrix")
  legends <- recordSliceSeriesLegends(
    uniformData,
    signedData,
    mfrow = c(2, 2),
    low = 0.2, # high is not given, so it comes from the histogram
    symmetric = TRUE,
    legend = "t-statistics"
  )
  expect_length(legends, 2)
  negative <- legends[[1]]$legend
  positive <- legends[[2]]$legend
  expect_length(positive, 2)
  expect_length(negative, 2)
  expect_equal(positive[1], 0.2)
  expect_gt(positive[2], 0.2)
  # the negative scale goes from -high (bottom) to -low (top)
  expect_equal(negative, -rev(positive))
})

test_that("the legend fills the side column without a locator (#197)", {
  skip_if_not_installed("plotrix")
  no_locator <- recordSliceSeriesLegends(
    uniformData,
    uniformData,
    mfrow = c(2, 2),
    low = 0.2,
    high = 0.8,
    legend = "values",
    locator = FALSE
  )
  expect_length(no_locator, 1)
  fig <- no_locator[[1]]$fig
  expect_equal(fig[4] - fig[3], 1) # full height
  expect_equal(fig[1], 2 / 3) # in the column to the right of the slices

  with_locator <- recordSliceSeriesLegends(
    uniformData,
    uniformData,
    mfrow = c(2, 2),
    low = 0.2,
    high = 0.8,
    legend = "values",
    locator = TRUE
  )
  fig <- with_locator[[1]]$fig
  expect_equal(fig[4] - fig[3], 1 / 2) # below the locator
  expect_equal(fig[1], 2 / 3)
})

test_that("slices are drawn with the voxel spacing (#216)", {
  vol <- uniformData
  attr(vol, "likeVolume") <- likeVol
  # file order (z, y, x); the mincArray order is (x, y, z)
  local_mocked_bindings(minc.separation.sizes = function(filename) c(3, 2, 1))

  expect_equal(RMINC:::getSlice(vol, 5, 1)$asp, 3 / 2)
  expect_equal(RMINC:::getSlice(vol, 5, 2)$asp, 3 / 1)
  expect_equal(RMINC:::getSlice(vol, 5, 3)$asp, 2 / 1)

  # one voxel in the plot is 3 times as high as it is wide
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  mincImage(vol, dimension = 2, slice = 5)
  usr <- par("usr")
  pin <- par("pin")
  expect_equal(
    (pin[2] / diff(usr[3:4])) / (pin[1] / diff(usr[1:2])),
    3,
    tolerance = 1e-6
  )

  # no file, or the option is off: square voxels
  expect_equal(RMINC:::getSlice(uniformData, 5, 2)$asp, 1)
  withr::local_options(RMINC_use_voxel_spacing = FALSE)
  expect_equal(RMINC:::getSlice(vol, 5, 2)$asp, 1)
})

test_that("slice ranges outside the volume are handled (#220)", {
  # the shiny slider starts at 55-280 of 340, before it knows the volume
  expect_equal(RMINC:::clampSliceRange(c(55, 280), 193), c(55, 193))
  expect_equal(RMINC:::clampSliceRange(c(0, 10), 15), c(1, 10))
  expect_equal(RMINC:::clampSliceRange(c(12, 3), 15), c(3, 12))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(
    mincPlotSliceSeries(uniformData, begin = 5, end = 20),
    "begin and end must give slices from 1 to 15"
  )
  expect_error(
    mincImage(uniformData, dimension = 2, slice = 16),
    "slice must be a number from 1 to 15"
  )
})

test_that("mesh colours and opacity are given per vertex (#274)", {
  skip_if_not_installed("rgl")
  withr::local_options(rgl.useNULL = TRUE)
  # two triangles that share the edge from vertex 2 to vertex 3
  bic <- list(
    vertex_matrix = matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0), nrow = 3),
    triangle_matrix = matrix(c(1, 2, 3, 2, 4, 3), nrow = 3)
  )
  palette <- c("red", "green", "blue", "white")
  mesh <-
    create_mesh(bic) %>%
    colour_mesh(
      colour_map = 1:4,
      colour_range = c(1, 4),
      symmetric = FALSE,
      palette = palette
    ) %>%
    add_opacity(a_map = 1:4, a_range = c(0.25, 1))

  expect_equal(mesh$material$color, palette)
  expect_equal(mesh$material$alpha, c(0.25, 0.5, 0.75, 1))

  rgl::open3d()
  on.exit(rgl::close3d(), add = TRUE)
  id <- rgl::shade3d(mesh, override = FALSE)
  drawn <- rgl::rgl.attrib(id, "colors")
  # rgl draws the triangles' corners in the order of mesh$it
  expected <- t(col2rgb(palette[mesh$it])) / 255
  expect_equal(unname(drawn[, 1:3]), unname(expected))
  # rgl stores alpha in steps of 1/255
  expect_equal(
    unname(drawn[, 4]),
    c(0.25, 0.5, 0.75, 1)[mesh$it],
    tolerance = 0.01
  )
})
