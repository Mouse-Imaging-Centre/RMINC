library(testthat)
library(dplyr)

context("visualizations")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")
likeVol <- file.path(dataPath, "testsummaryminc1.mnc")

uniformData <-
  runif(3375) %>%
  `dim<-`(c(15,15,15)) %>%
  mincArray

discreteData <-
  sample(1:10, 3375, replace = TRUE) %>%
  `dim<-`(c(15,15,15)) %>%
  mincArray


test_that("Slice Extraction Works", {
  expect_equal(RMINC:::getSlice(uniformData, 5, 1)$s, uniformData[5,,])
  expect_equal(RMINC:::getSlice(uniformData, 2, 2)$s, uniformData[,2,])
  expect_equal(RMINC:::getSlice(uniformData, 7, 3)$s, uniformData[,,7])
})
  
uniformSlice <- RMINC:::getSlice(uniformData, 13, 1)$s
discreteSlice <- RMINC:::getSlice(discreteData, 10, 2)$s

test_that("Slice Shifting/Clamping Works", {
  unifClamp <- RMINC:::scaleSlice(uniformSlice, low = .2, high = .7, 
                                  underTransparent = FALSE)
  
  expect_true(all(unifClamp <= abs(.7 - .2)))
  expect_true(all(unifClamp >= 0))
  
  discreteClamp <- RMINC:::scaleSlice(discreteSlice, low = 3, high = 4, 
                                      underTransparent = TRUE)
  expect_true(all(discreteClamp > 0, na.rm = TRUE))
  expect_true(all(discreteClamp <= 1, na.rm = TRUE))
  expect_true(any(discreteClamp == 1, na.rm = TRUE))
})

test_that("Slice Colour Scaling Works", {
  scaledUnif <- 
    RMINC:::scaleSliceToPalette(
      RMINC:::scaleSlice(uniformSlice, low = .5, high = .7, underTransparent = FALSE),
      low = .5,
      high = .7,
      palette = heat.colors(300)
    )
  
  expect_true(max(scaledUnif) <= 300)
  expect_true(min(scaledUnif) >= 0)
  
  scaledDiscrete <- 
    RMINC:::scaleSliceToPalette(
      RMINC:::scaleSlice(discreteSlice, low = 3, high = 4, underTransparent = TRUE),
      low = 3,
      high = 4,
      palette = heat.colors(300)
    )
  
  expect_true(all(scaledDiscrete <= 300, na.rm = TRUE))
  expect_true(any(scaledDiscrete == 300, na.rm = TRUE))
  expect_true(all(scaledDiscrete > 0, na.rm = TRUE))
})