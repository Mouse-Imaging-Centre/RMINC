library(testthat)

context("mincGetTagFile")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

matrix_3x4 <- mincGetTagFile(file.path(dataPath, "3rows_4cols.tag"))
matrix_4x3 <- mincGetTagFile(file.path(dataPath, "4rows_3cols_2comments.tag"))

test_that("mincGetTagFile works", {
  expect_equal(matrix_3x4, matrix(c(0.15, 0.175, 0.175, 10.15, 10.20, 10.20, 0.35, 0.325, 0.35 , 0, 0, 1), nrow=3, ncol=4))
  expect_equal(matrix_4x3, matrix(c(0.15, 0.175, 0.175, 0.2, 10.15, 10.2, 10.2, 0.25, 0.35, 0.325, 0.35, 5.75), nrow=4, ncol=3))
})

context("File Limits")

test_that("File limit check works", {
  expect_error(RMINC:::enoughAvailableFileDescriptors(Inf))
  expect_true(RMINC:::enoughAvailableFileDescriptors(0))
  expect_false(RMINC:::enoughAvailableFileDescriptors(Inf, error = FALSE))
})

test_that("File checking code all fails", {
  cur_ulim <- RMINC:::checkCurrentUlimit()
  skip_if_not(is.finite(cur_ulim))
  
  frame <- 
    data_frame(
      files = rep("", cur_ulim + 1)
      , group = rep(c("a", "b"), length.out = cur_ulim + 1))
  
  expect_error(mincLmer(files ~ (1 | group), mask = frame$files[1], data = frame), "file descriptors")
  expect_error(mincMean(frame$files), "file descriptors")
  expect_error(anatGetAll(frame$files, method = "labels"), "file descriptors")
  expect_error(vertexLm(files ~ group, data = frame), "file descriptors")
})