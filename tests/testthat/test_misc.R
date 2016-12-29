library(testthat)

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
      , group = c("a", "b"))
  
  expect_error(mincLmer(files ~ (1 | group), data = frame))
  expect_error(mincMean(frame$files))
  expect_error(anatGetAll(frame$files, method = "labels"))
  expect_error(vertexLm(files ~ group, data = frame))
})