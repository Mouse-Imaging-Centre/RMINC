library(testthat)

context("File Limits")

test_that("File limit check works", {
  expect_error(RMINC:::enoughAvailableFileDescriptors(Inf))
  expect_true(RMINC:::enoughAvailableFileDescriptors(0))
  expect_false(RMINC:::enoughAvailableFileDescriptors(Inf, error = FALSE))
})