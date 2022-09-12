library(testthat)
context("mincTable")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "RMINC-test-data-main/rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"), stringsAsFactors = TRUE)
maskfile <- file.path(dataPath, "testminc-mask.mnc")

ref <- sapply(gf$jacobians_fixed, function(f) mincGetVolume(f))
ref_masked <- ref[mincGetVolume(maskfile) > 0.5,]

mt <- mincTable(gf$jacobians_fixed)
mt_masked <- mincTable(gf$jacobians_fixed, mask = maskfile)

mt_back <- mincTable(gf$jacobians_fixed, file_backed = TRUE)
mt_back_masked <- mincTable(gf$jacobians_fixed, mask = maskfile, file_backed = TRUE)

test_that("Unmasked mincTable works", {
  expect_equal(ref, mt)
  expect_equal(ref, mt_back[,])
})

test_that("Masked mincTable works", {
  expect_equal(ref_masked, mt_masked)
  expect_equal(ref_masked, mt_back_masked[,])
})
