library(testthat)
context("mincWriteVolume")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(
  file.path(dataPath, "minc_summary_test_data.csv"),
  stringsAsFactors = TRUE
)
like_file <- as.character(gf$jacobians_0.2[1])
n_vox <- prod(minc.dimensions.sizes(like_file))

write_and_read <- function(values) {
  out <- tempfile(fileext = ".mnc")
  on.exit(unlink(out))
  mincWriteVolume(values, out, like.filename = like_file, clobber = TRUE)
  as.numeric(mincGetVolume(out))
}

test_that("written values keep their precision (#78)", {
  # A wide range with fractional values. When the output was unsigned
  # short, the error was about range / 65535 (here about 0.3).
  values <- seq(-1e4, 1e4, length.out = n_vox) + 0.123
  back <- write_and_read(values)
  expect_equal(back, values, tolerance = 1e-6)
  expect_lt(max(abs(back - values)), 1e-2)
})

test_that("integer labels are written exactly (#289)", {
  labels <- as.numeric(rep_len(0:1000, n_vox))
  expect_identical(write_and_read(labels), labels)
})

test_that("a constant volume can be written", {
  ones <- rep(1, n_vox)
  expect_identical(write_and_read(ones), ones)
})
