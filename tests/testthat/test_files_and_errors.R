context("Temporary files and read errors")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"), stringsAsFactors = TRUE)
gf$jacobians_fixed_2 <- as.character(gf$jacobians_fixed_2)

test_that("local parallel runs do not write to the working directory (#211)", {
  work <- tempfile("rminc_cwd_")
  dir.create(work)
  old_wd <- setwd(work)
  on.exit(
    {
      setwd(old_wd)
      Sys.chmod(work, "0755")
      unlink(work, recursive = TRUE)
    },
    add = TRUE
  )

  # Make the working directory read-only, as an installed package
  # directory often is. A root user can still write to it.
  Sys.chmod(work, "0555")
  skip_if(
    file.access(work, 2) == 0,
    "cannot make a read-only directory (root user?)"
  )

  # The example from the pMincApply documentation (#299)
  testFunc <- function(x) c(1, 2)
  verboseRun(
    pout <- pMincApply(gf$jacobians_fixed_2, testFunc, local = TRUE, cores = 2)
  )
  expect_equal(ncol(pout), 2)
  expect_true(all(pout[, 1] == 1))

  verboseRun(
    plm <- mincLm(jacobians_fixed_2 ~ Sex, gf, parallel = c("local", 2))
  )
  verboseRun(slm <- mincLm(jacobians_fixed_2 ~ Sex, gf))
  expect_equal(unclass(plm), unclass(slm), check.attributes = FALSE)

  expect_equal(list.files(work, all.files = TRUE, no.. = TRUE), character(0))
})

# Copy a MINC file and overwrite the end of the data with a fixed pattern.
# The header stays good, so the file opens, but the data does not read.
corrupt_copy <- function(src) {
  bad <- tempfile("rminc_corrupt_", fileext = ".mnc")
  file.copy(src, bad)
  n <- file.size(bad)
  bytes <- readBin(bad, "raw", n)
  idx <- seq(floor(n * 0.7), floor(n * 0.99))
  bytes[idx] <- as.raw(seq_along(idx) %% 256)
  writeBin(bytes, bad)
  bad
}

quietly <- function(expr) {
  capture.output(res <- expr)
  res
}

test_that("read errors name the file that has the problem (#53)", {
  bad <- corrupt_copy(gf$jacobians_fixed_2[1])
  on.exit(unlink(bad), add = TRUE)

  opens <- tryCatch(
    {
      minc.dimensions.sizes(bad)
      TRUE
    },
    error = function(e) FALSE
  )
  skip_if_not(opens, "the corrupt file does not open with this HDF5")
  reads <- tryCatch(
    {
      mincGetVolume(bad)
      TRUE
    },
    error = function(e) FALSE
  )
  skip_if(reads, "the corrupt file reads without error with this HDF5")

  expect_error(mincGetVolume(bad), basename(bad), fixed = TRUE)

  bad_gf <- gf
  bad_gf$jacobians_fixed_2[nrow(bad_gf)] <- bad
  expect_error(
    quietly(mincLm(jacobians_fixed_2 ~ Sex, bad_gf)),
    basename(bad),
    fixed = TRUE
  )
  expect_error(
    quietly(mincMean(bad_gf$jacobians_fixed_2)),
    basename(bad),
    fixed = TRUE
  )
})
