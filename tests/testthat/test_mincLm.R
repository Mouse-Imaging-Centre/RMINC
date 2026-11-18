library(testthat)
library(dplyr)
context("mincLm - two group test")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"), stringsAsFactors = TRUE)

voxel_left <- mincGetVoxel(gf$jacobians_fixed_2[1:10], 0, 0, 0)
voxel_right <- mincGetVoxel(gf$jacobians_fixed_2[11:20], 0, 0, 0)
Sex <- gf$Sex[1:10]
Scale <- gf$scale[1:10]
Coil <- as.factor(gf$coil[1:10])

gf$coil <- as.factor(gf$coil)
gftest <- gf[1:10, ]
gftest$voxel_right <- (gf$jacobians_fixed_2[11:20])
gftest$voxel_left_file <- gf$jacobians_fixed_2[1:10]

rmod <- lm(voxel_left ~ Sex)
rLm = summary(rmod)
# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm = verboseRun(
  "mincLm(voxel_left_file ~ Sex, gftest)",
  getOption("verbose")
)


test_that("mincLm Two Factors", {
  expect_equal(
    rmincLm[1, 1:6],
    with(
      rLm,
      c(
        fstatistic[1],
        r.squared[1],
        coefficients[, "Estimate"],
        coefficients[, "t value"]
      )
    ),
    check.attributes = FALSE
  )

  expect_that(attr(rmincLm, "df")[[2]], is_equivalent_to(rLm$df[2]))
})

test_that("Likelihood and information criteria are computed correctly", {
  expect_equal(as.numeric(rmincLm[1, "logLik"]), as.numeric(logLik(rmod)))
  expect_equal(as.numeric(AIC(rmincLm)[1]), as.numeric(AIC(rmod)))
  expect_equal(as.numeric(BIC(rmincLm)[1]), as.numeric(BIC(rmod)))
})

# now test that findPeaks works
context("mincFindPeaks - ensure we can find correct peaks")

has_find_peaks <- as.character(Sys.which("find_peaks")) != ""

if (has_find_peaks) {
  verboseRun({
    # get rid of NAs
    rmincLm[is.na(rmincLm)] <- 0
    # find peaks
    peaks <- mincFindPeaks(rmincLm, "tvalue-SexM", minDistance = 1)
    # find the min and max by hand
    rmincLmArray <- mincArray(rmincLm, "tvalue-SexM")
    maxPeak <- arrayInd(which.max(rmincLmArray), .dim = c(10, 10, 10))
    minPeak <- arrayInd(which.min(rmincLmArray), .dim = c(10, 10, 10))
    minPeakFromPeaks <- as.integer(arrange(peaks, value)[1, 1:3])
    maxPeakFromPeaks <- as.integer(peaks[1, 1:3])
  })
}

test_that("mincFindPeaks min and max", {
  skip_if_not(has_find_peaks, "find_peaks binary not installed")
  expect_that(maxPeak[1], is_equivalent_to(maxPeakFromPeaks[1]))
  expect_that(maxPeak[2], is_equivalent_to(maxPeakFromPeaks[2]))
  expect_that(maxPeak[3], is_equivalent_to(maxPeakFromPeaks[3]))
  expect_that(minPeak[1], is_equivalent_to(minPeakFromPeaks[1]))
  expect_that(minPeak[2], is_equivalent_to(minPeakFromPeaks[2]))
  expect_that(minPeak[3], is_equivalent_to(minPeakFromPeaks[3]))
})

context("vertexFindPeaks - ensure matches mincFindPeaks")
test_that("vertexFindPeaks matches mincFindPeaks", {
  skip_if_not(has_find_peaks, "find_peaks binary not installed")
  adj <- RMINC:::neighbour_list(10, 10, 10, 6)
  g <- igraph::graph_from_adj_list(lapply(adj, function(nebs) nebs + 1))
  pos_peaks <-
    vertexFindPeaks(
      rmincLm[, "tvalue-SexM"],
      g,
      direction = "positive",
      mindist = 1,
      output = "indices"
    )

  pos_peak <- arrayInd(
    pos_peaks[which.max(rmincLm[pos_peaks, "tvalue-SexM"])],
    .dim = rep(10, 3)
  )

  neg_peaks <-
    vertexFindPeaks(
      rmincLm[, "tvalue-SexM"],
      g,
      direction = "negative",
      mindist = 1,
      output = "indices"
    )

  neg_peak <- arrayInd(
    neg_peaks[which.min(rmincLm[neg_peaks, "tvalue-SexM"])],
    .dim = rep(10, 3)
  )

  expect_equal(pos_peak[1:3], maxPeak[1:3])
  expect_equal(neg_peak[1:3], minPeak[1:3])
})

context("mincLm - two group test with interaction")

# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm2 = verboseRun(
  "mincLm(voxel_left_file~Sex*Scale,gftest)",
  getOption("verbose")
)
gftest$voxel_left = voxel_left
rLm2 = summary(lm(voxel_left ~ Sex * Scale, gftest))

test_that("mincLm interaction", {
  expect_equal(
    rmincLm2[1, 1:10],
    with(
      rLm2,
      c(
        fstatistic[1],
        r.squared[1],
        coefficients[, "Estimate"],
        coefficients[, "t value"]
      )
    ),
    check.attributes = FALSE
  )

  expect_that(attr(rmincLm2, "df")[[2]], is_equivalent_to(rLm2$df[2]))
})


context("mincLm - three group test")

# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm3 <- verboseRun(
  "mincLm(voxel_left_file~coil,gftest)",
  getOption("verbose")
)

gftest$voxel_left = voxel_left
rLm3 = summary(lm(voxel_left ~ Coil, gftest))

test_that("mincLm Three Factors", {
  expect_equal(
    rmincLm3[1, 1:8],
    with(
      rLm3,
      c(
        fstatistic[1],
        r.squared[1],
        coefficients[, "Estimate"],
        coefficients[, "t value"]
      )
    ),
    check.attributes = FALSE
  )
  expect_that(attr(rmincLm3, "df")[[2]], is_equivalent_to(rLm3$df[2]))
})

context("mincLm - three group test with interaction")

# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm4 = verboseRun(
  "mincLm(voxel_left_file~Scale*Coil,gftest)",
  getOption("verbose")
)

gftest$voxel_left = voxel_left
rLm4 = summary(lm(voxel_left ~ Scale * Coil, gftest))

test_that("mincLm Three Factors Interaction", {
  expect_equal(
    rmincLm4[1, 1:14],
    with(
      rLm4,
      c(
        fstatistic[1],
        r.squared[1],
        coefficients[, "Estimate"],
        coefficients[, "t value"]
      )
    ),
    check.attributes = FALSE
  )
  expect_that(attr(rmincLm4, "df")[[2]], is_equivalent_to(rLm4$df[2]))
})

context("mincLm - Model Selection")

test_that("Model Selection Works", {
  comp1 <- compare_models(rmincLm, rmincLm2, rmincLm3, rmincLm4, metric = AIC)
  comp2 <- compare_models(rmincLm, rmincLm2, rmincLm3, rmincLm4)

  params <- sapply(list(rmincLm, rmincLm2, rmincLm3, rmincLm4), function(mod) {
    ncol(attr(mod, "model")) + 1
  })

  n <- nrow(gftest)

  aicc_corr <- 2 * (params + 1) * (params + 2) / (n - params - 2)
  aicc_corr_mat <- t(matrix(rep(aicc_corr, nrow(comp1)), ncol = nrow(comp1)))

  expect_equal(
    unclass(comp2),
    unclass(comp1 + aicc_corr_mat),
    check.attributes = FALSE
  )
  expect_equal(as.numeric(summary(comp2)$wins), c(100, 100, 119, 981))
})

context("mincLm - Parallel")
test_that("mincLm local multicore works", {
  skip_on_cran()
  skip_on_travis()

  if (Sys.getenv("TEST_Q_MINC") != "yes") {
    skip("qMinc tests disabled")
  }

  verboseRun(
    prlm <- mincLm(
      voxel_left_file ~ Scale * Coil,
      gftest,
      parallel = c("local", 3)
    )
  )

  expect_equal(rmincLm4, prlm, check.attributes = FALSE)
})

# test_that("mincLm queue parallel works", {
#   verboseRun(qrlm <- mincLm(voxel_left_file~Scale*Coil,gftest, parallel = c("sge", 2)))
#
#   expect_equal(rmincLm, qrlm, check.attributes = FALSE)
# })

## Need a much better test here.
context("mincLm - Randomize")
test_that("mincLm randomize works", {
  verboseRun(rlm <- mincLm(voxel_left_file ~ Sex, gftest))
  verboseRun(rand <- RMINC:::mincRandomize.mincLm(rlm, R = 10))

  thresh <- thresholds(rand)
  expect_equal(dim(thresholds(rand)), c(4, 2))
})

test_that("mincRandomize works with one column, by name or index (#318)", {
  verboseRun(rlm <- mincLm(voxel_left_file ~ Sex, gftest))

  verboseRun(
    rand_f <- mincRandomize(rlm, R = 3, columns = "F-statistic")
  )
  expect_equal(dim(rand_f$randomization_dist), c(3, 1))
  expect_equal(colnames(rand_f$randomization_dist), "F-statistic")
  expect_equal(colnames(thresholds(rand_f)), "F-statistic")
  expect_equal(ncol(rand_f$stats), 1)

  t_col <- grep("tvalue-", colnames(rlm))[2]
  verboseRun(rand_t <- mincRandomize(rlm, R = 3, columns = t_col))
  expect_equal(colnames(rand_t$randomization_dist), colnames(rlm)[t_col])

  expect_error(
    mincRandomize(rlm, R = 3, columns = "no-such-column"),
    "Unknown column"
  )
})

test_that("mincLm accepts a mask as voxel values (#258)", {
  like <- as.character(gftest$voxel_left_file[1])
  n_vox <- prod(minc.dimensions.sizes(like))
  mask_vals <- rep(c(0, 1), length.out = n_vox)
  mask_file <- tempfile(fileext = ".mnc")
  on.exit(unlink(mask_file))
  verboseRun(
    mincWriteVolume(mask_vals, mask_file, like.filename = like, clobber = TRUE)
  )

  verboseRun(no_mask <- mincLm(voxel_left_file ~ Sex, gftest))
  verboseRun(from_file <- mincLm(voxel_left_file ~ Sex, gftest, mask = mask_file))
  verboseRun(from_vec <- mincLm(voxel_left_file ~ Sex, gftest, mask = mask_vals))
  verboseRun(
    from_vol <- mincLm(
      voxel_left_file ~ Sex,
      gftest,
      mask = mincGetVolume(mask_file)
    )
  )
  verboseRun(
    from_par <- mincLm(
      voxel_left_file ~ Sex,
      gftest,
      mask = mask_vals,
      parallel = c("local", 2)
    )
  )

  expect_equal(unclass(from_vec), unclass(from_file), check.attributes = FALSE)
  expect_equal(unclass(from_vol), unclass(from_file), check.attributes = FALSE)
  expect_equal(
    unclass(from_par)[mask_vals > .5, ],
    unclass(from_file)[mask_vals > .5, ],
    check.attributes = FALSE
  )
  expect_equal(
    unclass(from_vec)[mask_vals > .5, ],
    unclass(no_mask)[mask_vals > .5, ],
    check.attributes = FALSE
  )

  # every column is zero outside the mask, not uninitialised memory
  expect_true(all(unclass(from_file)[mask_vals < .5, ] == 0))
  expect_true(all(unclass(from_vec)[mask_vals < .5, ] == 0))

  expect_error(
    mincLm(voxel_left_file ~ Sex, gftest, mask = mask_vals[-1]),
    "voxels"
  )
})

test_that("an unknown parallel backend is an error (#215)", {
  expect_error(
    mincLm(voxel_left_file ~ Sex, gftest, parallel = c("parallel", 2)),
    "Unknown parallel backend \"parallel\""
  )

  expect_true(RMINC:::parallel_is_local(c("local", 2)))
  expect_true(RMINC:::parallel_is_local(c("snowfall", 2)))
  expect_false(RMINC:::parallel_is_local(c("sge", 2)))
  expect_false(RMINC:::parallel_is_local(c("batchtools", 2)))
  expect_error(RMINC:::parallel_is_local(c("local", 0)), "positive whole")
  expect_error(RMINC:::parallel_is_local("local"), "two element")
})
