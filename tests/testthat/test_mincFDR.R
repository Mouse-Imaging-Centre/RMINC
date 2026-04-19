library(testthat)

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

rLm = summary(lm(voxel_left ~ Sex))


rmincLm <- verboseRun(
  "mincLm(voxel_left_file ~ Sex, gftest)",
  getOption("verbose")
)


rLmFDR1 = p.adjust(pt2(rmincLm[, 5], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 6], attr(rmincLm, "df")[[3]]), "fdr")

rmincFDR = verboseRun("mincFDR(rmincLm)", getOption("verbose"))

test_that("mincFDR Two Factors", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "mincLm(voxel_left_file~Sex*Scale,gftest)",
  getOption("verbose")
)

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left ~ Sex * Scale, gftest))

rLmFDR1 = p.adjust(pt2(rmincLm[, 7], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 8], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 9], attr(rmincLm, "df")[[4]]), "fdr")
rLmFDR4 = p.adjust(pt2(rmincLm[, 10], attr(rmincLm, "df")[[5]]), "fdr")

rmincFDR = verboseRun("mincFDR(rmincLm)", getOption("verbose"))


test_that("mincFDR interaction", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[1]), unname(rmincFDR[1, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[2]), unname(rmincFDR[2, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[3]), unname(rmincFDR[3, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[1]), unname(rmincFDR[1, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[2]), unname(rmincFDR[2, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[3]), unname(rmincFDR[3, 5]), ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "rmincLm = mincLm(voxel_left_file~coil,gftest)",
  getOption("verbose")
)

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left ~ Coil, gftest))


rLmFDR1 = p.adjust(pt2(rmincLm[, 6], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 7], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 8], attr(rmincLm, "df")[[4]]), "fdr")

rmincFDR = verboseRun("mincFDR(rmincLm)", getOption("verbose"))

test_that("mincFDR Three Factors", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[1]), unname(rmincFDR[1, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[2]), unname(rmincFDR[2, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[3]), unname(rmincFDR[3, 4]), ignore_attr = TRUE)
})


rmincLm <- verboseRun(
  "mincLm(voxel_left_file~Scale*coil,gftest)",
  getOption("verbose")
)

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left ~ Scale * Coil, gftest))


rLmFDR1 = p.adjust(pt2(rmincLm[, 9], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 10], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 11], attr(rmincLm, "df")[[4]]), "fdr")
rLmFDR4 = p.adjust(pt2(rmincLm[, 12], attr(rmincLm, "df")[[5]]), "fdr")
rLmFDR5 = p.adjust(pt2(rmincLm[, 13], attr(rmincLm, "df")[[6]]), "fdr")
rLmFDR6 = p.adjust(pt2(rmincLm[, 14], attr(rmincLm, "df")[[7]]), "fdr")


rmincFDR = verboseRun("mincFDR(rmincLm)", getOption("verbose"))

test_that("mincLm Three Factors Interaction", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[1]), unname(rmincFDR[1, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[2]), unname(rmincFDR[2, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[3]), unname(rmincFDR[3, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[1]), unname(rmincFDR[1, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[2]), unname(rmincFDR[2, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[3]), unname(rmincFDR[3, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR5[1]), unname(rmincFDR[1, 6]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR5[2]), unname(rmincFDR[2, 6]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR5[3]), unname(rmincFDR[3, 6]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR6[1]), unname(rmincFDR[1, 7]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR6[2]), unname(rmincFDR[2, 7]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR6[3]), unname(rmincFDR[3, 7]), ignore_attr = TRUE)
})
