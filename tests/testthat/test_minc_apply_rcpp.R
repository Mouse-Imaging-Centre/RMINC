library(testthat)
context("mincApplyRCPP")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"))

mm <- verboseRun("mincMean(gf$jacobians_fixed_2)",getOption("verbose"))
ma <- verboseRun(
    "mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(10,10,1))", 
    getOption("verbose"))

#Coerce numeric vector produced by mincApplyRCPP to mincMultiDim
dim(ma) <- c(length(ma), 1)
class(ma) <- c("mincMultiDim", "matrix")


test_that("mincApplyRCPP matches MincMean",{
  expect_equivalent(ma, mm)
})

setRMINCMaskedValue(val = 0)
mm_masked <- verboseRun("mincMean(gf$jacobians_fixed_2, 
                        mask = file.path(dataPath, 'testminc-mask.mnc'))", 
                        getOption("verbose"))
ma_masked <- verboseRun("mincApplyRCPP(gf$jacobians_fixed_2, mean, 
                        mask = file.path(dataPath, 'testminc-mask.mnc'),
                        slab_sizes = c(1,5,5))", 
                        getOption("verbose"))

dim(ma_masked) <- c(length(ma_masked), 1)
class(ma_masked) <- c("mincMultiDim", "matrix")


test_that("mincApplyRCPP masks like mincMean", {
  expect_equivalent(mm_masked, ma_masked)
  expect_equal(attr(mm_masked, "filenames"), attr(ma_masked, "filenames"))
})


# # #Generate Demo Mask
# demo_mask <- rep(1, 15^3)
# demo_mask[sample(1:(15^3), 100)] <- 0
# mincWriteVolume(demo_mask, output.filename = "/tmp/rminctestdata/mask15.mnc",
#                 like.filename = "/tmp/rminctestdata/absolute_jacobian_file_1.mnc")
# 
# test_files <-
#   list.files("/tmp/rminctestdata/", pattern = "testminc[0-9]+.*\\.mnc", full.names = TRUE)
# 
# test_mask <-
#   "/tmp/rminctestdata/testminc-mask.mnc"
# 
# 
# system.time({
#   rcpp_mean <-
#     mincApplyRCPP(test_files,
#                   mean,
#                   mask = test_mask,
#                   filter_masked = TRUE,
#                   slab_sizes = c(1,10,10),
#                   collate = unlist)
# })
# 
# system.time({
#   r_mean <-
#     mincApply(test_files, quote(mean(x)), mask = test_mask, reduce = TRUE)
# })
# 
# system.time({
#   orig_mean <-
#     mincMean(test_files,
#              mask = test_mask)
# })



