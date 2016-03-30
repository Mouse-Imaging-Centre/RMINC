library("testthat")
context("mincApplyRCPP")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")

mm <- verboseRun("mincMean(gf$jacobians_fixed_2)",getOption("verbose"))
ma <- verboseRun(
    "mincApplyRCPP(gf$jacobians_fixed_2, mean, collate = unlist, slab_sizes = c(10,10,1))", 
    getOption("verbose"))

#Coerce numeric vector produced by mincApplyRCPP to mincMultiDim
dim(ma) <- c(length(ma), 1)
class(ma) <- c("mincMultiDim", "matrix")


test_that("mincApplyRCPP matches MincMean",{
  expect_equal(ma, mm)
})

setRMINCMaskedValue(val = 0)
mm_masked <- verboseRun("mincMean(gf$jacobians_fixed_2, 
                        mask = '/tmp/rminctestdata/testminc-mask.mnc')", 
                        getOption("verbose"))
ma_masked <- verboseRun("mincApplyRCPP(gf$jacobians_fixed_2, mean, 
                        mask = '/tmp/rminctestdata/testminc-mask.mnc',
                        slab_sizes = c(1,5,5))", 
                        getOption("verbose"))

#Coerce numeric vector produced by mincApplyRCPP to mincMultiDim
ma_masked <- unlist(ma_masked)
mm_masked <- as.numeric(mm_masked)

test_that("mincApplyRCPP masks like mincMean", {
  expect_equal(mm_masked, ma_masked)
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



