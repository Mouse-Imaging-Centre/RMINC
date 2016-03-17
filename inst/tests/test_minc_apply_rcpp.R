requireNamespace("testthat")
context("mincApplyRCPP")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")

fuzzy_equal <-
  function(x, y, tol = 10^-10)
    (x - y) < tol

mm <- verboseRun("mincMean(gf$jacobians_fixed_2)",getOption("verbose"))
ma <- verboseRun(
    "mincApplyRCPP(gf$jacobians_fixed_2, mean, collate = unlist)", 
    getOption("verbose"))

#Coerce numeric vector produced by mincApplyRCPP to mincMultiDim
dim(ma) <- c(length(ma), 1)
class(ma) <- c("mincMultiDim", "matrix")


test_that("mincApplyRCPP matches MincMean",{
  expect_equal(ma, mm)
})


# # #Generate Demo Mask
# # demo_mask <- rep(1, 15^3)
# # demo_mask[sample(1:(15^3), 100)] <- 0
# # mincWriteVolume(demo_mask, output.filename = "/tmp/rminctestdata/mask15.mnc", 
# #                 like.filename = "/tmp/rminctestdata/absolute_jacobian_file_1.mnc")
# # 
# # test_files <-
# #   list.files("/tmp/rminctestdata/", pattern = "testminc[0-9]+.*\\.mnc", full.names = TRUE)
# # 
# # test_mask <-
# #   "/tmp/rminctestdata/testminc-mask.mnc"
# #   
# # 
# # system.time({
# #   rcpp_mean <-
# #     mincApplyRCPP(test_files, 
# #                   mean,
# #                   mask = test_mask,
# #                   filter_masked = TRUE,
# #                   collate = unlist) 
# # })
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
# 


