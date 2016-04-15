library(testthat)
context("test_qminc_apply")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")

ms <- verboseRun(
  "mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10))", 
  getOption("verbose"))

mq <- verboseRun(
  "qMincApply(gf$jacobians_fixed_2, mean, parallel_method = 'multicore', slab_sizes = c(5,1,10))",
  getOption("verbose")
)

mp <- verboseRun(
  "mcMincApply(gf$jacobians_fixed_2, mean, slab_sizes = c(10,10,10))",
  getOption("verbose")
)

test_that("Sequential mincApplyRCPP matches qMincApply",
          expect_equivalent(ms, mq))

test_that("Sequential mincApplyRCPP matches mcMincApply",
          expect_equivalent(ms, mp))

test_that("Attributes are set correctly", {
  expect_equal(attr(ms, "filenames"), attr(mq, "filenames"))
  expect_equal(attr(ms, "filenames"), attr(mp, "filenames"))
  expect_equal(attr(ms, "likeVolume"), attr(mq, "likeVolume"))
  expect_equal(attr(ms, "likeVolume"), attr(mp, "likeVolume"))
})

# if(getOption("RMINC_QUEUE") == "sge"){
#   gf2 <- gf
#   gf2$jacobians_fixed_2 <- sub("/tmp/", "/hpf/largeprojects/MICe/chammill/", gf$jacobians_fixed_2)
# 
#   ms_2 <-  verboseRun(
#     "mincApplyRCPP(gf2$jacobians_fixed_2, mean, slab_sizes = c(5,1,10))",
#     getOption("verbose"))
# 
#   mp_sge <- verboseRun("qMincApply(gf2$jacobians_fixed_2, mean, parallel_method = 'sge',
#                        slab_sizes = c(5,1,10), temp_dir = normalizePath(\"~/sge_reg\"))",
#                        getOption("verbose"))
# 
#   test_that("SGE parallel matches sequential", {
#     expect_equivalent(mp_sge, ms_2)
#     expect_equal(attr(mp_sge, "filenames"), attr(ms_2, "filenames"))
#     expect_equal(attr(mp_sge, "likeVolume"), attr(ms_2, "likeVolume"))
#   })
# }

