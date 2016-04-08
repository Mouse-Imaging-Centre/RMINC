library("testthat")
context("test_qminc_apply")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")

ms <- verboseRun(
  "mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10))", 
  getOption("verbose"))

mq <- verboseRun(
  "qMincApply(gf$jacobians_fixed_2, mean, queue = 'multicore', slab_sizes = c(5,1,10))",
  getOption("verbose")
)

mp <- verboseRun(
  "pMincApply2(gf$jacobians_fixed_2, mean, slab_sizes = c(10,10,10)",
  getOption("verbose")
)

test_that("Sequential mincApplyRCPP and Parallel qMincApply agree",
          expect_equal(ms, mq))

test_that("Sequential mincApplyRCPP matches pMincApply2",
          expect_equal(ms, mp))

# if(getOption("RMINC_QUEUE") == "sge"){
#   gf2 <- gf
#   gf2$jacobians_fixed_2 <- sub("/tmp/", "/hpf/largeprojects/MICe/chammill/", gf$jacobians_fixed_2) 
#   
#   ms_2 <-  verboseRun(
#     "mincApplyRCPP(gf2$jacobians_fixed_2, mean, slab_sizes = c(5,1,10))", 
#     getOption("verbose"))
#   
#   mp_sge <- verboseRun("qMincApply(gf2$jacobians_fixed_2, mean, queue = 'sge', 
#                        slab_sizes = c(5,1,10), temp_dir = normalizePath(\"~/sge_reg\"))",
#                        getOption("verbose"))
#   
#   test_that("SGE parallel matches sequential",
#             expect_equal(mp_sge, ms_2))
# }

