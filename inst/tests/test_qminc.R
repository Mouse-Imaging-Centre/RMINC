library("testthat")
context("test_qminc_apply")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")

ms <- verboseRun(
  "mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10))", 
  getOption("verbose"))

mp <- verboseRun(
  "qMincApply(gf$jacobians_fixed_2, mean, queue = 'multicore', slab_sizes = c(5,1,10))",
  getOption("verbose")
)

test_that("Sequential mincApplyRCPP and Parallel qMincApply agree",
          expect_equal(ms, mp))

if(getOption("RMINC_QUEUE") == "sge"){
  mp_sge <- verboseRun("qMincApply(gf$jacobians_fixed_2, mean, queue = 'sge', slab_sizes = c(5,1,10))",
                       getOption("verbose"))
  
  test_that("SGE parallel matches sequential",
            expect_equal(mp_sge, ms))
}

