library(testthat)
context("Parallel Functions")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")
 
gf <- read.csv(file.path(dataPath, "test_data_set.csv"))
mask_file <- file.path(dataPath, "testminc-mask.mnc")

test_that("Test sequential, multicore, and queue applies work", {
  
  skip_on_cran()
  skip_on_travis()
  
  if(Sys.getenv("TEST_Q_MINC") != "yes") 
    skip("qMinc tests disabled")
  
  options(BBmisc.ProgressBar.style = "off")
    
  m_sequential <- verboseRun(
    mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10)),
    getOption("verbose"))
  
  m_queue <- verboseRun(
    qMincApply(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10)),
    getOption("verbose")
  )
  
  m_multicore <- verboseRun(
    mcMincApply(gf$jacobians_fixed_2, mean, slab_sizes = c(10,10,10), cores = min(2, parallel::detectCores() - 1)),
    getOption("verbose")
  )
  
  expect_equivalent(m_sequential, m_queue)
  expect_equivalent(m_sequential, m_multicore)
  expect_equal(attr(m_sequential, "filenames"), attr(m_queue, "filenames"))
  expect_equal(attr(m_sequential, "filenames"), attr(m_multicore, "filenames"))
  expect_equal(attr(m_sequential, "likeVolume"), attr(m_queue, "likeVolume"))
  expect_equal(attr(m_sequential, "likeVolume"), attr(m_multicore, "likeVolume"))
})

test_that("Masking qMincApply behaves as expected", {
  
  skip_on_cran()
  skip_on_travis()
  
  if(Sys.getenv("TEST_Q_MINC") != "yes") 
    skip("qMinc tests disabled")
  
  m_sequential <- verboseRun(
    mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10), mask = mask_file),
    getOption("verbose"))
  
  m_queue <- verboseRun(
    qMincApply(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10), mask = mask_file),
    getOption("verbose")
  )
  
  expect_equivalent(m_sequential, m_queue)
})

setConfig(old_config)
options(old_options)

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

