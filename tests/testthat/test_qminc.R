library(testthat)
library(batchtools)
context("Parallel Functions")

## Write all tests as a function
## Skipping data download if on cran or travis
## downloading a dir local copy of RMINC test data
## so that queue jobs can be tested.
## not my favourite code ever.
test_parallel <-
  function(){

    old_options <- options()
    options(batchtools.verbose = getOption("verbose"))
    on.exit({ options(old_options) })
    
    if(identical(Sys.getenv("NOT_CRAN"), "true")
       && !identical(Sys.getenv("TRAVIS"), "true")){

      if(!exists("dataPath")){
        dataPath <- "."
      }
      
      if(!file.exists(file.path(dataPath, "RMINC-test-data-main/rminctestdata"))){
        
        if(getOption("verbose")){
          getRMINCTestData(dataPath)
        } else {
          capture.output(getRMINCTestData(dataPath), type = "message")
        }
        
        on.exit({
          unlink(file.path(dataPath, "RMINC-test-data-main/rminctestdata"), recursive = TRUE)
          unlink(file.path("rminctestdata.tar.gz"))
        }, add = TRUE)
        
      }
      
      dataPath <- file.path(dataPath, "RMINC-test-data-main/rminctestdata/")
      
      gf <- read.csv(file.path(dataPath, "test_data_set.csv"), stringsAsFactors = TRUE)
      mask_file <- file.path(dataPath, "testminc-mask.mnc")          
    }

    test_that("Test sequential, multicore, and queue applies work", {
      
      skip_on_cran()
      skip_on_travis()
      
      if(Sys.getenv("TEST_Q_MINC") != "yes") 
        skip("qMinc tests disabled")
        
        m_sequential <- verboseRun(
          mincApplyRCPP(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10)),
          getOption("verbose"))
        
        m_queue <- verboseRun(
          qMincApply(gf$jacobians_fixed_2, mean, slab_sizes = c(5,1,10)),
          getOption("verbose")
        )
        
        m_multicore <- verboseRun(
          mcMincApply(gf$jacobians_fixed_2, mean, slab_sizes = c(10,10,10)
                    , cores = min(2, parallel::detectCores() - 1)),
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
  }

################################
## Run the tests
test_parallel()
