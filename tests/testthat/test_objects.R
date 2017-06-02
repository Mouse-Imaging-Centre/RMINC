verboseRun({
  library(testthat)
  
  context("object tests")
  
  getRMINCTestData()
  dataPath <- file.path(tempdir(), "rminctestdata/")
  
  gf <- read.csv(file.path(dataPath, "test_data_set.csv"))
  
  
  a <- structure(1, class = "mincSingleDim")
  b <- structure(1, class = "mincMultiDim")
  c <- structure(1, class = "mincList")
  d <- 1
  e <- as.matrix(d)
  
  test_that("Minc class checks work",{
    expect_true(is.minc(a))
    expect_true(is.minc(b))
    expect_true(is.minc(c))
    expect_true(!is.minc(d))
    expect_true(!is.minc(e))
  })
  
  volume <- mincGetVolume(gf$jacobians_fixed_2[1])
  volume_numeric <- as.numeric(volume)
  
  volume_reconstituted <- as.minc(volume_numeric)
  
  test_that("conversion to and from mincSingleDim",{
    expect_equal(class(volume_reconstituted), class(volume))
    expect_equivalent(volume_reconstituted, volume)
  })
  
  volume_rec_with_attrs <-
    setMincAttributes(volume_reconstituted, mincAttributes(volume))
  
  test_that("minc attributes can be copied",
            expect_equal(volume_rec_with_attrs, volume))
  
  multidim <- mincLm(jacobians_fixed_2 ~ Sex, data = gf)
  
  multidim_matrix <- as.numeric(multidim) #as.matrix does nothing, because it is one
  dim(multidim_matrix) <- dim(multidim)
  
  multidim_reconstituted <- as.minc(multidim_matrix)
  class(multidim_reconstituted) <- c("mincLm", class(multidim_reconstituted))
  multidim_reconstituted <- setMincAttributes(multidim_reconstituted, mincAttributes(multidim))
  
  
  test_that("conversion to and from mincMultiDim", {
    # Ensure values work
    expect_equal(multidim, multidim_reconstituted, check.attributes = FALSE)
    
    # Ensure attributes are equal up to reordering
    expect_equal(attributes(multidim), 
                 attributes(multidim_reconstituted)[names(attributes(multidim))])
  })
  
  
  single_list <- list(1, 2, getOption("RMINC_MASKED_VALUE"))
  single_list_minc <- simplify2minc(single_list)
  
  test_that("singleDim collapse with masked works", {
    expect_true(inherits(single_list_minc, "mincSingleDim"))
    expect_equal(length(single_list_minc), 3)
  })
  
  multi_list <- list(getOption("RMINC_MASKED_VALUE"), c(1:4), c(1:4))
  multi_list_minc <- simplify2minc(multi_list)
  
  test_that("multiDim collapse with masked works",{
    expect_true(inherits(multi_list_minc, "mincMultiDim"))
    expect_equal(dim(multi_list_minc), c(3,4))
  })
  
  list_list <- list(list(1), list(2), getOption("RMINC_MASKED_VALUE"))
  list_list_minc <- simplify2minc(list_list)
  
  test_that("list collapse with masked works", {
    expect_true(inherits(list_list_minc, "mincList"))
    expect_equal(length(list_list_minc), 3)
  })
})


