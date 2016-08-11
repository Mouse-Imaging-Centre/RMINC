library(testthat)
context("mincApply")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"))

mm <- verboseRun("mincMean(gf$jacobians_fixed_2)",getOption("verbose"))
ma <- verboseRun("mincApply(gf$jacobians_fixed_2,quote(mean(x)))",getOption("verbose"))

test_that("mincApply one output",{
    expect_equivalent(as.numeric(mm), as.numeric(ma))
})


# Need to define global variable when running tests, but normally do not...
testFunc <<- function (x) { return(c(1,2))}
ma <- verboseRun("mincApply(gf$jacobians_fixed_2, quote(testFunc(x)))",getOption("verbose"))

test_that("mincApply two output",{
 	   expect_equal(ma[,1], rep(1, length(ma[,1]))) 
	   expect_equal(ma[,2], rep(2, length(ma[,2])))
})

# Core setting prevents memory explosion on many core machines 
ma <- verboseRun("pMincApply(gf$jacobians_fixed_2, mean, cores = min(2, parallel::detectCores() - 1))",getOption("verbose"))
dim(ma) <- c(length(ma), 1)
class(ma) <- "mincMultiDim"

test_that("pmincapply local",{
  expect_equivalent(ma, mm)
})

ma <- verboseRun("pMincApply(gf$jacobians_fixed_2, testFunc, cores = min(2, parallel::detectCores() - 1))",
                 getOption("verbose"))

test_that("pmincApply snowfall two output",{
  expect_equal(ma[,1], rep(1, length(ma[,1]))) 
  expect_equal(ma[,2], rep(2, length(ma[,2])))
})

