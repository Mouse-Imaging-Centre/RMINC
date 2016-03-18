requireNamespace("testthat")
context("mincApply")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")


mm <- verboseRun("mincMean(gf$jacobians_fixed_2)",getOption("verbose"))
ma <- verboseRun("mincApply(gf$jacobians_fixed_2,quote(mean(x)))",getOption("verbose"))

test_that("mincApply one output",{
    for (nVox in 1:length(mm)) {
 	   expect_equal(ma[nVox], mm[nVox]) }
})


# Need to define global variable when running tests, but normally do not...
testFunc <<- function (x) { return(c(1,2))}
ma <- verboseRun("mincApply(gf$jacobians_fixed_2,quote(testFunc(x)))",getOption("verbose"))

test_that("mincApply two output",{
    for (nVox in 1:dim(ma)[1]) {
 	   expect_equal(ma[nVox,1], 1) 
	   expect_equal(ma[nVox,2], 2)}
})


# Add a test to check if snowfall is installed.
# Do this through a tryCatch block
result = tryCatch({
	library(snowfall)
	result = TRUE
}, error = function(e) {
	result = FALSE
})


test_that("snowfall is installed",{
	expect_true(result)
})

# Add a test to check if RSGE is installed.
# Do this through a tryCatch block
result = tryCatch({
	library(Rsge)
	result = TRUE
}, error = function(e) {
	result = FALSE
})


test_that("Rsge is installed",{
	expect_true(result)
})





if("package:snowfall" %in% search()) {

  ma <- verboseRun("pMincApply(gf$jacobians_fixed_2,quote(mean(x)),global = 'gf',packages = 'car')",getOption("verbose"))

  test_that("pmincapply snowfall",{
    for (nVox in 1:length(mm)) {
      expect_that(ma[nVox],equals(mm[nVox],tolerance = 0.00001))
    }
  })

  ma <- verboseRun("pMincApply(gf$jacobians_fixed_2,quote(testFunc(x)),global = c('testFunc','gf'),packages = c('car','stats'))",getOption("verbose"))

  test_that("pmincApply snowfall two output",{
    for (nVox in 1:dim(ma)[1]) {
      expect_equal(ma[nVox,1], 1) 
      expect_equal(ma[nVox,2], 2)}
  })
}
