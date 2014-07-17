context("mincApply")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")

sink("/dev/null"); mm <- mincMean(gf$jacobians_fixed_2); sink();
sink("/dev/null"); ma <- mincApply(gf$jacobians_fixed_2,quote(mean(x))); sink();

test_that("mincApply one output",{
    for (nVox in 1:length(mm)) {
 	   expect_equal(ma[nVox], mm[nVox]) }
})


testFunc <<- function (x) { return(c(1,2))}
sink("/dev/null"); ma <- mincApply(gf$jacobians_fixed_2,quote(testFunc(x))); sink();

test_that("mincApply two output",{
    for (nVox in 1:dim(ma)[1]) {
 	   expect_equal(ma[nVox,1], 1) 
	   expect_equal(ma[nVox,2], 2)}
})


library(snowfall)
sink("/dev/null"); ma <- pMincApply(gf$jacobians_fixed_2,quote(mean(x))); sink();

test_that("pmincapply snowfall",{
    for (nVox in 1:length(mm)) {
 	   expect_equal(ma[nVox], mm[nVox]) }
})


sink("/dev/null"); ma <- pMincApply(gf$jacobians_fixed_2,quote(testFunc(x)),global = "testFunc"); sink();
test_that("pmincApply snowfall two output",{
    for (nVox in 1:dim(ma)[1]) {
 	   expect_equal(ma[nVox,1], 1) 
	   expect_equal(ma[nVox,2], 2)}
})

gf <- read.csv("/projects/moush/matthijs/2013-08-test-csv-RMINC/test_data_set.csv")

