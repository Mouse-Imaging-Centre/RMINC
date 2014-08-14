context("mincApply")

gf <<- read.csv("/tmp/rminctestdata/test_data_set.csv")

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


sink("/dev/null"); library(snowfall); sink();

test_that("snowfall is installed",{
	expect_true("package:snowfall" %in% search())
})


if("package:snowfall" %in% search()) {

  sink("/dev/null")
  ma <- pMincApply(gf$jacobians_fixed_2,quote(mean(x)),global = 'gf')
  sink()
  test_that("pmincapply snowfall",{
    for (nVox in 1:length(mm)) {
      expect_that(ma[nVox],equals(mm[nVox],tolerance = 0.00001))
    }
  })

  sink("/dev/null")
  ma <- pMincApply(gf$jacobians_fixed_2,quote(testFunc(x)),global = c('testFunc','gf'))
  sink()
  test_that("pmincApply snowfall two output",{
    for (nVox in 1:dim(ma)[1]) {
      expect_equal(ma[nVox,1], 1) 
      expect_equal(ma[nVox,2], 2)}
  })
  
  sink()
  sink()
}
