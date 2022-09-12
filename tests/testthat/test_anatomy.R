library(testthat)
context("anatApply")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "RMINC-test-data-main/rminctestdata/")

gf <- read.csv(file.path(dataPath, "CIVET_TEST.csv"), stringsAsFactors = TRUE)
gf <- civet.getAllFilenames(gf,"ID","TEST",file.path(dataPath, "CIVET"),"TRUE","1.1.12")
gf <- civet.readAllCivetFiles(file.path(dataPath, "AAL.csv"),gf)

mm <- subset(gf$lobeThickness,gf$Primary.Diagnosis=="ADHD")
mm <- mean(mm[,1])

test_that("anatApply one output", {
  ma <- verboseRun("anatApply(gf$lobeThickness,gf$Primary.Diagnosis)",getOption("verbose"))
  expect_equal(ma[1], mm) 
})

context("anatMean")

#Calculate mean
test_that("anatMean", {
  vm <- verboseRun("anatMean(gf$lobeThickness)",getOption("verbose"))
  expect_equal(colMeans(gf$lobeThickness), vm)
})

context("anatSum")

#Calculate sum
test_that("anatSum", {
  vs <- verboseRun("anatSum(gf$lobeThickness)",getOption("verbose"))
  expect_equal(colSums(gf$lobeThickness), vs)
})

context("anatVar")

#Calculate variance
test_that("anatVar", {
  vv <- verboseRun("anatVar(gf$lobeThickness)",getOption("verbose"))
  expect_equal(apply(gf$lobeThickness, 2, var), vv)
})


context("anatSd")

#Calculate standard deviation
test_that("anatSd", {
  vsd <- verboseRun("anatSd(gf$lobeThickness)",getOption("verbose"))
  expect_equal(apply(gf$lobeThickness, 2, sd), vsd)
})
