library(testthat)
context("anatApply")

if(!exists("dataPath")) {
  dataPath <- tempdir()
}
  
getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "CIVET_TEST.csv"))
gf <- civet.getAllFilenames(gf,"ID","TEST",file.path(dataPath, "CIVET"),"TRUE","1.1.12")
gf <- civet.readAllCivetFiles(file.path(dataPath, "AAL.csv"),gf)

mm = subset(gf$lobeThickness,gf$Primary.Diagnosis=="ADHD")
mm = mean(mm[,1])

ma <- verboseRun("anatApply(gf$lobeThickness,gf$Primary.Diagnosis)",getOption("verbose"))

test_that("anatApply one output", {
       expect_equal(ma[1], mm) 
})
