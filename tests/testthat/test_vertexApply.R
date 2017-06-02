library(testthat)
context("vertexApply")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gftest <- read.csv(file.path(dataPath, "subject.csv"))

subjectFile = matrix(data=NA,nrow=10,1)
subjectFile[1,1] = file.path(dataPath, "vertex2.txt")
subjectFile[2,1] = file.path(dataPath, "vertex3.txt")
subjectFile[3,1] = file.path(dataPath, "vertex4.txt")
subjectFile[4,1] = file.path(dataPath, "vertex3.txt")
subjectFile[5,1] = file.path(dataPath, "vertex1.txt")
subjectFile[6,1] = file.path(dataPath, "vertex2.txt")
subjectFile[7,1] = file.path(dataPath, "vertex4.txt")
subjectFile[8,1] = file.path(dataPath, "vertex2.txt")
subjectFile[9,1] = file.path(dataPath, "vertex3.txt")
subjectFile[10,1] = file.path(dataPath, "vertex1.txt")
gftest$testFilesLeft <- (subjectFile)

mm <- verboseRun("vertexMean(gftest$testFilesLeft)",getOption("verbose"))
ma <- verboseRun("vertexApply(gftest$testFilesLeft,mean)",getOption("verbose"))

test_that("vertexApply one output",{
  expect_equal(ma, mm, check.attributes = FALSE) 
})
 
# Need to define global variable when running tests, but normally do not...
testFunc <- function (x) { return(c(1,2))}
ma <- verboseRun("vertexApply(gftest$testFilesLeft, testFunc)", getOption("verbose"))

test_that("vertexApply two output",{
  expect_equal(ma[,1], rep(1, nrow(ma))) 
  expect_equal(ma[,2], rep(2, nrow(ma)))
})

