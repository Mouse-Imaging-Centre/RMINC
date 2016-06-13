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
ma <- verboseRun("vertexApply(gftest$testFilesLeft,quote(mean(x)))",getOption("verbose"))

test_that("vertexApply one output",{
     for (nVox in 1:length(mm)) {
           expect_equal(ma[nVox], mm[nVox]) }
 })
 
# Need to define global variable when running tests, but normally do not...
testFunc <<- function (x) { return(c(1,2))}
ma <- verboseRun("vertexApply(gftest$testFilesLeft,quote(testFunc(x)))",getOption("verbose"))

test_that("vertexApply two output",{
    for (nVox in 1:dim(ma)[1]) {
 	   expect_equal(ma[nVox,1], 1) 
	   expect_equal(ma[nVox,2], 2)}
})

