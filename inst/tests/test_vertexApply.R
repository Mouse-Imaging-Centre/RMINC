requireNamespace("testthat")
context("vertexApply")

gftest <- read.csv('/tmp/rminctestdata/subject.csv')
subjectFile = matrix(data=NA,nrow=10,1)
subjectFile[1,1]  = '/tmp/rminctestdata/vertex2.txt'
subjectFile[2,1]  = '/tmp/rminctestdata/vertex3.txt'
subjectFile[3,1]  = '/tmp/rminctestdata/vertex4.txt'
subjectFile[4,1]  = '/tmp/rminctestdata/vertex3.txt'
subjectFile[5,1]  = '/tmp/rminctestdata/vertex1.txt'
subjectFile[6,1]  = '/tmp/rminctestdata/vertex2.txt'
subjectFile[7,1]  = '/tmp/rminctestdata/vertex4.txt'
subjectFile[8,1]  = '/tmp/rminctestdata/vertex2.txt'
subjectFile[9,1]  = '/tmp/rminctestdata/vertex3.txt'
subjectFile[10,1] = '/tmp/rminctestdata/vertex1.txt'
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

