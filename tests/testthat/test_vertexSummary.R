library(testthat)
#testthat test script for vertex summary functions
# vertexMean, vertexSd, vertexVar, vertexSum 

if(!exists("dataPath")) 
    dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

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
gftest$testLeft <- t(vertexTable(gftest$testFilesLeft))

# test ability to work with .csv.gz files
# convert txt file into a .csv file with a column thk, and a dummy column 'dummy'
for(i in c("vertex1.txt","vertex2.txt","vertex3.txt","vertex4.txt")) {
    o<-sub('.txt','.csv.gz',i)
    dummy_data<-read.table(file.path(dataPath, i))
    colnames(dummy_data)<-'thk'
    dummy_data$dummy<-100
    readr::write_csv(dummy_data,file.path(dataPath, o))
}

gftest2 <- read.csv(file.path(dataPath, "subject.csv"))
subjectFile2 = matrix(data=NA,nrow=10,1)
subjectFile2[1,1] = file.path(dataPath, "vertex2.csv.gz")
subjectFile2[2,1] = file.path(dataPath, "vertex3.csv.gz")
subjectFile2[3,1] = file.path(dataPath, "vertex4.csv.gz")
subjectFile2[4,1] = file.path(dataPath, "vertex3.csv.gz")
subjectFile2[5,1] = file.path(dataPath, "vertex1.csv.gz")
subjectFile2[6,1] = file.path(dataPath, "vertex2.csv.gz")
subjectFile2[7,1] = file.path(dataPath, "vertex4.csv.gz")
subjectFile2[8,1] = file.path(dataPath, "vertex2.csv.gz")
subjectFile2[9,1] = file.path(dataPath, "vertex3.csv.gz")
subjectFile2[10,1] = file.path(dataPath, "vertex1.csv.gz")
gftest2$testFilesLeft <- (subjectFile2)
gftest2$testLeft <- t(vertexTable(gftest2$testFilesLeft))

context("vertexTable")

table1<-vertexTable(gftest$testFilesLeft ,column=1)
table2<-vertexTable(gftest2$testFilesLeft ,column=1)
test_that("vertexTable", {
    expect_equal(table1,table2)
})


context("vertexMean")

#Calculate mean

vm <- verboseRun("vertexMean(gftest$testFilesLeft)",getOption("verbose"))

test_that("vertexMean", {
    expect_equal(mean(gftest$testLeft[,1]), vm[1])
    expect_equal(mean(gftest$testLeft[,2]), vm[2])
    expect_equal(mean(gftest$testLeft[,3]), vm[3])	
})

vm2 <- verboseRun("vertexMean(gftest2$testFilesLeft)",getOption("verbose"))
test_that("vertexMean.csv", {
    expect_equal(mean(gftest$testLeft[,1]), vm2[1])
    expect_equal(mean(gftest$testLeft[,2]), vm2[2])
    expect_equal(mean(gftest$testLeft[,3]), vm2[3])	
})


context("writeVertexMean")
writeVertex(vm,file.path(dataPath, "test_mean_no_col_names.txt"),col.names=F,header=F)
writeVertex(vm,file.path(dataPath, "test_mean_with_col_names.txt"),col.names=T,header=F)
writeVertex(vm,file.path(dataPath, "test_mean_with_col_names.csv"),col.names=T,header=F)
writeVertex(vm,file.path(dataPath, "test_mean_with_header.txt"),col.names=F,header=T)
writeVertex(vm,file.path(dataPath, "test_mean_with_col_names_gz.csv.gz"),col.names=T,header=F)

test_that("writeVertexMean", {
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_mean_no_col_names.txt")), "7cd5da918fb95bb99ed0a30b5d794ddf")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_mean_with_col_names.txt")), "0eb3b6a54bddc917f0f8e677a2905a57")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_mean_with_col_names.csv")), "0eb3b6a54bddc917f0f8e677a2905a57")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_mean_with_header.txt")), "19f6023459df8f953fe086ec32177485")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_mean_with_col_names_gz.csv.gz")), "2f4e3b87f7bb2a01f0615f1a949823bd")
})



context("vertexSum")

#Calculate sum

vs <- verboseRun("vertexSum(gftest$testFilesLeft)",getOption("verbose"))

test_that("vertexSum", {
    expect_equal(sum(gftest$testLeft[,1]), vs[1])
    expect_equal(sum(gftest$testLeft[,2]), vs[2])
    expect_equal(sum(gftest$testLeft[,3]), vs[3])	
})

vs2 <- verboseRun("vertexSum(gftest2$testFilesLeft)",getOption("verbose"))
test_that("vertexSum.csv", {
    expect_equal(sum(gftest$testLeft[,1]), vs2[1])
    expect_equal(sum(gftest$testLeft[,2]), vs2[2])
    expect_equal(sum(gftest$testLeft[,3]), vs2[3])	
})

context("vertexVar")

#Calculate variance
vv <- verboseRun("vertexVar(gftest$testFilesLeft)",getOption("verbose"))

test_that("vertexVar", {
    expect_equal(var(gftest$testLeft[,1]), vv[1])
    expect_equal(var(gftest$testLeft[,2]), vv[2])
    expect_equal(var(gftest$testLeft[,3]), vv[3])	
})
vv2 <- verboseRun("vertexVar(gftest2$testFilesLeft)",getOption("verbose"))

test_that("vertexVar", {
    expect_equal(var(gftest$testLeft[,1]), vv2[1])
    expect_equal(var(gftest$testLeft[,2]), vv2[2])
    expect_equal(var(gftest$testLeft[,3]), vv2[3])	
})


context("vertexSd")

#Calculate standard deviation
vsd <- verboseRun("vertexSd(gftest$testFilesLeft)",getOption("verbose"))

test_that("vertexSd", {
    expect_equal(sd(gftest$testLeft[,1]), vsd[1])
    expect_equal(sd(gftest$testLeft[,2]), vsd[2])
    expect_equal(sd(gftest$testLeft[,3]), vsd[3])	
})
#Calculate standard deviation

vsd2 <- verboseRun("vertexSd(gftest2$testFilesLeft)",getOption("verbose"))

test_that("vertexSd", {
    expect_equal(sd(gftest$testLeft[,1]), vsd2[1])
    expect_equal(sd(gftest$testLeft[,2]), vsd2[2])
    expect_equal(sd(gftest$testLeft[,3]), vsd2[3])	
})


# Test saving results of the vertex summary into file