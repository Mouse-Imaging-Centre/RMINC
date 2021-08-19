library(testthat)
context("vertexLm")

dataPath="/home/vfonov/src/RMINC/tests/testthat/dataPath"
#if(!exists("dataPath"))
#  dataPath <- tempdir()

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

rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Age,gftest)",getOption("verbose"))
rmincLm2 <- verboseRun("vertexLm(testFilesLeft ~ Age,gftest2,column='thk')",getOption("verbose"))


gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
gftest2$testLeft = t(vertexTable(gftest2$testFilesLeft))

rmod <- lm(testLeft[,1]~Age,gftest)
rLm = summary(rmod)

test_that("vertexLm Two Factors",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})


test_that("vertexLm Two Factors for .csv.gz file",{
	expect_that(rmincLm[1,1],is_equivalent_to(rmincLm2[1,1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rmincLm2[1,2]))
	expect_that(rmincLm[1,3],is_equivalent_to(rmincLm2[1,3]))
	expect_that(rmincLm[1,4],is_equivalent_to(rmincLm2[1,4]))
	expect_that(rmincLm[1,5],is_equivalent_to(rmincLm2[1,5]))
	expect_that(rmincLm[1,6],is_equivalent_to(rmincLm2[1,6]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(attr(rmincLm2,"df")[[2]]))
})


context("writeVertexLm")
writeVertex(rmincLm,file.path(dataPath, "test_lm_no_col_names.txt"),col.names=F,header=F)
writeVertex(rmincLm,file.path(dataPath, "test_lm_with_col_names.txt"),col.names=T,header=F)
writeVertex(rmincLm,file.path(dataPath, "test_lm_with_col_names.csv"),col.names=T,header=F)
writeVertex(rmincLm,file.path(dataPath, "test_lm_with_header.txt"),col.names=F,header=T)
writeVertex(rmincLm,file.path(dataPath, "test_lm_with_col_names_gz.csv.gz"),col.names=T,header=F)

test_that("writeVertexLm", {
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_lm_no_col_names.txt")), "e959ba23a46866f5b9ffafb514c7eb93")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_lm_with_col_names.txt")), "733e1cf1fc8ca06012f7799ae8bc8f06")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_lm_with_col_names.csv")), "1ac7012c8dc8f717a88cc9e2e3f38ba9")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_lm_with_header.txt")), "b79cc94435bffd6a11a8a4b659408d04")
    expect_equivalent(tools::md5sum(file.path(dataPath, "test_lm_with_col_names_gz.csv.gz")), "f294190f168801e74a4e756efd6446d5")
})




test_that("Likelihood and information criteria are computed correctly", {
  expect_equal(as.numeric(rmincLm[1,"logLik"]), as.numeric(logLik(rmod)))
  expect_equal(as.numeric(AIC(rmincLm)[1]), as.numeric(AIC(rmod)))
  expect_equal(as.numeric(BIC(rmincLm)[1]), as.numeric(BIC(rmod)))
})

rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Age*Sex,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Age*Sex,gftest))

test_that("vertexLm Interaction",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[3,1]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[4,1]))
	expect_that(rmincLm[1,7],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,8],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(rmincLm[1,9],is_equivalent_to(rLm$coefficients[3,3]))	
	expect_that(rmincLm[1,10],is_equivalent_to(rLm$coefficients[4,3]))	
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})

rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Group,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Group,gftest))

test_that("vertexLm Three Factors",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[3,1]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,7],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(rmincLm[1,8],is_equivalent_to(rLm$coefficients[3,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})


rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Age*Group,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Age*Group,gftest))

test_that("vertexLm Three Factors Interaction",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[3,1]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[4,1]))
	expect_that(rmincLm[1,7],is_equivalent_to(rLm$coefficients[5,1]))
	expect_that(rmincLm[1,8],is_equivalent_to(rLm$coefficients[6,1]))
	expect_that(rmincLm[1,9],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,10],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(rmincLm[1,11],is_equivalent_to(rLm$coefficients[3,3]))
	expect_that(rmincLm[1,12],is_equivalent_to(rLm$coefficients[4,3]))
	expect_that(rmincLm[1,13],is_equivalent_to(rLm$coefficients[5,3]))
	expect_that(rmincLm[1,14],is_equivalent_to(rLm$coefficients[6,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})
