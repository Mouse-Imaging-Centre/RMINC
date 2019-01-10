library(testthat)
context("vertexFDR")

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
rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Sex,gftest) ",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Sex,gftest))
rLmFDR1 = p.adjust( pt2(rmincLm[,5],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,6],attr(rmincLm,"df")[[3]]),"fdr")

rmincFDR <- verboseRun("vertexFDR(rmincLm)",getOption("verbose"))

test_that("vertexFDR Two Factors",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
})
rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Age*Sex,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Age*Sex,gftest))

rLmFDR1 = p.adjust( pt2(rmincLm[,7],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,8],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,9],attr(rmincLm,"df")[[4]]),"fdr")
rLmFDR4 = p.adjust( pt2(rmincLm[,10],attr(rmincLm,"df")[[5]]),"fdr")

rmincFDR <- verboseRun("vertexFDR(rmincLm)",getOption("verbose"))


test_that("vertexFDR Interaction",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
	expect_that(rLmFDR3[1],is_equivalent_to(rmincFDR[1,4]))
	expect_that(rLmFDR3[2],is_equivalent_to(rmincFDR[2,4]))
	expect_that(rLmFDR3[3],is_equivalent_to(rmincFDR[3,4]))
	expect_that(rLmFDR4[1],is_equivalent_to(rmincFDR[1,5]))
	expect_that(rLmFDR4[2],is_equivalent_to(rmincFDR[2,5]))
	expect_that(rLmFDR4[3],is_equivalent_to(rmincFDR[3,5]))
})

rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Group,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Group,gftest))

rLmFDR1 = p.adjust( pt2(rmincLm[,6],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,7],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,8],attr(rmincLm,"df")[[4]]),"fdr")

rmincFDR <- verboseRun("vertexFDR(rmincLm)",getOption("verbose"))

test_that("vertexFDR Three Factors",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
	expect_that(rLmFDR3[1],is_equivalent_to(rmincFDR[1,4]))
	expect_that(rLmFDR3[2],is_equivalent_to(rmincFDR[2,4]))
	expect_that(rLmFDR3[3],is_equivalent_to(rmincFDR[3,4]))
})

rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Age*Group,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[,1]~Age*Group,gftest))

rLmFDR1 = p.adjust( pt2(rmincLm[,9],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,10],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,11],attr(rmincLm,"df")[[4]]),"fdr")
rLmFDR4 = p.adjust( pt2(rmincLm[,12],attr(rmincLm,"df")[[5]]),"fdr")
rLmFDR5 = p.adjust( pt2(rmincLm[,13],attr(rmincLm,"df")[[6]]),"fdr")
rLmFDR6 = p.adjust( pt2(rmincLm[,14],attr(rmincLm,"df")[[7]]),"fdr")


rmincFDR <- verboseRun("vertexFDR(rmincLm)",getOption("verbose"))

test_that("vertexLm Three Factors Interaction",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
	expect_that(rLmFDR3[1],is_equivalent_to(rmincFDR[1,4]))
	expect_that(rLmFDR3[2],is_equivalent_to(rmincFDR[2,4]))
	expect_that(rLmFDR3[3],is_equivalent_to(rmincFDR[3,4]))
	expect_that(rLmFDR4[1],is_equivalent_to(rmincFDR[1,5]))
	expect_that(rLmFDR4[2],is_equivalent_to(rmincFDR[2,5]))
	expect_that(rLmFDR4[3],is_equivalent_to(rmincFDR[3,5]))
	expect_that(rLmFDR5[1],is_equivalent_to(rmincFDR[1,6]))
	expect_that(rLmFDR5[2],is_equivalent_to(rmincFDR[2,6]))
	expect_that(rLmFDR5[3],is_equivalent_to(rmincFDR[3,6]))
	expect_that(rLmFDR6[1],is_equivalent_to(rmincFDR[1,7]))
	expect_that(rLmFDR6[2],is_equivalent_to(rmincFDR[2,7]))
	expect_that(rLmFDR6[3],is_equivalent_to(rmincFDR[3,7]))
})
