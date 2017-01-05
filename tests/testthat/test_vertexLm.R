library(testthat)
context("vertexLm")

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

rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Age,gftest)",getOption("verbose"))


gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
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

test_that("Likelihood and information criteria are computed correctly", {
  expect_equal(as.numeric(rmincLm[1,"logLik"]), as.numeric(logLik(rmod)))
  expect_equal(as.numeric(AIC(rmincLm)[1]), as.numeric(AIC(rmod)))
  expect_equal(as.numeric(BIC(rmincLm)[1]), as.numeric(BIC(rmod)))
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
