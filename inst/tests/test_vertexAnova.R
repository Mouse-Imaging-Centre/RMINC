library("testthat")
context("vertexAnova")

gftest <- read.csv('/tmp/rminctestdata/subject.csv')
subjectFile = matrix(data=NA,nrow=10,1)
subjectFile[1,1] = '/tmp/rminctestdata/vertex2.txt'
subjectFile[2,1] = '/tmp/rminctestdata/vertex3.txt'
subjectFile[3,1] = '/tmp/rminctestdata/vertex4.txt'
subjectFile[4,1] = '/tmp/rminctestdata/vertex3.txt'
subjectFile[5,1] = '/tmp/rminctestdata/vertex1.txt'
subjectFile[6,1] = '/tmp/rminctestdata/vertex2.txt'
subjectFile[7,1] = '/tmp/rminctestdata/vertex4.txt'
subjectFile[8,1] = '/tmp/rminctestdata/vertex2.txt'
subjectFile[9,1] = '/tmp/rminctestdata/vertex3.txt'
subjectFile[10,1] = '/tmp/rminctestdata/vertex1.txt'
gftest$testFilesLeft <- (subjectFile)


rmincAnova <- verboseRun("vertexAnova(testFilesLeft ~ Sex,gftest)",getOption("verbose"))
gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[,1]~Sex,gftest))

test_that("vertexAnova Two Factors",{
	expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
})

rmincAnova <- verboseRun("vertexAnova(testFilesLeft ~ Age*Sex,gftest)",getOption("verbose"))
gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[,1]~Age*Sex,gftest))

test_that("vertexAnova Interaction",{
	expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
	expect_that(rmincAnova[1,2],is_equivalent_to(rAnova$F[2]))
	expect_that(rmincAnova[1,3],is_equivalent_to(rAnova$F[3]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[4]))
	expect_that(attr(rmincAnova,"df")[[2]][1],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[2]][2],is_equivalent_to(rAnova$Df[4]))
	expect_that(attr(rmincAnova,"df")[[3]][1],is_equivalent_to(rAnova$Df[3]))
	expect_that(attr(rmincAnova,"df")[[3]][2],is_equivalent_to(rAnova$Df[4]))
})

rmincAnova <- verboseRun("vertexAnova(testFilesLeft ~ Group,gftest)",getOption("verbose"))
gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[,1]~Group,gftest))

test_that("vertexAnova Three Factors",{
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
})
rmincAnova <- verboseRun("vertexAnova(testFilesLeft ~ Age*Group,gftest)",getOption("verbose"))

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[,1]~Age*Group,gftest))

test_that("vertexAnova Three Factors Interaction",{
         expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
         expect_that(rmincAnova[1,2],is_equivalent_to(rAnova$F[2]))
         expect_that(rmincAnova[1,3],is_equivalent_to(rAnova$F[3]))
         expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
         expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[4]))
         expect_that(attr(rmincAnova,"df")[[2]][1],is_equivalent_to(rAnova$Df[2]))
         expect_that(attr(rmincAnova,"df")[[2]][2],is_equivalent_to(rAnova$Df[4]))
         expect_that(attr(rmincAnova,"df")[[3]][1],is_equivalent_to(rAnova$Df[3]))
         expect_that(attr(rmincAnova,"df")[[3]][2],is_equivalent_to(rAnova$Df[4])) 
})





