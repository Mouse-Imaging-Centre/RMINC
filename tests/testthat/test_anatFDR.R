library(testthat)
context("anatFDR")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "RMINC-test-data-main/rminctestdata/")

gf <- read.csv(file.path(dataPath, "CIVET_TEST.csv"), stringsAsFactors = TRUE)
gf <- civet.getAllFilenames(gf,"ID","TEST",file.path(dataPath, "CIVET"),"TRUE","1.1.12")
gf <- civet.readAllCivetFiles(file.path(dataPath, "AAL.csv"),gf)

rmincLm <- verboseRun("anatLm(~ Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness~Sex))

rLmFDR1 = p.adjust( pt2(rmincLm[,5],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,6],attr(rmincLm,"df")[[3]]),"fdr")


rmincFDR = verboseRun("anatFDR(rmincLm)",getOption("verbose"))


test_that("anatFDR Two Factors",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
})

rmincLm <- verboseRun("anatLm(~ Age*Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness~Age*Sex))

rLmFDR1 = p.adjust( pt2(rmincLm[,7],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,8],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,9],attr(rmincLm,"df")[[4]]),"fdr")
rLmFDR4 = p.adjust( pt2(rmincLm[,10],attr(rmincLm,"df")[[5]]),"fdr")

rmincFDR = verboseRun("anatFDR(rmincLm)",getOption("verbose"))


test_that("anatFDR Interaction",{
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

rmincLm <- verboseRun("anatLm(~ Primary.Diagnosis,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness~Primary.Diagnosis))

rLmFDR1 = p.adjust( pt2(rmincLm[,6],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,7],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,8],attr(rmincLm,"df")[[4]]),"fdr")

rmincFDR = verboseRun("anatFDR(rmincLm)",getOption("verbose"))

test_that("anatFDR Three Factors",{
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

rmincLm <- verboseRun("anatLm(~Primary.Diagnosis*Age,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness~Primary.Diagnosis*Age))

rLmFDR1 = p.adjust( pt2(rmincLm[,9],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,10],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,11],attr(rmincLm,"df")[[4]]),"fdr")
rLmFDR4 = p.adjust( pt2(rmincLm[,12],attr(rmincLm,"df")[[5]]),"fdr")
rLmFDR5 = p.adjust( pt2(rmincLm[,13],attr(rmincLm,"df")[[6]]),"fdr")
rLmFDR6 = p.adjust( pt2(rmincLm[,14],attr(rmincLm,"df")[[7]]),"fdr")


rmincFDR = verboseRun("anatFDR(rmincLm)",getOption("verbose"))

test_that("anatFDR Three Factors Interaction",{
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

test_that("rownames are preserved", {
  expect_equal(rownames(rmincFDR), rownames(rmincLm))
})