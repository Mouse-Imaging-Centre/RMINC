requireNamespace("testthat")
context("anatLm")

gf <- read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf <- civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf <- civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)



rmincLm = verboseRun("anatLm(~ Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness~Sex))

test_that("anatLm Two Factors",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})

rmincLm = verboseRun("anatLm(~ Age*Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness~Age*Sex))

test_that("anatLm Interaction",{
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

rmincLm = verboseRun("anatLm(~ Primary.Diagnosis,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness~Primary.Diagnosis))


test_that("anatLm Three Factors",{
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

rmincLm = verboseRun("anatLm(~Primary.Diagnosis*Age,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness~Primary.Diagnosis*Age))

test_that("anatLm Three Factors Interaction",{
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
