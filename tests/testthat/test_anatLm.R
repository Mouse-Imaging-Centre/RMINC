library(testthat)
context("anatLm")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "CIVET_TEST.csv"))
gf <- civet.getAllFilenames(gf,"ID","TEST",file.path(dataPath, "CIVET"),"TRUE","1.1.12")
gf <- civet.readAllCivetFiles(file.path(dataPath, "AAL.csv"),gf)

rmincLm = verboseRun("anatLm(~ Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rmod <- lm(lobeThickness~Sex)
rLm = summary(rmod)

test_that("anatLm Two Factors",{
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

context("weighted anatLm")

test_that("Weighted anatLm works", {
  w <- runif(20)
  y <- matrix(rnorm(60), ncol = 3)
  x <- data.frame(a = runif(20), b = rnorm(20), c = rgamma(20, 1))
  
  verboseRun(alm <- anatLm(~ a + b + c, data = x, anat = y, w = w))
  lmods <- apply(y, 2, function(col) lm(col ~ a + b + c, data = x, weights = w))
  expect_equivalent(as.numeric(t(sapply(lmods, function(m) summary(m)$coefficients[, "t value"])))
                    , as.numeric(alm[, grepl("tvalue", colnames(alm))]))
  expect_equivalent(as.numeric(t(sapply(lmods, function(m) summary(m)$r.squared)))
                    , as.numeric(alm[, "R-squared"]))
  expect_equivalent(as.numeric(t(sapply(lmods, function(m) summary(m)$fstatistic["value"])))
                    , as.numeric(alm[, "F-statistic"]))
  expect_equivalent(as.numeric(t(sapply(lmods, coef))), as.numeric(alm[, grepl("beta", colnames(alm))]))
  expect_equivalent(sapply(lmods, logLik), alm[,"logLik"])
  expect_equivalent(sapply(lmods, AIC), AIC(alm))
})

