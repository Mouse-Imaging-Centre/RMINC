library(testthat)

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "CIVET_TEST.csv"), stringsAsFactors = TRUE)
gf <- civet.getAllFilenames(
  gf,
  "ID",
  "TEST",
  file.path(dataPath, "CIVET"),
  "TRUE",
  "1.1.12"
)
gf <- civet.readAllCivetFiles(file.path(dataPath, "AAL.csv"), gf)

rmincLm = verboseRun("anatLm(~ Sex,gf,gf$lobeThickness)", getOption("verbose"))

lobeThickness = gf$lobeThickness[, 1]
Age = gf$Age
Sex = gf$Sex
rmod <- lm(lobeThickness ~ Sex)
rLm = summary(rmod)

test_that("anatLm Two Factors", {
  expect_equal(rmincLm[1, 1], rLm$fstatistic[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 2], rLm$r.squared[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 3], rLm$coefficients[1, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 4], rLm$coefficients[2, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 5], rLm$coefficients[1, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 6], rLm$coefficients[2, 3], ignore_attr = TRUE)
  expect_equal(attr(rmincLm, "df")[[2]], rLm$df[2], ignore_attr = TRUE)
})

test_that("Likelihood and information criteria are computed correctly", {
  expect_equal(as.numeric(rmincLm[1, "logLik"]), as.numeric(logLik(rmod)))
  expect_equal(as.numeric(AIC(rmincLm)[1]), as.numeric(AIC(rmod)))
  expect_equal(as.numeric(BIC(rmincLm)[1]), as.numeric(BIC(rmod)))
})

rmincLm = verboseRun(
  "anatLm(~ Age*Sex,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness ~ Age * Sex))

test_that("anatLm Interaction", {
  expect_equal(rmincLm[1, 1], rLm$fstatistic[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 2], rLm$r.squared[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 3], rLm$coefficients[1, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 4], rLm$coefficients[2, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 5], rLm$coefficients[3, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 6], rLm$coefficients[4, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 7], rLm$coefficients[1, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 8], rLm$coefficients[2, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 9], rLm$coefficients[3, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 10], rLm$coefficients[4, 3], ignore_attr = TRUE)
  expect_equal(attr(rmincLm, "df")[[2]], rLm$df[2], ignore_attr = TRUE)
})

rmincLm = verboseRun(
  "anatLm(~ Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness ~ Primary.Diagnosis))


test_that("anatLm Three Factors", {
  expect_equal(rmincLm[1, 1], rLm$fstatistic[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 2], rLm$r.squared[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 3], rLm$coefficients[1, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 4], rLm$coefficients[2, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 5], rLm$coefficients[3, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 6], rLm$coefficients[1, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 7], rLm$coefficients[2, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 8], rLm$coefficients[3, 3], ignore_attr = TRUE)
  expect_equal(attr(rmincLm, "df")[[2]], rLm$df[2], ignore_attr = TRUE)
})

rmincLm = verboseRun(
  "anatLm(~Primary.Diagnosis*Age,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness ~ Primary.Diagnosis * Age))

test_that("anatLm Three Factors Interaction", {
  expect_equal(rmincLm[1, 1], rLm$fstatistic[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 2], rLm$r.squared[1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 3], rLm$coefficients[1, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 4], rLm$coefficients[2, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 5], rLm$coefficients[3, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 6], rLm$coefficients[4, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 7], rLm$coefficients[5, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 8], rLm$coefficients[6, 1], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 9], rLm$coefficients[1, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 10], rLm$coefficients[2, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 11], rLm$coefficients[3, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 12], rLm$coefficients[4, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 13], rLm$coefficients[5, 3], ignore_attr = TRUE)
  expect_equal(rmincLm[1, 14], rLm$coefficients[6, 3], ignore_attr = TRUE)
  expect_equal(attr(rmincLm, "df")[[2]], rLm$df[2], ignore_attr = TRUE)
})


test_that("Weighted anatLm works", {
  w <- runif(20)
  y <- matrix(rnorm(60), ncol = 3)
  x <- data.frame(a = runif(20), b = rnorm(20), c = rgamma(20, 1))

  verboseRun(alm <- anatLm(~ a + b + c, data = x, anat = y, w = w))
  lmods <- apply(y, 2, function(col) lm(col ~ a + b + c, data = x, weights = w))
  expect_equal(
    as.numeric(t(sapply(lmods, function(m) {
      summary(m)$coefficients[, "t value"]
    }))),
    as.numeric(alm[, grepl("tvalue", colnames(alm))]),
    ignore_attr = TRUE
  )
  expect_equal(
    as.numeric(t(sapply(lmods, function(m) summary(m)$r.squared))),
    as.numeric(alm[, "R-squared"]),
    ignore_attr = TRUE
  )
  expect_equal(
    as.numeric(t(sapply(lmods, function(m) summary(m)$fstatistic["value"]))),
    as.numeric(alm[, "F-statistic"]),
    ignore_attr = TRUE
  )
  expect_equal(
    as.numeric(t(sapply(lmods, coef))),
    as.numeric(alm[, grepl("beta", colnames(alm))]),
    ignore_attr = TRUE
  )
  expect_equal(sapply(lmods, logLik), alm[, "logLik"], ignore_attr = TRUE)
  expect_equal(sapply(lmods, AIC), AIC(alm), ignore_attr = TRUE)
})


test_that("Test that passing a matrix on the RHS fails", {
  d <- data.frame(y = rnorm(10))
  d$x <- scale(rnorm(10))

  expect_error(anatLm(~x, data = d, anat = as.matrix(d$y)), "matrix")
})
