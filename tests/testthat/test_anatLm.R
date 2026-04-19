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
  expect_equal(unname(rmincLm[1, 1]), unname(rLm$fstatistic[1]))
  expect_equal(unname(rmincLm[1, 2]), unname(rLm$r.squared[1]))
  expect_equal(unname(rmincLm[1, 3]), unname(rLm$coefficients[1, 1]))
  expect_equal(unname(rmincLm[1, 4]), unname(rLm$coefficients[2, 1]))
  expect_equal(unname(rmincLm[1, 5]), unname(rLm$coefficients[1, 3]))
  expect_equal(unname(rmincLm[1, 6]), unname(rLm$coefficients[2, 3]))
  expect_equal(unname(attr(rmincLm, "df")[[2]]), unname(rLm$df[2]))
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
  expect_equal(unname(rmincLm[1, 1]), unname(rLm$fstatistic[1]))
  expect_equal(unname(rmincLm[1, 2]), unname(rLm$r.squared[1]))
  expect_equal(unname(rmincLm[1, 3]), unname(rLm$coefficients[1, 1]))
  expect_equal(unname(rmincLm[1, 4]), unname(rLm$coefficients[2, 1]))
  expect_equal(unname(rmincLm[1, 5]), unname(rLm$coefficients[3, 1]))
  expect_equal(unname(rmincLm[1, 6]), unname(rLm$coefficients[4, 1]))
  expect_equal(unname(rmincLm[1, 7]), unname(rLm$coefficients[1, 3]))
  expect_equal(unname(rmincLm[1, 8]), unname(rLm$coefficients[2, 3]))
  expect_equal(unname(rmincLm[1, 9]), unname(rLm$coefficients[3, 3]))
  expect_equal(unname(rmincLm[1, 10]), unname(rLm$coefficients[4, 3]))
  expect_equal(unname(attr(rmincLm, "df")[[2]]), unname(rLm$df[2]))
})

rmincLm = verboseRun(
  "anatLm(~ Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness ~ Primary.Diagnosis))


test_that("anatLm Three Factors", {
  expect_equal(unname(rmincLm[1, 1]), unname(rLm$fstatistic[1]))
  expect_equal(unname(rmincLm[1, 2]), unname(rLm$r.squared[1]))
  expect_equal(unname(rmincLm[1, 3]), unname(rLm$coefficients[1, 1]))
  expect_equal(unname(rmincLm[1, 4]), unname(rLm$coefficients[2, 1]))
  expect_equal(unname(rmincLm[1, 5]), unname(rLm$coefficients[3, 1]))
  expect_equal(unname(rmincLm[1, 6]), unname(rLm$coefficients[1, 3]))
  expect_equal(unname(rmincLm[1, 7]), unname(rLm$coefficients[2, 3]))
  expect_equal(unname(rmincLm[1, 8]), unname(rLm$coefficients[3, 3]))
  expect_equal(unname(attr(rmincLm, "df")[[2]]), unname(rLm$df[2]))
})

rmincLm = verboseRun(
  "anatLm(~Primary.Diagnosis*Age,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness ~ Primary.Diagnosis * Age))

test_that("anatLm Three Factors Interaction", {
  expect_equal(unname(rmincLm[1, 1]), unname(rLm$fstatistic[1]))
  expect_equal(unname(rmincLm[1, 2]), unname(rLm$r.squared[1]))
  expect_equal(unname(rmincLm[1, 3]), unname(rLm$coefficients[1, 1]))
  expect_equal(unname(rmincLm[1, 4]), unname(rLm$coefficients[2, 1]))
  expect_equal(unname(rmincLm[1, 5]), unname(rLm$coefficients[3, 1]))
  expect_equal(unname(rmincLm[1, 6]), unname(rLm$coefficients[4, 1]))
  expect_equal(unname(rmincLm[1, 7]), unname(rLm$coefficients[5, 1]))
  expect_equal(unname(rmincLm[1, 8]), unname(rLm$coefficients[6, 1]))
  expect_equal(unname(rmincLm[1, 9]), unname(rLm$coefficients[1, 3]))
  expect_equal(unname(rmincLm[1, 10]), unname(rLm$coefficients[2, 3]))
  expect_equal(unname(rmincLm[1, 11]), unname(rLm$coefficients[3, 3]))
  expect_equal(unname(rmincLm[1, 12]), unname(rLm$coefficients[4, 3]))
  expect_equal(unname(rmincLm[1, 13]), unname(rLm$coefficients[5, 3]))
  expect_equal(unname(rmincLm[1, 14]), unname(rLm$coefficients[6, 3]))
  expect_equal(unname(attr(rmincLm, "df")[[2]]), unname(rLm$df[2]))
})


test_that("Weighted anatLm works", {
  w <- runif(20)
  y <- matrix(rnorm(60), ncol = 3)
  x <- data.frame(a = runif(20), b = rnorm(20), c = rgamma(20, 1))

  verboseRun(alm <- anatLm(~ a + b + c, data = x, anat = y, w = w))
  lmods <- apply(y, 2, function(col) lm(col ~ a + b + c, data = x, weights = w))
  expect_equal(unname(as.numeric(t(sapply(lmods, function(m) {
      summary(m)$coefficients[, "t value"]
    })))), unname(as.numeric(alm[, grepl("tvalue", colnames(alm))])))
  expect_equal(unname(as.numeric(t(sapply(lmods, function(m) summary(m)$r.squared)))), unname(as.numeric(alm[, "R-squared"])))
  expect_equal(unname(as.numeric(t(sapply(lmods, function(m) summary(m)$fstatistic["value"])))), unname(as.numeric(alm[, "F-statistic"])))
  expect_equal(unname(as.numeric(t(sapply(lmods, coef)))), unname(as.numeric(alm[, grepl("beta", colnames(alm))])))
  expect_equal(unname(sapply(lmods, logLik)), unname(alm[, "logLik"]))
  expect_equal(unname(sapply(lmods, AIC)), unname(AIC(alm)))
})


test_that("Test that passing a matrix on the RHS fails", {
  d <- data.frame(y = rnorm(10))
  d$x <- scale(rnorm(10))

  expect_error(anatLm(~x, data = d, anat = as.matrix(d$y)), "matrix")
})
