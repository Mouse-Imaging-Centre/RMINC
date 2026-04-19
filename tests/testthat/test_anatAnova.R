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

rmincAnova = verboseRun(
  "anatAnova(~ Sex,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Age = gf$Age
Sex = gf$Sex
rAnova = anova(lm(lobeThickness ~ Sex))

test_that("anatAnova Two Factors", {
  expect_equal(unname(rmincAnova[1, 1]), unname(rAnova$F[1]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][2]), unname(rAnova$Df[2]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][1]), unname(rAnova$Df[1]))
})

rmincAnova = verboseRun(
  "anatAnova(~ Age*Sex,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Age = gf$Age
Sex = gf$Sex
rAnova = anova(lm(lobeThickness ~ Age * Sex))

test_that("anatAnova Interaction", {
  expect_equal(unname(rmincAnova[1, 1]), unname(rAnova$F[1]))
  expect_equal(unname(rmincAnova[1, 2]), unname(rAnova$F[2]))
  expect_equal(unname(rmincAnova[1, 3]), unname(rAnova$F[3]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][1]), unname(rAnova$Df[1]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][2]), unname(rAnova$Df[4]))
  expect_equal(unname(attr(rmincAnova, "df")[[2]][1]), unname(rAnova$Df[2]))
  expect_equal(unname(attr(rmincAnova, "df")[[2]][2]), unname(rAnova$Df[4]))
  expect_equal(unname(attr(rmincAnova, "df")[[3]][1]), unname(rAnova$Df[3]))
  expect_equal(unname(attr(rmincAnova, "df")[[3]][2]), unname(rAnova$Df[4]))
})
rmincAnova = verboseRun(
  "anatAnova(~ Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rAnova = anova(lm(lobeThickness ~ Primary.Diagnosis))


test_that("anatAnova Three Factors", {
  expect_equal(unname(attr(rmincAnova, "df")[[1]][2]), unname(rAnova$Df[2]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][1]), unname(rAnova$Df[1]))
})


rmincAnova = verboseRun(
  "anatAnova(~Age*Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = as.factor(gf$Primary.Diagnosis)
rAnova = anova(lm(lobeThickness ~ Age * Primary.Diagnosis))

test_that("anatAnova Three Factors Interaction", {
  expect_equal(unname(rmincAnova[1, 1]), unname(rAnova$F[1]))
  expect_equal(unname(rmincAnova[1, 2]), unname(rAnova$F[2]))
  expect_equal(unname(rmincAnova[1, 3]), unname(rAnova$F[3]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][1]), unname(rAnova$Df[1]))
  expect_equal(unname(attr(rmincAnova, "df")[[1]][2]), unname(rAnova$Df[4]))
  expect_equal(unname(attr(rmincAnova, "df")[[2]][1]), unname(rAnova$Df[2]))
  expect_equal(unname(attr(rmincAnova, "df")[[2]][2]), unname(rAnova$Df[4]))
  expect_equal(unname(attr(rmincAnova, "df")[[3]][1]), unname(rAnova$Df[3]))
  expect_equal(unname(attr(rmincAnova, "df")[[3]][2]), unname(rAnova$Df[4]))
})
