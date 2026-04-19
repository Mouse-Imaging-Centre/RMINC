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
  expect_equal(rmincAnova[1, 1], rAnova$F[1], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][2], rAnova$Df[2], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][1], rAnova$Df[1], ignore_attr = TRUE)
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
  expect_equal(rmincAnova[1, 1], rAnova$F[1], ignore_attr = TRUE)
  expect_equal(rmincAnova[1, 2], rAnova$F[2], ignore_attr = TRUE)
  expect_equal(rmincAnova[1, 3], rAnova$F[3], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][1], rAnova$Df[1], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][2], rAnova$Df[4], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[2]][1], rAnova$Df[2], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[2]][2], rAnova$Df[4], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[3]][1], rAnova$Df[3], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[3]][2], rAnova$Df[4], ignore_attr = TRUE)
})
rmincAnova = verboseRun(
  "anatAnova(~ Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rAnova = anova(lm(lobeThickness ~ Primary.Diagnosis))


test_that("anatAnova Three Factors", {
  expect_equal(attr(rmincAnova, "df")[[1]][2], rAnova$Df[2], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][1], rAnova$Df[1], ignore_attr = TRUE)
})


rmincAnova = verboseRun(
  "anatAnova(~Age*Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = as.factor(gf$Primary.Diagnosis)
rAnova = anova(lm(lobeThickness ~ Age * Primary.Diagnosis))

test_that("anatAnova Three Factors Interaction", {
  expect_equal(rmincAnova[1, 1], rAnova$F[1], ignore_attr = TRUE)
  expect_equal(rmincAnova[1, 2], rAnova$F[2], ignore_attr = TRUE)
  expect_equal(rmincAnova[1, 3], rAnova$F[3], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][1], rAnova$Df[1], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][2], rAnova$Df[4], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[2]][1], rAnova$Df[2], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[2]][2], rAnova$Df[4], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[3]][1], rAnova$Df[3], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[3]][2], rAnova$Df[4], ignore_attr = TRUE)
})
