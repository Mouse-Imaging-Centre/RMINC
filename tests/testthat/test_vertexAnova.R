library(testthat)

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gftest <- read.csv(file.path(dataPath, "subject.csv"), stringsAsFactors = TRUE)

subjectFile = matrix(data = NA, nrow = 10, 1)
subjectFile[1, 1] = file.path(dataPath, "vertex2.txt")
subjectFile[2, 1] = file.path(dataPath, "vertex3.txt")
subjectFile[3, 1] = file.path(dataPath, "vertex4.txt")
subjectFile[4, 1] = file.path(dataPath, "vertex3.txt")
subjectFile[5, 1] = file.path(dataPath, "vertex1.txt")
subjectFile[6, 1] = file.path(dataPath, "vertex2.txt")
subjectFile[7, 1] = file.path(dataPath, "vertex4.txt")
subjectFile[8, 1] = file.path(dataPath, "vertex2.txt")
subjectFile[9, 1] = file.path(dataPath, "vertex3.txt")
subjectFile[10, 1] = file.path(dataPath, "vertex1.txt")
gftest$testFilesLeft <- (subjectFile)


rmincAnova <- verboseRun(
  "vertexAnova(testFilesLeft ~ Sex,gftest)",
  getOption("verbose")
)
gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[, 1] ~ Sex, gftest))

test_that("vertexAnova Two Factors", {
  expect_equal(rmincAnova[1, 1], rAnova$F[1], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][2], rAnova$Df[2], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][1], rAnova$Df[1], ignore_attr = TRUE)
})

rmincAnova <- verboseRun(
  "vertexAnova(testFilesLeft ~ Age*Sex,gftest)",
  getOption("verbose")
)
gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[, 1] ~ Age * Sex, gftest))

test_that("vertexAnova Interaction", {
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

rmincAnova <- verboseRun(
  "vertexAnova(testFilesLeft ~ Group,gftest)",
  getOption("verbose")
)
gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[, 1] ~ Group, gftest))

test_that("vertexAnova Three Factors", {
  expect_equal(attr(rmincAnova, "df")[[1]][2], rAnova$Df[2], ignore_attr = TRUE)
  expect_equal(attr(rmincAnova, "df")[[1]][1], rAnova$Df[1], ignore_attr = TRUE)
})
rmincAnova <- verboseRun(
  "vertexAnova(testFilesLeft ~ Age*Group,gftest)",
  getOption("verbose")
)

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rAnova = anova(lm(testLeft[, 1] ~ Age * Group, gftest))

test_that("vertexAnova Three Factors Interaction", {
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
