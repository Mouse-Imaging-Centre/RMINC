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

#rmincLm <- verboseRun("vertexLm(testFilesLeft ~ Sex,gftest) ",getOption("verbose"))
rmincLm <- vertexLm(testFilesLeft ~ Sex, gftest, column = 1)

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[, 1] ~ Sex, gftest))
rLmFDR1 = p.adjust(pt2(rmincLm[, 5], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 6], attr(rmincLm, "df")[[3]]), "fdr")

rmincFDR <- verboseRun("vertexFDR(rmincLm)", getOption("verbose"))

test_that("vertexFDR Two Factors", {
  expect_equal(rLmFDR1[1], rmincFDR[1, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[2], rmincFDR[2, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[3], rmincFDR[3, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR2[1], rmincFDR[1, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[2], rmincFDR[2, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[3], rmincFDR[3, 3], ignore_attr = TRUE)
})
rmincLm <- verboseRun(
  "vertexLm(testFilesLeft ~ Age*Sex,gftest)",
  getOption("verbose")
)

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[, 1] ~ Age * Sex, gftest))

rLmFDR1 = p.adjust(pt2(rmincLm[, 7], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 8], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 9], attr(rmincLm, "df")[[4]]), "fdr")
rLmFDR4 = p.adjust(pt2(rmincLm[, 10], attr(rmincLm, "df")[[5]]), "fdr")

rmincFDR <- verboseRun("vertexFDR(rmincLm)", getOption("verbose"))


test_that("vertexFDR Interaction", {
  expect_equal(rLmFDR1[1], rmincFDR[1, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[2], rmincFDR[2, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[3], rmincFDR[3, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR2[1], rmincFDR[1, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[2], rmincFDR[2, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[3], rmincFDR[3, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR3[1], rmincFDR[1, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR3[2], rmincFDR[2, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR3[3], rmincFDR[3, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR4[1], rmincFDR[1, 5], ignore_attr = TRUE)
  expect_equal(rLmFDR4[2], rmincFDR[2, 5], ignore_attr = TRUE)
  expect_equal(rLmFDR4[3], rmincFDR[3, 5], ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "vertexLm(testFilesLeft ~ Group,gftest)",
  getOption("verbose")
)

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[, 1] ~ Group, gftest))

rLmFDR1 = p.adjust(pt2(rmincLm[, 6], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 7], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 8], attr(rmincLm, "df")[[4]]), "fdr")

rmincFDR <- verboseRun("vertexFDR(rmincLm)", getOption("verbose"))

test_that("vertexFDR Three Factors", {
  expect_equal(rLmFDR1[1], rmincFDR[1, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[2], rmincFDR[2, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[3], rmincFDR[3, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR2[1], rmincFDR[1, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[2], rmincFDR[2, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[3], rmincFDR[3, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR3[1], rmincFDR[1, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR3[2], rmincFDR[2, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR3[3], rmincFDR[3, 4], ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "vertexLm(testFilesLeft ~ Age*Group,gftest)",
  getOption("verbose")
)

gftest$testLeft = t(vertexTable(gftest$testFilesLeft))
rLm = summary(lm(testLeft[, 1] ~ Age * Group, gftest))

rLmFDR1 = p.adjust(pt2(rmincLm[, 9], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 10], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 11], attr(rmincLm, "df")[[4]]), "fdr")
rLmFDR4 = p.adjust(pt2(rmincLm[, 12], attr(rmincLm, "df")[[5]]), "fdr")
rLmFDR5 = p.adjust(pt2(rmincLm[, 13], attr(rmincLm, "df")[[6]]), "fdr")
rLmFDR6 = p.adjust(pt2(rmincLm[, 14], attr(rmincLm, "df")[[7]]), "fdr")


rmincFDR <- verboseRun("vertexFDR(rmincLm)", getOption("verbose"))

test_that("vertexLm Three Factors Interaction", {
  expect_equal(rLmFDR1[1], rmincFDR[1, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[2], rmincFDR[2, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR1[3], rmincFDR[3, 2], ignore_attr = TRUE)
  expect_equal(rLmFDR2[1], rmincFDR[1, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[2], rmincFDR[2, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR2[3], rmincFDR[3, 3], ignore_attr = TRUE)
  expect_equal(rLmFDR3[1], rmincFDR[1, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR3[2], rmincFDR[2, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR3[3], rmincFDR[3, 4], ignore_attr = TRUE)
  expect_equal(rLmFDR4[1], rmincFDR[1, 5], ignore_attr = TRUE)
  expect_equal(rLmFDR4[2], rmincFDR[2, 5], ignore_attr = TRUE)
  expect_equal(rLmFDR4[3], rmincFDR[3, 5], ignore_attr = TRUE)
  expect_equal(rLmFDR5[1], rmincFDR[1, 6], ignore_attr = TRUE)
  expect_equal(rLmFDR5[2], rmincFDR[2, 6], ignore_attr = TRUE)
  expect_equal(rLmFDR5[3], rmincFDR[3, 6], ignore_attr = TRUE)
  expect_equal(rLmFDR6[1], rmincFDR[1, 7], ignore_attr = TRUE)
  expect_equal(rLmFDR6[2], rmincFDR[2, 7], ignore_attr = TRUE)
  expect_equal(rLmFDR6[3], rmincFDR[3, 7], ignore_attr = TRUE)
})
