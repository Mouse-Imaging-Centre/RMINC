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

rmincLm <- verboseRun("anatLm(~ Sex,gf,gf$lobeThickness)", getOption("verbose"))

lobeThickness = gf$lobeThickness[, 1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness ~ Sex))

rLmFDR1 = p.adjust(pt2(rmincLm[, 5], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 6], attr(rmincLm, "df")[[3]]), "fdr")


rmincFDR = verboseRun("anatFDR(rmincLm)", getOption("verbose"))


test_that("anatFDR Two Factors", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "anatLm(~ Age*Sex,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Age = gf$Age
Sex = gf$Sex
rLm = summary(lm(lobeThickness ~ Age * Sex))

rLmFDR1 = p.adjust(pt2(rmincLm[, 7], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 8], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 9], attr(rmincLm, "df")[[4]]), "fdr")
rLmFDR4 = p.adjust(pt2(rmincLm[, 10], attr(rmincLm, "df")[[5]]), "fdr")

rmincFDR = verboseRun("anatFDR(rmincLm)", getOption("verbose"))


test_that("anatFDR Interaction", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[1]), unname(rmincFDR[1, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[2]), unname(rmincFDR[2, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[3]), unname(rmincFDR[3, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[1]), unname(rmincFDR[1, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[2]), unname(rmincFDR[2, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[3]), unname(rmincFDR[3, 5]), ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "anatLm(~ Primary.Diagnosis,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness ~ Primary.Diagnosis))

rLmFDR1 = p.adjust(pt2(rmincLm[, 6], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 7], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 8], attr(rmincLm, "df")[[4]]), "fdr")

rmincFDR = verboseRun("anatFDR(rmincLm)", getOption("verbose"))

test_that("anatFDR Three Factors", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[1]), unname(rmincFDR[1, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[2]), unname(rmincFDR[2, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[3]), unname(rmincFDR[3, 4]), ignore_attr = TRUE)
})

rmincLm <- verboseRun(
  "anatLm(~Primary.Diagnosis*Age,gf,gf$lobeThickness)",
  getOption("verbose")
)

lobeThickness = gf$lobeThickness[, 1]
Primary.Diagnosis = gf$Primary.Diagnosis
rLm = summary(lm(lobeThickness ~ Primary.Diagnosis * Age))

rLmFDR1 = p.adjust(pt2(rmincLm[, 9], attr(rmincLm, "df")[[2]]), "fdr")
rLmFDR2 = p.adjust(pt2(rmincLm[, 10], attr(rmincLm, "df")[[3]]), "fdr")
rLmFDR3 = p.adjust(pt2(rmincLm[, 11], attr(rmincLm, "df")[[4]]), "fdr")
rLmFDR4 = p.adjust(pt2(rmincLm[, 12], attr(rmincLm, "df")[[5]]), "fdr")
rLmFDR5 = p.adjust(pt2(rmincLm[, 13], attr(rmincLm, "df")[[6]]), "fdr")
rLmFDR6 = p.adjust(pt2(rmincLm[, 14], attr(rmincLm, "df")[[7]]), "fdr")


rmincFDR = verboseRun("anatFDR(rmincLm)", getOption("verbose"))

test_that("anatFDR Three Factors Interaction", {
  expect_equal(unname(rLmFDR1[1]), unname(rmincFDR[1, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[2]), unname(rmincFDR[2, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR1[3]), unname(rmincFDR[3, 2]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[1]), unname(rmincFDR[1, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[2]), unname(rmincFDR[2, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR2[3]), unname(rmincFDR[3, 3]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[1]), unname(rmincFDR[1, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[2]), unname(rmincFDR[2, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR3[3]), unname(rmincFDR[3, 4]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[1]), unname(rmincFDR[1, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[2]), unname(rmincFDR[2, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR4[3]), unname(rmincFDR[3, 5]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR5[1]), unname(rmincFDR[1, 6]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR5[2]), unname(rmincFDR[2, 6]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR5[3]), unname(rmincFDR[3, 6]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR6[1]), unname(rmincFDR[1, 7]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR6[2]), unname(rmincFDR[2, 7]), ignore_attr = TRUE)
  expect_equal(unname(rLmFDR6[3]), unname(rmincFDR[3, 7]), ignore_attr = TRUE)
})

test_that("rownames are preserved", {
  expect_equal(rownames(rmincFDR), rownames(rmincLm))
})
