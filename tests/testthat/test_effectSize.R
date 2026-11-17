library(testthat)
context("EffectSize")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

# vertex data
gfv <- read.csv(file.path(dataPath, "subject.csv"), stringsAsFactors = TRUE)
gfv$testFilesLeft <- file.path(
  dataPath,
  paste0("vertex", c(2, 3, 4, 3, 1, 2, 4, 2, 3, 1), ".txt")
)

# voxel data
gf <- read.csv(file.path(dataPath, "test_data_set.csv"), stringsAsFactors = TRUE)
gftest <- gf[1:10, ]

test_that("vertexEffectSize works when the formula is a local variable (#305)", {
  verboseRun(direct <- vertexEffectSize(vertexLm(testFilesLeft ~ Sex, gfv)))

  fit_in_function <- function(formula_string, d) {
    form <- as.formula(formula_string)
    model <- vertexLm(form, data = d)
    vertexEffectSize(model)
  }
  verboseRun(local_es <- fit_in_function("testFilesLeft ~ Sex", gfv))

  expect_equal(unclass(local_es), unclass(direct), check.attributes = FALSE)
  expect_true(any(grepl("^hedgesg-", colnames(local_es))))
})

test_that("mincEffectSize works on a model that is not a global variable (#305)", {
  verboseRun(direct <- mincEffectSize(mincLm(jacobians_fixed_2 ~ Sex, gftest)))

  fit_in_function <- function(d) {
    local_model <- mincLm(jacobians_fixed_2 ~ Sex, d)
    mincEffectSize(local_model)
  }
  verboseRun(local_es <- fit_in_function(gftest))

  expect_s3_class(local_es, "mincEffectSize")
  expect_equal(unclass(local_es), unclass(direct), check.attributes = FALSE)
})
