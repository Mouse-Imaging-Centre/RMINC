library(testthat)
context("anatApply")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"))
first_file <- gf$jacobians_0.2[1]
segmentation <- file.path(dataPath, "test_segmentation.mnc")

### Stub