library(testthat)
context("vertexApply")

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

mm <- verboseRun("vertexMean(gftest$testFilesLeft)", getOption("verbose"))
ma <- verboseRun("vertexApply(gftest$testFilesLeft,mean)", getOption("verbose"))

test_that("vertexApply one output", {
  expect_equal(ma, mm, check.attributes = FALSE)
})

# Need to define global variable when running tests, but normally do not...
testFunc <- function(x) {
  return(c(1, 2))
}
ma <- verboseRun(
  "vertexApply(gftest$testFilesLeft, testFunc)",
  getOption("verbose")
)

test_that("vertexApply two output", {
  expect_equal(ma[, 1], rep(1, nrow(ma)))
  expect_equal(ma[, 2], rep(2, nrow(ma)))
})

test_that("transposed vertexApply masks vertices, not subjects (#217)", {
  vt <- vertexTable(gftest$testFilesLeft)
  n_vert <- nrow(vt)
  msk <- as.numeric(seq_len(n_vert) <= n_vert / 2)
  keep <- msk > .5

  expected <- colMeans(vt[keep, , drop = FALSE])

  verboseRun(
    tm <- vertexApply(gftest$testFilesLeft, mean, transpose = TRUE, mask = msk)
  )
  expect_equal(as.numeric(tm), unname(expected))

  # the mask can also be a file name or a one-column data frame
  mask_file <- tempfile(fileext = ".txt")
  on.exit(unlink(mask_file))
  writeLines(as.character(msk), mask_file)
  verboseRun(
    tm_file <- vertexApply(
      gftest$testFilesLeft, mean, transpose = TRUE, mask = mask_file
    )
  )
  expect_equal(as.numeric(tm_file), unname(expected))
  expect_equal(
    as.numeric(vertexApply(
      gftest$testFilesLeft, mean, transpose = TRUE, mask = data.frame(msk)
    )),
    unname(expected)
  )

  expect_error(
    vertexApply(gftest$testFilesLeft, mean, transpose = TRUE, mask = msk[-1]),
    "vertices"
  )
})

test_that("vertexAtlasApply applies the mask to vertices and atlas (#217)", {
  vt <- vertexTable(gftest$testFilesLeft)
  n_vert <- nrow(vt)
  atlas <- rep_len(1:3, n_vert)
  msk <- as.numeric(seq_len(n_vert) %% 4 != 0)
  keep <- msk > .5

  expected <- apply(
    vt[keep, , drop = FALSE], 2,
    function(x) tapply(x, atlas[keep], mean)
  )

  verboseRun(
    vam <- vertexAtlasApply(gftest$testFilesLeft, atlas, mean, mask = msk)
  )
  expect_equal(dim(vam), dim(expected))
  expect_equal(as.numeric(vam), as.numeric(expected))
})
