library(testthat)
context("anatModel subsetting and conversion")

set.seed(1)
make_group <- function(level, m1, m2) {
  data.frame(
    vol = rnorm(30, m1, 800),
    vol2 = rnorm(30, m2, 400),
    IV = factor(level, levels = c("Control", "Treatment"))
  )
}
groups <- rbind(
  make_group("Control", 10000, 5000),
  make_group("Treatment", 10600, 4800)
)
anat2 <- as.matrix(groups[c("vol", "vol2")])
anat1 <- as.matrix(groups["vol"])

verboseRun(fit2 <- anatLm(~IV, data = groups, anat = anat2))
verboseRun(fit1 <- anatLm(~IV, data = groups, anat = anat1))

test_that("one row of an anatModel is still an anatModel (#187)", {
  row1 <- fit2[1, ]
  expect_s3_class(row1, "anatModel")
  expect_equal(dim(row1), c(1, ncol(fit2)))
  expect_equal(rownames(row1), "vol")
  expect_equal(colnames(row1), colnames(fit2))
  expect_equal(attr(row1, "stat-type"), attr(fit2, "stat-type"))
  expect_equal(attr(row1, "df"), attr(fit2, "df"))

  # selecting one column still gives a vector, as before
  expect_null(dim(fit2[, "tvalue-IVTreatment"]))
  expect_equal(names(fit2[, "tvalue-IVTreatment"]), c("vol", "vol2"))
  expect_null(dim(fit2[1, "tvalue-IVTreatment"]))
  # unless drop = FALSE
  expect_equal(dim(fit2[, "tvalue-IVTreatment", drop = FALSE]), c(2, 1))
})

test_that("anatLm works with an anatomy matrix of one structure (#282)", {
  expect_s3_class(fit1, "anatModel")
  expect_equal(dim(fit1), c(1, ncol(fit2)))
  expect_equal(unclass(fit1)[1, ], unclass(fit2)[1, ], check.attributes = FALSE)

  out <- capture.output(print(fit1))
  expect_true(any(grepl("tvalue-IVTreatment", out)))
  expect_true(any(grepl("^vol ", out)))

  es1 <- anatEffectSize(fit1)
  es2 <- anatEffectSize(fit2)
  expect_equal(dim(es1), c(1, ncol(es2)))
  expect_equal(unclass(es1)[1, ], unclass(es2)[1, ], check.attributes = FALSE)

  verboseRun(q1 <- anatFDR(fit1))
  expect_equal(nrow(q1), 1)
})

test_that("as.data.frame can keep the structure names in a column (#125)", {
  df <- as.data.frame(fit2)
  expect_equal(rownames(df), c("vol", "vol2"))
  expect_equal(names(df), colnames(fit2))

  dfs <- as.data.frame(fit2, structure_column = "Structure")
  expect_equal(names(dfs), c("Structure", colnames(fit2)))
  expect_equal(dfs$Structure, c("vol", "vol2"))

  # the names stay after a filter that drops row names
  kept <- dfs[dfs$`tvalue-IVTreatment` > 0, , drop = FALSE]
  rownames(kept) <- NULL
  expect_equal(kept$Structure, "vol")

  expect_equal(
    as.data.frame(fit1, structure_column = "Structure")$Structure,
    "vol"
  )
  expect_error(as.data.frame(fit2, structure_column = 1), "one string")
})
