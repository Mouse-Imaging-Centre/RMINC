library(testthat)
#testthat test script for functions that call mincSummary
# mincMean, mincSd, mincVar, mincSum
context("mincSummary (Mean, Sd, Var, Sum,t-test,correlation,wilcoxon)")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"), stringsAsFactors = TRUE)
gf$vox <- mincGetVoxel(gf$jacobians_0.2, 0, 0, 0)

#Calculate mean, sd, variance, and sum
mm <- verboseRun("mincMean(gf$jacobians_0.2)", getOption("verbose"))
ms <- verboseRun("mincSd(gf$jacobians_0.2)", getOption("verbose"))
mv <- verboseRun("mincVar(gf$jacobians_0.2)", getOption("verbose"))
ms2 <- verboseRun("mincSum(gf$jacobians_0.2)", getOption("verbose"))

#Calculate mean, sd, variance and sum with a factor...
mms <- verboseRun("mincMean(gf$jacobians_0.2, gf$Strain)", getOption("verbose"))
mss <- verboseRun("mincSd(gf$jacobians_0.2, gf$Strain)", getOption("verbose"))
mvs <- verboseRun("mincVar(gf$jacobians_0.2, gf$Strain)", getOption("verbose"))
ms2s <- verboseRun("mincSum(gf$jacobians_0.2, gf$Strain)", getOption("verbose"))

#...and verify with tapply commands
mt <- tapply(gf$vox, gf$Strain, mean)
st <- tapply(gf$vox, gf$Strain, sd)
vt <- tapply(gf$vox, gf$Strain, var)
s2t <- tapply(gf$vox, gf$Strain, sum)

test_that("mincSummary functions", {
  expect_equal(mean(gf$vox), mm[[1, 1]])
  expect_equal(sd(gf$vox), ms[[1, 1]])
  expect_equal(var(gf$vox), mv[[1, 1]])
  expect_equal(sum(gf$vox), ms2[[1, 1]])
  expect_equal(mt[1], mms[1, 1])
  expect_equal(mt[2], mms[1, 2])
  expect_equal(st[1], mss[1, 1])
  expect_equal(st[2], mss[1, 2])
  expect_equal(vt[1], mvs[1, 1])
  expect_equal(vt[2], mvs[1, 2])
  expect_equal(s2t[1], ms2s[1, 1])
  expect_equal(s2t[2], ms2s[1, 2])
})


test_that("summary functions count groups from the largest group code (#241)", {
  # Group codes 2, 1, 1, ..., 0, 0. The old C code counted only 2 groups
  # for this order and wrote past the end of its buffers.
  n <- nrow(gf)
  n_b <- floor((n - 1) / 2)
  g3 <- factor(c("c", rep("b", n_b), rep("a", n - 1 - n_b)))

  mm3 <- mincMean(gf$jacobians_0.2, g3)
  ms3 <- mincSum(gf$jacobians_0.2, g3)
  mv3 <- mincVar(gf$jacobians_0.2, g3)

  expect_equal(dim(mm3)[2], 3)
  expect_equal(colnames(mm3), c("a", "b", "c"))
  expect_equal(unname(mm3[1, ]), unname(c(tapply(gf$vox, g3, mean))))
  expect_equal(unname(ms3[1, ]), unname(c(tapply(gf$vox, g3, sum))))
  # group "c" has one file, so only check the variance of "a" and "b"
  expect_equal(unname(mv3[1, 1:2]), unname(c(tapply(gf$vox, g3, var))[1:2]))
})

test_that("mincSummary always returns a one-row-per-voxel matrix (#309)", {
  n_vox <- prod(minc.dimensions.sizes(gf$jacobians_0.2[1]))

  m_mean <- mincSummary(gf$jacobians_0.2, method = "mean")
  expect_equal(dim(m_mean), c(n_vox, 1))
  expect_s3_class(m_mean, "mincMultiDim")

  m_t <- mincSummary(gf$jacobians_0.2, gf$Strain, method = "t-test")
  expect_equal(dim(m_t), c(n_vox, 1))
  expect_s3_class(m_t, "mincMultiDim")
  expect_equal(attr(m_t, "likeVolume"), as.character(gf$jacobians_0.2[1]))

  # median is not computed by the C code; before, it returned values that
  # were never set
  expect_error(
    mincSummary(gf$jacobians_0.2, method = "median"),
    "mincApply"
  )
})

test_that("mincSummary checks the grouping", {
  n <- nrow(gf)
  # character grouping is converted to a factor
  chr <- rep(c("x", "y"), length.out = n)
  mchr <- mincSummary(gf$jacobians_0.2, chr, method = "mean")
  expect_equal(colnames(mchr), c("x", "y"))
  expect_equal(unname(mchr[1, ]), unname(c(tapply(gf$vox, chr, mean))))

  # unused levels are dropped, so colnames match the columns
  unused <- factor(chr, levels = c("x", "y", "z"))
  expect_equal(colnames(mincSummary(gf$jacobians_0.2, unused)), c("x", "y"))

  expect_error(
    mincSummary(gf$jacobians_0.2, replace(chr, 1, NA)),
    "must not contain NA"
  )
  expect_error(
    mincSummary(gf$jacobians_0.2, chr[-1]),
    "one value for each file"
  )
})


mtt <- verboseRun("mincTtest(gf$jacobians_0.2,gf$Strain)", getOption("verbose"))
ttt <- t.test(vox ~ Strain, data = gf)

test_that("ttest", {
  expect_equivalent(ttt$statistic, mtt[1])
})


gf_paired <- gf[1:20, ]
mptt <- verboseRun(
  "mincPairedTtest(gf_paired$jacobians_0.2,gf_paired$Strain)",
  getOption("verbose")
) # To Do: Ask case where unequal lengths
pttt <- t.test(
  gf_paired$vox[gf_paired$Strain == levels(gf_paired$Strain)[1]],
  gf_paired$vox[gf_paired$Strain == levels(gf_paired$Strain)[2]],
  paired = TRUE
)

test_that("paired ttest", {
  expect_equivalent(pttt$statistic, mptt[1])
})


mc <- verboseRun(
  "mincCorrelation(gf$jacobians_0.2,gf$Weight)",
  getOption("verbose")
)
tc <- cor(gf$Weight, gf$vox)

test_that("correlation", {
  expect_equivalent(tc, mc[1])
})


gf_paired$vox_round <- round(gf_paired$vox, 8)
suppressWarnings(
  tw <- wilcox.test(vox_round ~ Strain, data = gf_paired)
)
mw <- verboseRun(
  "mincWilcoxon(gf_paired$jacobians_0.2,gf_paired $Strain)",
  getOption("verbose")
)
test_that("wilcoxon-ties", {
  expect_equivalent(tw[[1]], mw[1])
})

gf$vox <- mincGetVoxel(gf$jacobians_0.2, 5, 5, 5)
gf_paired <- gf[1:20, ]
mw <- verboseRun(
  "mincWilcoxon(gf_paired$jacobians_0.2,gf_paired $Strain)",
  getOption("verbose")
)
tw <- wilcox.test(vox ~ Strain, data = gf_paired)
test_that("wilcoxon", {
  expect_equivalent(tw[[1]], mw[15 * 15 * 5 + 15 * 5 + 6])
})

# mincFDR

# t-test
rttest = p.adjust(pt2(mtt[, 1], attr(mtt, "df")), "fdr")
verboseRun("rmincFDR = mincFDR(mtt)", getOption("verbose"))
test_that("mincTtest works with mincFDR", {
  expect_equal(as.numeric(rmincFDR), rttest, tolerance = 0.0001)
})

# paired t-test
rttest = p.adjust(pt2(mptt[, 1], attr(mptt, "df")), "fdr")
verboseRun("rmincFDR = mincFDR(mptt)", getOption("verbose"))
test_that("mincPairedTtest works with mincFDR", {
  expect_equal(as.numeric(rmincFDR), rttest, tolerance = 0.0001)
})

# wilcox
# uses pwilcox to compute p-values
rWilcoxFDR = p.adjust(
  1 - pwilcox(mw[, 1], attr(mw, "m"), attr(mw, "n"), lower.tail = FALSE),
  "fdr"
)
verboseRun("rmincWilcoxFDR = mincFDR(mw)", getOption("verbose"))
test_that("mincWilcoxon works with mincFDR", {
  expect_equal(as.numeric(rmincWilcoxFDR), rWilcoxFDR, tolerance = 0.0001)
})

ithreshold = which(mw[, 1] == attr(rmincWilcoxFDR, "thresholds")[[1]])
test_that("mincWilcoxon is thresholded properly by mincFDR", {
  expect_lt(rWilcoxFDR[ithreshold[1]], 0.1)
})
