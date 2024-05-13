library(testthat)
#testthat test script for functions that call mincSummary
# mincMean, mincSd, mincVar, mincSum 
context("mincSummary (Mean, Sd, Var, Sum,t-test,correlation,wilcoxon)")

if(!exists("dataPath")) {
  dataPath <- tempdir()
}
  
getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"))
gf$vox <- mincGetVoxel(gf$jacobians_0.2, 0, 0, 0)

#Calculate mean, sd, variance, and sum 
mm <- verboseRun("mincMean(gf$jacobians_0.2)",getOption("verbose"))
ms <- verboseRun("mincSd(gf$jacobians_0.2)",getOption("verbose"))
mv <- verboseRun("mincVar(gf$jacobians_0.2)",getOption("verbose"))
ms2 <- verboseRun("mincSum(gf$jacobians_0.2)",getOption("verbose"))

#Calculate mean, sd, variance and sum with a factor...
mms <- verboseRun("mincMean(gf$jacobians_0.2, gf$Strain)",getOption("verbose"))
mss <- verboseRun("mincSd(gf$jacobians_0.2, gf$Strain)",getOption("verbose"))
mvs <- verboseRun("mincVar(gf$jacobians_0.2, gf$Strain)",getOption("verbose"))
ms2s <- verboseRun("mincSum(gf$jacobians_0.2, gf$Strain)",getOption("verbose"))

#...and verify with tapply commands
mt <- tapply(gf$vox, gf$Strain, mean)
st <- tapply(gf$vox, gf$Strain, sd)
vt <- tapply(gf$vox, gf$Strain, var)
s2t <- tapply(gf$vox, gf$Strain, sum)

test_that("mincSummary functions", {
    expect_equal(mean(gf$vox), mm[[1,1]])
    expect_equal(sd(gf$vox), ms[[1,1]])
    expect_equal(var(gf$vox), mv[[1,1]])
    expect_equal(sum(gf$vox), ms2[[1,1]])
    expect_equal(mt[1], mms[1,1])
    expect_equal(mt[2], mms[1,2])
    expect_equal(st[1], mss[1,1])
    expect_equal(st[2], mss[1,2])
    expect_equal(vt[1], mvs[1,1])
    expect_equal(vt[2], mvs[1,2])
    expect_equal(s2t[1], ms2s[1,1])
    expect_equal(s2t[2], ms2s[1,2])
    
})


mtt <- verboseRun("mincTtest(gf$jacobians_0.2,gf$Strain)",getOption("verbose"))
ttt <- t.test(vox~Strain,data=gf)

test_that("ttest", {
    expect_equivalent(ttt$statistic, mtt[1])
})


gf_paired <- gf[1:20,];
mptt <- verboseRun("mincPairedTtest(gf_paired$jacobians_0.2,gf_paired$Strain)",getOption("verbose")) # To Do: Ask case where unequal lengths
pttt <- t.test(vox~Strain,data=gf_paired ,paired=TRUE)

test_that("paired ttest", {
    expect_equivalent(pttt$statistic, mptt[1])
})


mc <- verboseRun("mincCorrelation(gf$jacobians_0.2,gf$Weight)",getOption("verbose"))
tc <- cor(gf$Weight,gf$vox)

test_that("correlation", {
    expect_equivalent(tc, mc[1])
})


gf_paired$vox_round <- round(gf_paired$vox,8)
suppressWarnings(
  tw <- wilcox.test(vox_round~Strain,data=gf_paired)
)
mw <- verboseRun("mincWilcoxon(gf_paired$jacobians_0.2,gf_paired $Strain)",getOption("verbose"))
test_that("wilcoxon-ties", {
    expect_equivalent(tw[[1]], mw[1])
})

gf$vox <- mincGetVoxel(gf$jacobians_0.2, 5, 5, 5)
gf_paired <- gf[1:20,];
mw <- verboseRun("mincWilcoxon(gf_paired$jacobians_0.2,gf_paired $Strain)",getOption("verbose"))
tw <- wilcox.test(vox~Strain,data=gf_paired)
test_that("wilcoxon", {
    expect_equivalent(tw[[1]], mw[15*15*5 + 15*5 + 6])
})

# mincFDR

# t-test
rttest = p.adjust( pt2(mtt[,1],attr(mtt,"df")),"fdr")
verboseRun("rmincFDR = mincFDR(mtt)", getOption("verbose"))
test_that("mincTtest works with mincFDR",{
  expect_equal(as.numeric(rmincFDR), rttest,tolerance = 0.0001) 
})

# paired t-test
rttest = p.adjust( pt2(mptt[,1],attr(mptt,"df")),"fdr")
verboseRun("rmincFDR = mincFDR(mptt)", getOption("verbose"))
test_that("mincPairedTtest works with mincFDR",{
  expect_equal(as.numeric(rmincFDR), rttest,tolerance = 0.0001) 
})

# wilcox 
# uses pwilcox to compute p-values 
rWilcoxFDR = p.adjust(1 - pwilcox(mw[,1],attr(mw,"m"),attr(mw,"n"),lower.tail = FALSE),"fdr")
verboseRun("rmincWilcoxFDR = mincFDR(mw)", getOption("verbose"))
test_that("mincWilcoxon works with mincFDR",{
  expect_equal(as.numeric(rmincWilcoxFDR), rWilcoxFDR,tolerance = 0.0001) 
})

ithreshold = which(mw[,1] == attr(rmincWilcoxFDR,"thresholds")[[1]])
test_that("mincWilcoxon is thresholded properly by mincFDR",{
	expect_lt(rWilcoxFDR[ithreshold[1]], 0.1) 
})

