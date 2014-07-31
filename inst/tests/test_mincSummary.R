#testthat test script for functions that call mincSummary
# mincMean, mincSd, mincVar, mincSum 
context("mincSummary (Mean, Sd, Var, Sum,t-test,correlation,wilcoxon)")

gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") 
gf$vox <- mincGetVoxel(gf$jacobians_0.2, 0, 0, 0)

#Calculate mean, sd, variance, and sum 
sink("/dev/null"); mm <- mincMean(gf$jacobians_0.2); sink();
sink("/dev/null"); ms <- mincSd(gf$jacobians_0.2); sink();
sink("/dev/null"); mv <- mincVar(gf$jacobians_0.2); sink();
sink("/dev/null"); ms2 <- mincSum(gf$jacobians_0.2); sink();

#Calculate mean, sd, variance and sum with a factor...
sink("/dev/null"); mms <- mincMean(gf$jacobians_0.2, gf$Strain); sink();
sink("/dev/null"); mss <- mincSd(gf$jacobians_0.2, gf$Strain); sink();
sink("/dev/null"); mvs <- mincVar(gf$jacobians_0.2, gf$Strain); sink();
sink("/dev/null"); ms2s <- mincSum(gf$jacobians_0.2, gf$Strain); sink();

#...and verify with tapply commands
mt <- tapply(gf$vox, gf$Strain, mean)
st <- tapply(gf$vox, gf$Strain, sd)
vt <- tapply(gf$vox, gf$Strain, var)
s2t <- tapply(gf$vox, gf$Strain, sum)

test_that("mincSummary functions", {
    expect_equal(mean(gf$vox), mm[1,1])
    expect_equal(sd(gf$vox), ms[1,1])
    expect_equal(var(gf$vox), mv[1,1])
    expect_equal(sum(gf$vox), ms2[1,1])
    expect_equal(mt[1], mms[1,1])
    expect_equal(mt[2], mms[1,2])
    expect_equal(st[1], mss[1,1])
    expect_equal(st[2], mss[1,2])
    expect_equal(vt[1], mvs[1,1])
    expect_equal(vt[2], mvs[1,2])
    expect_equal(s2t[1], ms2s[1,1])
    expect_equal(s2t[2], ms2s[1,2])
    
})


sink("/dev/null"); mtt <- mincTtest(gf$jacobians_0.2,gf$Strain); sink();
ttt <- t.test(vox~Strain,data=gf)

test_that("ttest", {
    expect_equivalent(ttt$statistic, mtt[1])
})


gf_paired = gf[1:20,];
sink("/dev/null"); mtt <- mincPairedTtest(gf_paired $jacobians_0.2,gf_paired $Strain); sink(); # To Do: Ask case where unequal lengths
ttt <- t.test(vox~Strain,data=gf_paired ,paired=TRUE)

test_that("paired ttest", {
    expect_equivalent(ttt$statistic, mtt[1])
})


sink("/dev/null"); mc <- mincCorrelation(gf$jacobians_0.2,gf$Weight); sink();
tc <- cor(gf$Weight,gf$vox)

test_that("correlation", {
    expect_equivalent(tc, mc[1])
})

gf$vox <- mincGetVoxel(gf$jacobians_0.2, 5, 5, 5)
gf_paired = gf[1:20,];
sink("/dev/null"); mw <- mincWilcoxon(gf_paired$jacobians_0.2,gf_paired $Strain); sink();
tw <- wilcox.test(vox~Strain,data=gf_paired)
gf$vox <- mincGetVoxel(gf$jacobians_0.2, 5, 5, 5)
gf_paired = gf[1:20,];
test_that("wilcoxon", {
    expect_equivalent(100-tw[[1]], mw[15*15*5 + 15*5 + 6])
})

