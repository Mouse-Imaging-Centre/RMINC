#testthat test script for functions that call mincSummary
# mincMean, mincSd, mincVar, mincSum 
context("mincSummary (Mean, Sd, Var, Sum)")

gf <- read.csv("/tmp/rminctestdata/minc_summary_test_data.csv") # Note: this will change. Need standard data set. 
gf$vox <- mincGetVoxel(gf$jacobians_0.2, 0, 0, 0) # I think we should make sure this is a voxel where we are sure something is happening

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



