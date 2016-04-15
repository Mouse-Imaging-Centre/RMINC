library(testthat)
context("mincFDR")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
voxel_left <- mincGetVoxel(gf$jacobians_fixed_2[1:10], 0,0,0)
voxel_right <- mincGetVoxel(gf$jacobians_fixed_2[11:20], 0,0,0)
Sex <- gf$Sex[1:10]
Scale <- gf$scale[1:10]
Coil <- as.factor(gf$coil[1:10])
gf$coil <- as.factor(gf$coil)
gftest <- gf[1:10,]
gftest$voxel_right <- (gf$jacobians_fixed_2[11:20])
gftest$voxel_left_file <- gf$jacobians_fixed_2[1:10]

rLm = summary(lm(voxel_left  ~ Sex))


rmincLm <- verboseRun("mincLm(voxel_left_file ~ Sex, gftest)",getOption("verbose"))


rLmFDR1 = p.adjust( pt2(rmincLm[,5],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,6],attr(rmincLm,"df")[[3]]),"fdr")

rmincFDR = verboseRun("mincFDR(rmincLm)",getOption("verbose"))

test_that("mincFDR Two Factors",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
})

rmincLm <- verboseRun("mincLm(voxel_left_file~Sex*Scale,gftest)",getOption("verbose"))

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left~Sex*Scale,gftest))

rLmFDR1 = p.adjust( pt2(rmincLm[,7],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,8],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,9],attr(rmincLm,"df")[[4]]),"fdr")
rLmFDR4 = p.adjust( pt2(rmincLm[,10],attr(rmincLm,"df")[[5]]),"fdr")

rmincFDR = verboseRun("mincFDR(rmincLm)",getOption("verbose"))


test_that("mincFDR interaction",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
	expect_that(rLmFDR3[1],is_equivalent_to(rmincFDR[1,4]))
	expect_that(rLmFDR3[2],is_equivalent_to(rmincFDR[2,4]))
	expect_that(rLmFDR3[3],is_equivalent_to(rmincFDR[3,4]))
	expect_that(rLmFDR4[1],is_equivalent_to(rmincFDR[1,5]))
	expect_that(rLmFDR4[2],is_equivalent_to(rmincFDR[2,5]))
	expect_that(rLmFDR4[3],is_equivalent_to(rmincFDR[3,5]))
})

rmincLm <- verboseRun("rmincLm = mincLm(voxel_left_file~coil,gftest)",getOption("verbose"))

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left~Coil,gftest))


rLmFDR1 = p.adjust( pt2(rmincLm[,6],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,7],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,8],attr(rmincLm,"df")[[4]]),"fdr")

rmincFDR = verboseRun("mincFDR(rmincLm)",getOption("verbose"))

test_that("mincFDR Three Factors",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
	expect_that(rLmFDR3[1],is_equivalent_to(rmincFDR[1,4]))
	expect_that(rLmFDR3[2],is_equivalent_to(rmincFDR[2,4]))
	expect_that(rLmFDR3[3],is_equivalent_to(rmincFDR[3,4]))
})


rmincLm <- verboseRun("mincLm(voxel_left_file~Scale*coil,gftest)",getOption("verbose"))

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left~Scale*Coil,gftest))


rLmFDR1 = p.adjust( pt2(rmincLm[,9],attr(rmincLm,"df")[[2]]),"fdr")
rLmFDR2 = p.adjust( pt2(rmincLm[,10],attr(rmincLm,"df")[[3]]),"fdr")
rLmFDR3 = p.adjust( pt2(rmincLm[,11],attr(rmincLm,"df")[[4]]),"fdr")
rLmFDR4 = p.adjust( pt2(rmincLm[,12],attr(rmincLm,"df")[[5]]),"fdr")
rLmFDR5 = p.adjust( pt2(rmincLm[,13],attr(rmincLm,"df")[[6]]),"fdr")
rLmFDR6 = p.adjust( pt2(rmincLm[,14],attr(rmincLm,"df")[[7]]),"fdr")


rmincFDR = verboseRun("mincFDR(rmincLm)",getOption("verbose"))

test_that("mincLm Three Factors Interaction",{
	expect_that(rLmFDR1[1],is_equivalent_to(rmincFDR[1,2]))
	expect_that(rLmFDR1[2],is_equivalent_to(rmincFDR[2,2]))
	expect_that(rLmFDR1[3],is_equivalent_to(rmincFDR[3,2]))
	expect_that(rLmFDR2[1],is_equivalent_to(rmincFDR[1,3]))
	expect_that(rLmFDR2[2],is_equivalent_to(rmincFDR[2,3]))
	expect_that(rLmFDR2[3],is_equivalent_to(rmincFDR[3,3]))
	expect_that(rLmFDR3[1],is_equivalent_to(rmincFDR[1,4]))
	expect_that(rLmFDR3[2],is_equivalent_to(rmincFDR[2,4]))
	expect_that(rLmFDR3[3],is_equivalent_to(rmincFDR[3,4]))
	expect_that(rLmFDR4[1],is_equivalent_to(rmincFDR[1,5]))
	expect_that(rLmFDR4[2],is_equivalent_to(rmincFDR[2,5]))
	expect_that(rLmFDR4[3],is_equivalent_to(rmincFDR[3,5]))
	expect_that(rLmFDR5[1],is_equivalent_to(rmincFDR[1,6]))
	expect_that(rLmFDR5[2],is_equivalent_to(rmincFDR[2,6]))
	expect_that(rLmFDR5[3],is_equivalent_to(rmincFDR[3,6]))
	expect_that(rLmFDR6[1],is_equivalent_to(rmincFDR[1,7]))
	expect_that(rLmFDR6[2],is_equivalent_to(rmincFDR[2,7]))
	expect_that(rLmFDR6[3],is_equivalent_to(rmincFDR[3,7]))
})
