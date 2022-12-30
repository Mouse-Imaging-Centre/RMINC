library(testthat)
context("mincAnova")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "RMINC-test-data-main/rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"), stringsAsFactors = TRUE)

voxel_left <- mincGetVoxel(gf$jacobians_fixed_2[1:10], 0,0,0)
voxel_right <- mincGetVoxel(gf$jacobians_fixed_2[11:20], 0,0,0)
Sex <- gf$Sex[1:10]
Scale <- gf$scale[1:10]
Coil <- as.factor(gf$coil[1:10])
gftest <- gf[1:10,]
gftest$voxel_right <- (gf$jacobians_fixed_2[11:20])
gftest$voxel_left_file <- gf$jacobians_fixed_2[1:10]
gftest$Coil <- Coil

rAnova = anova(lm(voxel_left  ~ Sex))
rmincAnova = verboseRun("mincAnova(voxel_left_file ~ Sex, gftest)",getOption("verbose"))

test_that("mincAnova Two Factors",{
	expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
})

rmincAnova = verboseRun("mincAnova(voxel_left_file~Sex*Scale,gftest)",getOption("verbose"))


gftest$voxel_left <- voxel_left
rAnova = anova(lm(voxel_left~Sex*Scale,gftest))

test_that("mincAnova interaction",{
	expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
	expect_that(rmincAnova[1,2],is_equivalent_to(rAnova$F[2]))
	expect_that(rmincAnova[1,3],is_equivalent_to(rAnova$F[3]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[4]))
	expect_that(attr(rmincAnova,"df")[[2]][1],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[2]][2],is_equivalent_to(rAnova$Df[4]))
	expect_that(attr(rmincAnova,"df")[[3]][1],is_equivalent_to(rAnova$Df[3]))
	expect_that(attr(rmincAnova,"df")[[3]][2],is_equivalent_to(rAnova$Df[4]))
})


rmincAnova = verboseRun("mincAnova(voxel_left_file~Coil,gftest)",getOption("verbose"))


gftest$voxel_left = voxel_left
rAnova = anova(lm(voxel_left~Coil,gftest))

test_that("mincAnova Three Factors",{
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
})


gftest$voxel_left = voxel_left
rAnova = anova(lm(voxel_left~Scale*Coil,gftest))

rmincAnova = verboseRun("mincAnova(voxel_left_file~Scale*Coil,gftest)",getOption("verbose"))

test_that("mincAnova Three Factors Interaction",{
         expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
         expect_that(rmincAnova[1,2],is_equivalent_to(rAnova$F[2]))
         expect_that(rmincAnova[1,3],is_equivalent_to(rAnova$F[3]))
         expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
         expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[4]))
         expect_that(attr(rmincAnova,"df")[[2]][1],is_equivalent_to(rAnova$Df[2]))
         expect_that(attr(rmincAnova,"df")[[2]][2],is_equivalent_to(rAnova$Df[4]))
         expect_that(attr(rmincAnova,"df")[[3]][1],is_equivalent_to(rAnova$Df[3]))
         expect_that(attr(rmincAnova,"df")[[3]][2],is_equivalent_to(rAnova$Df[4])) 
})

test_that("mincAnova local multicore works", {
  skip_on_cran()
  skip_on_travis()
  
  if(Sys.getenv("TEST_Q_MINC") != "yes") 
    skip("qMinc tests disabled")
  
  verboseRun(pra <- mincAnova(voxel_left_file~Scale*Coil,gftest, parallel = c("local", 3)))
  
  expect_equal(rmincAnova, pra, check.attributes = FALSE)
})
