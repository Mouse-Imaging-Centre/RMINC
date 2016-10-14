library(testthat)
context("mincLm - two group test")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"))

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
# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm = verboseRun("mincLm(voxel_left_file ~ Sex, gftest)",getOption("verbose"))


test_that("mincLm Two Factors",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})

context("mincLm - two group test with interaction")

# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm = verboseRun("mincLm(voxel_left_file~Sex*Scale,gftest)",getOption("verbose"))
gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left~Sex*Scale,gftest))

test_that("mincLm interaction",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[3,1]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[4,1]))
	expect_that(rmincLm[1,7],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,8],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(rmincLm[1,9],is_equivalent_to(rLm$coefficients[3,3]))	
	expect_that(rmincLm[1,10],is_equivalent_to(rLm$coefficients[4,3]))	
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})

context("mincLm - three group test")

# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm = verboseRun("mincLm(voxel_left_file~coil,gftest)",getOption("verbose"))

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left~Coil,gftest))

test_that("mincLm Three Factors",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[3,1]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,7],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(rmincLm[1,8],is_equivalent_to(rLm$coefficients[3,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})

context("mincLm - three group test with interaction")

# silence the output of mincLm, in order to make the test output information is more clear to read
rmincLm = verboseRun("mincLm(voxel_left_file~Scale*Coil,gftest)",getOption("verbose"))

gftest$voxel_left = voxel_left
rLm = summary(lm(voxel_left~Scale*Coil,gftest))

test_that("mincLm Three Factors Interaction",{
	expect_that(rmincLm[1,1],is_equivalent_to(rLm$fstatistic[1]))
	expect_that(rmincLm[1,2],is_equivalent_to(rLm$r.squared[1]))
	expect_that(rmincLm[1,3],is_equivalent_to(rLm$coefficients[1,1]))
	expect_that(rmincLm[1,4],is_equivalent_to(rLm$coefficients[2,1]))
	expect_that(rmincLm[1,5],is_equivalent_to(rLm$coefficients[3,1]))
	expect_that(rmincLm[1,6],is_equivalent_to(rLm$coefficients[4,1]))
	expect_that(rmincLm[1,7],is_equivalent_to(rLm$coefficients[5,1]))
	expect_that(rmincLm[1,8],is_equivalent_to(rLm$coefficients[6,1]))
	expect_that(rmincLm[1,9],is_equivalent_to(rLm$coefficients[1,3]))
	expect_that(rmincLm[1,10],is_equivalent_to(rLm$coefficients[2,3]))
	expect_that(rmincLm[1,11],is_equivalent_to(rLm$coefficients[3,3]))
	expect_that(rmincLm[1,12],is_equivalent_to(rLm$coefficients[4,3]))
	expect_that(rmincLm[1,13],is_equivalent_to(rLm$coefficients[5,3]))
	expect_that(rmincLm[1,14],is_equivalent_to(rLm$coefficients[6,3]))
	expect_that(attr(rmincLm,"df")[[2]],is_equivalent_to(rLm$df[2]))
})

test_that("mincLm local multicore works", {
  skip_on_cran()
  skip_on_travis()
  
  if(Sys.getenv("TEST_Q_MINC") != "yes") 
    skip("qMinc tests disabled")
  
  verboseRun(prlm <- mincLm(voxel_left_file~Scale*Coil,gftest, parallel = c("local", 3)))
  
  expect_equal(rmincLm, prlm, check.attributes = FALSE)
})

# test_that("mincLm queue parallel works", {
#   verboseRun(qrlm <- mincLm(voxel_left_file~Scale*Coil,gftest, parallel = c("sge", 2)))
#   
#   expect_equal(rmincLm, qrlm, check.attributes = FALSE)
# })

