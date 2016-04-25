library(testthat)
library(lme4)

context("mincLmer - basic test")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "test_data_set.csv"))
maskfile <- file.path(dataPath, "testminc-mask.mnc")

# pick a voxel inside the mask
voxelIndex <- 453 # for later comparisons
gf$v <- mincGetVoxel(gf$jacobians_fixed_2, 4, 5, 2)
gf$coil <- as.factor(gf$coil)
l <- lmer(v ~ Sex + (1|coil), gf)

vsreml <- verboseRun("mincLmer(jacobians_fixed_2 ~ Sex + (1|coil), gf, mask=maskfile)",getOption("verbose"))

test_that("mincLmer basics", {
  expect_that(vsreml[voxelIndex,1], is_equivalent_to(fixef(l)[1]))
  expect_that(vsreml[voxelIndex,2], is_equivalent_to(fixef(l)[2]))
})
            
context("mincLmer - maximum likelihood test")
l <- lmer(v ~ Sex + (1|coil), gf, REML=F)

vsml <- verboseRun("mincLmer(jacobians_fixed_2 ~ Sex + (1|coil), gf, REML=F, mask=maskfile)",getOption("verbose"))

test_that("mincLmer basics", {
  expect_that(vsml[voxelIndex,1], is_equivalent_to(fixef(l)[1]))
  expect_that(vsml[voxelIndex,2], is_equivalent_to(fixef(l)[2]))
})

context("mincLmer - log likelihood ratios")

vsml2 <-  verboseRun("mincLmer(jacobians_fixed_2 ~ 1 + (1|coil), gf, REML=F, mask=maskfile)",getOption("verbose"))

l2 <- lmer(v ~ 1 + (1|coil), gf, REML=F)

test_that("logLikRatio", {
  expect_that( mincLogLikRatio(vsreml, vsml), throws_error() )
  expect_that( mincLogLikRatio(vsml, vsml2)[voxelIndex],
              is_equivalent_to(anova(l,l2)[2,6]) )
})

context("mincLmer - estimate DF")
test_that("empty DF by default", {
  expect_that( attr(vsreml, "df"), is_equivalent_to(NULL))
  expect_that( mincFDR(vsreml), throws_error() )
})

vsreml <- mincLmerEstimateDF(vsreml)
df <- attr(vsreml, "df")
test_that("DF within reasonable range", {
  expect_that( df[[2]], is_less_than(nrow(gf)+1))
  expect_that( df[[2]], is_more_than(1))
})
          

#TODO: test parallel variants of mincLmer

