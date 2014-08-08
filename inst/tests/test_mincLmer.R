context("mincLmer - basic test")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
library(lme4)
gf$v <- mincGetVoxel(gf$jacobians_fixed_2, 0, 0, 0)
gf$coil <- as.factor(gf$coil)
l <- lmer(v ~ Sex + (1|coil), gf)
sink("/dev/null")
vsreml <- mincLmer(jacobians_fixed_2 ~ Sex + (1|coil), gf)
sink()
test_that("mincLmer basics", {
  expect_that(vsreml[1,1], is_equivalent_to(fixef(l)[1]))
  expect_that(vsreml[1,2], is_equivalent_to(fixef(l)[2]))
})
            
context("mincLmer - maximum likelihood test")
l <- lmer(v ~ Sex + (1|coil), gf, REML=F)
sink("/dev/null")
vsml <- mincLmer(jacobians_fixed_2 ~ Sex + (1|coil), gf, REML=F)
sink()
test_that("mincLmer basics", {
  expect_that(vsml[1,1], is_equivalent_to(fixef(l)[1]))
  expect_that(vsml[1,2], is_equivalent_to(fixef(l)[2]))
})

context("mincLmer - log likelihood ratios")
vsml2 <-  mincLmer(jacobians_fixed_2 ~ 1 + (1|coil), gf, REML=F)
l2 <- lmer(v ~ 1 + (1|coil), gf, REML=F)

test_that("logLikRatio", {
  expect_that( mincLogLikRatio(vsreml, vsml), throws_error() )
  expect_that( mincLogLikRatio(vsml, vsml2)[1],
              is_equivalent_to(anova(l,l2)[2,6]) )
})
