context("mincLme - basic test")

gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
library(lme4)
gf$v <- mincGetVoxel(gf$jacobians_fixed_2, 0, 0, 0)
gf$coil <- as.factor(gf$coil)
l <- lmer(v ~ Sex + (1|coil), gf)
sink("/dev/null")
vs <- mincLmer(jacobians_fixed_2 ~ Sex + (1|coil), gf)
sink()
test_that("mincLme basics", {
  expect_that(vs[1,1], is_equivalent_to(fixef(l)[1]))
  expect_that(vs[1,2], is_equivalent_to(fixef(l)[2]))
})
            
