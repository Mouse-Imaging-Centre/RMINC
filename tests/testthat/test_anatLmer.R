suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(lmerTest)
})

context("anatGet")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"))
first_file <- gf$jacobians_0.2[1]
segmentation <- file.path(dataPath, "test_segmentation.mnc")
labels <- file.path(dataPath, "test_defs.csv")
label_frame <- read.csv(labels)
known_labels <- with(label_frame, union(left.label, right.label))

test_that("anatLmer works", {
  
  has_mincstuffs <- Sys.which("label_volumes_from_jacobians") != ""
  skip_if_not(has_mincstuffs)
  
  jacobians <-
    anatGetAll(gf$jacobians_0.2, 
               atlas = segmentation, 
               method = "jacobians",
               defs = labels)
  
  lmer_res_nlhs <- anatLmer( ~ Pain.sensitivity + (1|Genotype), 
                             data = gf, anat = jacobians)
  
  lmer_res_lhs <- anatLmer(dummy ~ Pain.sensitivity + (1|Genotype), 
                           data = gf, anat = jacobians)
  
  lmer_ref <-
    apply(jacobians, 2, function(col){
      gf2 <- gf %>% mutate(resp = col)
      mod <- lmer(resp ~ Pain.sensitivity + (1|Genotype), data = gf2)
      RMINC:::mincLmerExtractVariables(mod)
    }) %>%
    t
  
  expect_equal(unclass(lmer_res_nlhs), lmer_ref, tolerance = 10e-5, check.attributes = FALSE)
  expect_equal(lmer_res_lhs, lmer_res_nlhs, tolerance = 10e-5, check.attributes = FALSE)
})