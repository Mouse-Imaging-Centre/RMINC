suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(lmerTest)
})

context("anatLmer")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"))
first_file <- gf$jacobians_0.2[1]
segmentation <- file.path(dataPath, "test_segmentation.mnc")
labels <- file.path(dataPath, "test_defs.csv")
label_frame <- read.csv(labels)
known_labels <- with(label_frame, union(left.label, right.label))

anat_env <- new.env()

test_that("anatLmer works", {
  evalq({
    jacobians <-
      anatGetAll(gf$jacobians_0.2, 
                 atlas = segmentation, 
                 method = "jacobians",
                 defs = labels)
    
    suppressWarnings(
      lmer_res_nlhs <- anatLmer( ~ Pain.sensitivity + (1|Genotype), 
                                 data = gf, anat = jacobians))
    
    suppressWarnings(
      lmer_res_lhs <- anatLmer(dummy ~ Pain.sensitivity + (1|Genotype), 
                               data = gf, anat = jacobians))
    
    suppressWarnings({
      lmer_ref <-
        apply(jacobians, 2, function(col){
          gf2 <- gf %>% mutate(resp = col)
          mod <- lmer(resp ~ Pain.sensitivity + (1|Genotype), data = gf2)
          RMINC:::fixef_summary(mod)
        }) %>%
        t
    })
    
    expect_equal(unclass(lmer_res_nlhs), lmer_ref, tolerance = 10e-5, check.attributes = FALSE)
    expect_equal(lmer_res_lhs, lmer_res_nlhs, tolerance = 10e-5, check.attributes = FALSE)
  }, envir = anat_env)
})

test_that("anatLmer estimate DF returns sensible results", {
  evalq({
    expect_warning(verboseRun(with_dfs <- anatLmerEstimateDF(lmer_res_nlhs)),
                   regex = "Unable to estimate")
    dfs <- attr(with_dfs, "df")
    expect_true(between(dfs[1], 1, 3))
    expect_true(all(between(dfs[2:4], 30, 50)))            
  }, envir = anat_env)
})

test_that("anatLmer exotic formulae work", {
  evalq({
    verboseRun(
      exotic_lmer <- anatLmer(~ I(factor(as.numeric(Pain.sensitivity) - 1)) + (1 | Genotype)
                              , data = gf, anat = jacobians)
    )
    
    expect_warning(verboseRun(with_dfs <- anatLmerEstimateDF(exotic_lmer)),
                   regex = "Unable to estimate")
    dfs <- attr(with_dfs, "df")
    expect_true(between(dfs[1], 1, 3))
    expect_true(all(between(dfs[2:4], 30, 50)))
  }, envir = anat_env)
})