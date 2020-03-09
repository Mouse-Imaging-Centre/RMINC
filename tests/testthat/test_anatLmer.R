suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(lme4)
})

handle_conv_warnings <- function(expr){
  withCallingHandlers(expr, warning = function(w){
    if(grepl("converge", w)) invokeRestart("muffleWarning")
  })
}

context("anatLmer")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

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
    
    handle_conv_warnings(
      lmer_res_nlhs <- anatLmer( ~ Pain.sensitivity + (1|Genotype), 
                                 data = gf, anat = jacobians))
    
    handle_conv_warnings(
      lmer_res_lhs <- anatLmer(dummy ~ Pain.sensitivity + (1|Genotype), 
                               data = gf, anat = jacobians))
    
    handle_conv_warnings({
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
    handle_conv_warnings({
      verboseRun(with_dfs <- anatLmerEstimateDF(lmer_res_nlhs))
    })
    dfs <- attr(with_dfs, "df")
    expect_true(between(dfs[1], 1, 3))
    expect_true(all(between(dfs[2:4], 20, 31)))            
  }, envir = anat_env)
})

test_that("anatLmer exotic formulae work", {
  evalq({
    handle_conv_warnings({
      verboseRun(
        exotic_lmer <- anatLmer(~ I(factor(as.numeric(Pain.sensitivity) - 1)) + (1 | Genotype)
                              , data = gf, anat = jacobians)
      )
    
      verboseRun(with_dfs <- anatLmerEstimateDF(exotic_lmer))
    })
    dfs <- attr(with_dfs, "df")
    expect_true(between(dfs[1], 1, 3))
    expect_true(all(between(dfs[2:4], 20, 31)))
  }, envir = anat_env)
})

test_that("anatLmer exotic formulae work", {
  evalq({
    handle_conv_warnings({
      verboseRun(
        exotic_lmer <- anatLmer(~ I(factor(as.numeric(Pain.sensitivity) - 1)) + (1 | Genotype)
                              , data = gf, anat = jacobians)
      )
    
      verboseRun(with_dfs <- anatLmerEstimateDF(exotic_lmer))
    })
    dfs <- attr(with_dfs, "df")
    expect_true(between(dfs[1], 1, 3))
    expect_true(all(between(dfs[2:4], 20, 31)))
  }, envir = anat_env)
})

test_that("weighted lmer works", {
    verboseRun({
        
        d <- data.frame(x = rnorm(20)
                      , w = runif(20)
                      , g = rep(c("A", "B"), length.out = 20))

        a <- matrix(rnorm(60), ncol = 3)

        handle_conv_warnings({
          weighted_lmer <- anatLmer( ~ x + (1 | g), data = d, anat = a, weights = d$w)
          
          lmer_ref <-
            apply(a, 2, function(col){
              d <- d %>% mutate(resp = col)
              mod <- lmer(resp ~ x + (1 | g), data = d, weights = d$w)
              RMINC:::fixef_summary(mod)
            }) %>%
            t
        })

        expect_equivalent(unclass(weighted_lmer), lmer_ref)
    })
})
      
