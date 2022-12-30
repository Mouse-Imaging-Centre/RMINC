library(testthat)
library(lme4)
context("vertexLmer")

if(!exists("dataPath"))
  dataPath <- tempdir()

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "RMINC-test-data-main/rminctestdata/")

gftest <- read.csv(file.path(dataPath, "subject.csv"), stringsAsFactors = TRUE)

subjectFile <- matrix(data = NA,nrow = 10,1)
subjectFile[1,1] <- file.path(dataPath, "vertex2.txt")
subjectFile[2,1] <- file.path(dataPath, "vertex3.txt")
subjectFile[3,1] <- file.path(dataPath, "vertex4.txt")
subjectFile[4,1] <- file.path(dataPath, "vertex3.txt")
subjectFile[5,1] <- file.path(dataPath, "vertex1.txt")
subjectFile[6,1] <- file.path(dataPath, "vertex2.txt")
subjectFile[7,1] <- file.path(dataPath, "vertex4.txt")
subjectFile[8,1] <- file.path(dataPath, "vertex2.txt")
subjectFile[9,1] <- file.path(dataPath, "vertex3.txt")
subjectFile[10,1] <- file.path(dataPath, "vertex1.txt")
gftest$testFilesLeft <- subjectFile

handle_conv_warnings <- function(expr){
  withCallingHandlers(expr, warning = function(w){
    if(grepl("converge", w)) invokeRestart("muffleWarning")
  })
}

test_env <- new.env()

test_that("Vertex REML Lmer Works", {
  evalq({
    vertex_table <- vertexTable(gftest$testFilesLeft)

    handle_conv_warnings({
      slow_lmer <-
        lapply(seq_len(nrow(vertex_table)), function(row){
          vals <- vertex_table[row,]
          lmer(vals ~ Age + (1 | Sex), data = gftest)
        })
    
      fast_lmer <-
        vertexLmer(testFilesLeft ~ Age + (1|Sex), data = gftest)
    })
    expect_equal(fast_lmer[,1:2], t(sapply(slow_lmer, fixef)), check.attributes = FALSE)
    expect_equal(fast_lmer[,3:4], t(sapply(slow_lmer, function(x) coefficients(summary(x))[,"t value"]))
                 , check.attributes = FALSE)
  }, envir = test_env)
})

test_that("Vertex ML Lmer Works", {
  evalq({
    handle_conv_warnings({
      fast_lmer2 <-
        vertexLmer(testFilesLeft ~ Age + (1|Sex), data = gftest, REML = FALSE)
    
      slow_lmer2 <-
        lapply(seq_len(nrow(vertex_table)), function(row){
          vals <- vertex_table[row,]
          lmer(vals ~ Age + (1 | Sex), data = gftest, REML = FALSE)
        })
    })
    
    expect_equal(fast_lmer2[,1:2], t(sapply(slow_lmer2, fixef)), check.attributes = FALSE, tolerance = 10e-5)
    expect_equal(fast_lmer2[,3:4], t(sapply(slow_lmer2, function(x) coefficients(summary(x))[,"t value"]))
                 , check.attributes = FALSE)
  }, envir = test_env)
})

test_that("Likelihood Ratio Tests for vertexLmer Work", {
  evalq({
    handle_conv_warnings({
      slow_lmer3 <-
        lapply(seq_len(nrow(vertex_table)), function(row){
          vals <- vertex_table[row,]
          lmer(vals ~ (1 | Sex), data = gftest, REML = FALSE)
        })
    
      fast_lmer3 <-
        vertexLmer(testFilesLeft ~ (1|Sex), data = gftest, REML = FALSE)
    })
    
    expect_error(mincLogLikRatio(fast_lmer, fast_lmer2), regexp = "REML=FALSE")
    expect_equal(mincLogLikRatio(fast_lmer2, fast_lmer3)[2]
                 , anova(slow_lmer2[[2]], slow_lmer3[[2]])[2,6]
                 , check.attributes = FALSE)
  }, envir = test_env)
})


context("vertexLmer - estimate DF")
test_that("empty DF by default", {
  evalq({
    expect_that( attr(fast_lmer, "df"), is_equivalent_to(NULL))
    expect_that( vertexFDR(fast_lmer), throws_error() )
  }, envir = test_env)
})


test_that("DF within reasonable range", {
  evalq({
    handle_conv_warnings({
      verboseRun(fast_lmer_df <- vertexLmerEstimateDF(fast_lmer))
    })
    df <- attr(fast_lmer_df, "df")
    expect_that( df[[2]], is_less_than(nrow(gftest)+1))
    expect_that( df[[2]], is_more_than(1))
  }, envir = test_env)
})

test_that("vertexLmer works with NAs", {
  verboseRun({
    gf_missing <- gftest
    gf_missing[1, "Sex"] <- NA

    handle_conv_warnings({
      missing <- vertexLmer(testFilesLeft ~ Age + (1|Sex), gf_missing)
      missing_dfs <- vertexLmerEstimateDF(missing)
    })
    df <- attr(missing_dfs, "df")
  })
  
  expect_that( df[[2]], is_less_than(nrow(attr(missing, "data")) + 1))
  expect_that( df[[2]], is_more_than(1))
})

