library(testthat)
library(lme4)
context("vertexLmer")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gftest <- read.csv(file.path(dataPath, "subject.csv"))

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

test_env <- new.env()

test_that("Vertex REML Lmer Works", {
  evalq({
    vertex_table <- vertexTable(gftest$testFilesLeft)
    
    slow_lmer <-
      lapply(seq_len(nrow(vertex_table)), function(row){
        vals <- vertex_table[row,]
        lmer(vals ~ Age + (1 | Sex), data = gftest)
      })
    
    fast_lmer <-
      vertexLmer(testFilesLeft ~ Age + (1|Sex), data = gftest)
    
    expect_equal(fast_lmer[,1:2], t(sapply(slow_lmer, fixef)), check.attributes = FALSE)
    expect_equal(fast_lmer[,3:4], t(sapply(slow_lmer, function(x) coefficients(summary(x))[,"t value"]))
                 , check.attributes = FALSE)
  }, envir = test_env)
})

test_that("Vertex ML Lmer Works", {
  evalq({
    fast_lmer2 <-
      vertexLmer(testFilesLeft ~ Age + (1|Sex), data = gftest, REML = FALSE)
    
    slow_lmer2 <-
      lapply(seq_len(nrow(vertex_table)), function(row){
        vals <- vertex_table[row,]
        lmer(vals ~ Age + (1 | Sex), data = gftest, REML = FALSE)
      })
    
    expect_equal(fast_lmer2[,1:2], t(sapply(slow_lmer2, fixef)), check.attributes = FALSE, tolerance = 10e-5)
    expect_equal(fast_lmer2[,3:4], t(sapply(slow_lmer2, function(x) coefficients(summary(x))[,"t value"]))
                 , check.attributes = FALSE)
  }, envir = test_env)
})

test_that("Likelihood Ratio Tests for vertexLmer Work", {
  evalq({
    slow_lmer3 <-
      lapply(seq_len(nrow(vertex_table)), function(row){
        vals <- vertex_table[row,]
        lmer(vals ~ (1 | Sex), data = gftest, REML = FALSE)
      })
    
    fast_lmer3 <-
      vertexLmer(testFilesLeft ~ (1|Sex), data = gftest, REML = FALSE)
    
    expect_error(mincLogLikRatio(fast_lmer, fast_lmer2), regexp = "REML=FALSE")
    expect_equal(mincLogLikRatio(fast_lmer2, fast_lmer3)[2]
                 , anova(slow_lmer2[[2]], slow_lmer3[[2]])[2,6]
                 , check.attributes = FALSE)
  }, envir = test_env)
})


