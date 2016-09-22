suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
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

test_env <- new.env()

test_that("anatGetAll works", {
  evalq({
    has_mincstuffs <- as.character(Sys.which("label_volumes_from_jacobians")) != ""
    skip_if_not(has_mincstuffs)
    
    vox_vol <- prod(minc.separation.sizes(segmentation))
    
    label_counts <-
      anatGetAll(segmentation,
                 atlas = segmentation, 
                 method = "labels",
                 defs = labels)
    
    label_counts_ref <- 
      mincGetVolume(segmentation) %>%
      table %>%
      `*`(vox_vol) %>%
      as.data.frame %>%
      setNames(c("label", "count")) %>%
      filter(label %in% known_labels)
    
    expect_equal(label_counts[1,], label_counts_ref$count, check.attributes = FALSE)
    
    label_sums <-
      anatGetAll(gf$jacobians_0.2, 
                 atlas = segmentation, 
                 method = "sums",
                 defs = labels)
    
    seg_vol <- mincGetVolume(segmentation)
    
    ref_sums <-
      sapply(gf$jacobians_0.2, mincGetVolume) %>%
      apply(2, function(vol) tapply(vol, seg_vol, sum)) %>%
      t %>%
      .[,colnames(.) %in% known_labels]
    
    expect_equal(unclass(label_sums), ref_sums, tolerance = 10e-5, check.attributes = FALSE)
    
    label_means <-
      anatGetAll(gf$jacobians_0.2, 
                 atlas = segmentation, 
                 method = "means",
                 defs = labels)
    
    expect_equal(sweep(label_sums, 2, label_counts / vox_vol, FUN = `/`), 
                 label_means, tol = 10e-5, ignore.attributes = TRUE)
    #curious to note computing counts then sums and sweeping dividing
    #out the volume is faster than anatGetAll with means...
    
    jacobians <-
      anatGetAll(gf$jacobians_0.2, 
                 atlas = segmentation, 
                 method = "jacobians",
                 defs = labels)
    
    ref_jacobians <-
      sapply(gf$jacobians_0.2, mincGetVolume) %>%
      apply(2, function(vol) tapply(vol, seg_vol, function(j) sum(exp(j) * vox_vol))) %>%
      t %>%
      .[,colnames(.) %in% known_labels]
    
    expect_equal(ref_jacobians, 
                 unclass(jacobians), tol = 10e-5, check.attributes = FALSE)
  }, envir = test_env)
})

test_that("AnatGetAll2 works", {
  evalq({
    new_jacobians <- anatGetAll2(gf$jacobians_0.2, defs = labels, atlas = segmentation, method = "jacobians")
    expect_equal(new_jacobians, jacobians, check.attributes = FALSE, tolerance = 10e-4)
    
    new_label_counts <- anatGetAll2(gf$jacobians_0.2)
  }, envir = test_env)
})



