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

test_that("anatGetAll_old works", {
  evalq({
    has_mincstuffs <- as.character(Sys.which("label_volumes_from_jacobians")) != ""
    skip_if_not(has_mincstuffs)
    
    vox_vol <- prod(minc.separation.sizes(segmentation))
    
    label_counts <-
      anatGetAll_old(segmentation,
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
      anatGetAll_old(gf$jacobians_0.2, 
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
      anatGetAll_old(gf$jacobians_0.2, 
                     atlas = segmentation, 
                     method = "means",
                     defs = labels)
    
    expect_equal(sweep(label_sums, 2, label_counts / vox_vol, FUN = `/`), 
                 label_means, tol = 10e-5, ignore.attributes = TRUE)
    #curious to note computing counts then sums and sweeping dividing
    #out the volume is faster than anatGetAll with means...
    
    jacobians <-
      anatGetAll_old(gf$jacobians_0.2, 
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

test_that("AnatGetAll works", {
  evalq({
    skip_if_not(has_mincstuffs)
    
    verboseRun(
      new_label_counts <- anatGetAll(segmentation, defs = labels, method = "labels")
    )
    
    expect_equal(new_label_counts[1,], label_counts_ref$count, check.attributes = FALSE)
    expect_equal(colnames(new_label_counts), colnames(label_counts))
    
    verboseRun(
      new_jacobians <- anatGetAll(gf$jacobians_0.2, defs = labels, atlas = segmentation, method = "jacobians")
    )
    
    expect_equal(new_jacobians, jacobians, check.attributes = FALSE, tolerance = 10e-4)
    expect_equal(colnames(new_jacobians), colnames(jacobians))
    
    verboseRun(
      new_sums <- anatGetAll(gf$jacobians_0.2, defs = labels, atlas = segmentation, method = "sums")
    )
    
    expect_equal(new_sums, label_sums, check.attributes = FALSE, tolerance = 10e-4)
    expect_equal(colnames(new_sums), colnames(label_sums))
    
    verboseRun(
      new_means <- anatGetAll(gf$jacobians_0.2, defs = labels, atlas = segmentation, method = "means")
    )
    
    expect_equal(new_means, label_means, check.attributes = FALSE, tolerance = 10e-4)
    expect_equal(colnames(new_means), colnames(label_means))
    
  }, envir = test_env)
})

test_that("AnatGetAll Flags Garbage", {
  gf2 <- gf
  gf2[1,"jacobians_0.2"] <- NA
  expect_error(anatGetAll(gf2$jacobians_0.2, defs = labels, atlas = segmentation)
               , regex = "could not be read")
})

test_that("Multires Works", {
  xfm <- file.path(dataPath, "scale10.xfm")
  cat("MNI Transform File"
    , "Transform_Type = Linear;"
    , "Linear_Transform ="
    , "10 0 0 0"
    , "0 10 0 0"
    , "0 0 10 0;"
    , sep = "\n", file = xfm)


  first_file <- gf$jacobians_0.2[1]
  transformed_file <- file.path(dataPath, "scaled_jacobians.mnc")
  sprintf("mincresample %s -tfm_input_sampling -transformation %s -step 1 1 1 %s"
        , first_file
        , xfm
        , transformed_file) %>%
    system(ignore.stdout = TRUE, ignore.stderr = TRUE)

  vols <-
    anatGetAll(c(as.character(gf$jacobians[1:5]), transformed_file)
             , defs = labels, atlas = segmentation, strict = FALSE)

  expect_equivalent(vols[1,] * 1000, vols[6,])
})


context("anatomy hierarchy summarization")
## Test hierarchy
test_hier1 <- label_frame
test_hier1$hierarchy <- paste0("struc", rep(1:11, length.out = 34))
test_hier1_file <- file.path(dataPath, "test_hierarchy1.csv")
write.csv(test_hier1, test_hier1_file, row.names = FALSE)

test_that("Anatomy hierarchy works", {
  evalq({
    skip_if_not(has_mincstuffs)
    hier_res <- anatSummarize(anat = new_jacobians, summarize_by = "hierarchy", defs = test_hier1_file)
    expect_equal(ncol(hier_res), 24)
  }, envir = test_env)
})


## Test hierarchy
test_hier2 <- label_frame
test_hier2$other.colname <- paste0("struc", rep(1:8, length.out = 34))
test_hier2_file <- file.path(dataPath, "test_hierarchy2.csv")
write.csv(test_hier2, test_hier2_file, row.names = FALSE)

test_that("anatomy hierarchy works", {
  evalq({
    skip_if_not(has_mincstuffs)
    hier_res <- anatSummarize(anat = new_jacobians, summarize_by = "other.colname", defs = test_hier2_file)
    expect_equal(ncol(hier_res), 17)
  }, envir = test_env)
})
