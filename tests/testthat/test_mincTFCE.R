context("Single dim mincTFCE")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"))
first_file <- gf$jacobians_0.2[1]

first_vol <- mincGetVolume(first_file)
local_env <- environment()

#mincPlotSliceSeries(mincArray(tfce_vol))
#mincPlotSliceSeries(mincArray(`likeVolume<-`(vtfce, first_file)))

test_voxels <-
  c(604L, 1188L, 488L, 262L, 1778L, 979L, 1399L, 1421L, 263L, 1434L, 
    562L, 828L, 247L, 1556L, 816L, 1552L, 833L, 563L, 638L, 1180L, 
    1925L, 725L, 1084L, 573L, 1443L)

example_values <-
  c(-0.000171213208313786, 8.36697609919128e-05, -0.000144912590568723, 
    -0.00027140487386863, 0.000210153896640721, 8.36697609919128e-05, 
    -3.16195020479375e-05, 5.47700425465171e-05, -0.000144912590568723, 
    -4.47209250138795e-05, 0.000171216272042315, -4.47209250138795e-05, 
    5.47700425465171e-05, 3.16225657764666e-05, -8.94433818920235e-05, 
    8.36697609919128e-05, -8.94433818920235e-05, 4.47239887424086e-05, 
    -0.000137840429415466, 9.48646336008707e-05, -7.74581871265381e-05, 
    6.32435996886687e-05, -4.47209250138795e-05, -7.07119320172097e-05, 
    -5.47751264690862e-05)
  
test_that("minc single dim TFCE works",{
  skip_on_cran()
  skip_on_travis()
  
  evalq({
    tfce_vol <- mincTFCE(first_vol)
    expect_equal(example_values, tfce_vol[test_voxels])
  }, envir = local_env)
})

context("Multi dim TFCE")

test_multidim <-
  cbind(v1 = first_vol, v2 = first_vol, v3 = rep(0, length(first_vol))) %>%
  as.minc %>%
  `likeVolume<-`(first_file) 



test_that("minc multi dim TFCE works", {
  skip_on_cran()
  skip_on_travis()
  
  evalq({
    md_tfce <- mincTFCE(test_multidim)
    expect_equal(md_tfce[test_voxels,1], example_values)
    expect_equal(md_tfce[test_voxels,2], example_values)
    expect_equal(md_tfce[test_voxels,3], rep(0, length(example_values)))
  }, envir = local_env)
})

context("MincLm Randomization TFCE")

## Really need a better test here
test_that("lm randomization works", {
  evalq({
    skip_on_cran()
    skip_on_travis()
    
    verboseRun(lmod <- mincLm(jacobians_0.2 ~ Genotype, gf))
    verboseRun(randomization_results <- mincTFCE(lmod, R = 10))
    
    thresh <- thresholds(randomization_results)["1%", ]
    maxes <- apply(randomization_results$tfce, 2, max)
    expect_true(between(thresh[1], 40, 100)) #Reasonable bounds??
    expect_true(between(thresh[2], .5, 5))    #Reasonable bounds??
    expect_lt(maxes[1], thresh[1])
    expect_gt(maxes[2], thresh[2])
  }, envir = local_env)
})

context("VertexTFCE numeric")

test_that("Vertex TFCE approximately matches mincTFCE", {
  skip_on_cran()
  skip_on_travis()
  
  tfce_vol <- mincTFCE(first_vol, d = .0027)
  vtfce <- vertexTFCE(first_vol, RMINC:::neighbour_list(15,15,15,6), weights = .001)
  expect_false(any(abs(vtfce - tfce_vol) > 1e-4))
})


# mincPlotSliceSeries(mincArray(tfce_vol))
# mincPlotSliceSeries(mincArray(`likeVolume<-`(vtfce, first_file)))
