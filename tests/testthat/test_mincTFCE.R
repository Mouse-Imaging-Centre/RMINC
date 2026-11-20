context("Single dim mincTFCE")

if (!exists("dataPath")) {
  dataPath <- tempdir()
}

getRMINCTestData(dataPath)
dataPath <- file.path(dataPath, "rminctestdata/")

gf <- read.csv(file.path(dataPath, "minc_summary_test_data.csv"), stringsAsFactors = TRUE)
first_file <- gf$jacobians_0.2[1]

first_vol <- mincGetVolume(first_file)
local_env <- environment()

# mincTFCE runs the TFCE program from minc-stuffs
has_tfce <- function() {
  nzchar(Sys.which("TFCE")) &&
    system("TFCE --help", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
}

#mincPlotSliceSeries(mincArray(tfce_vol))
#mincPlotSliceSeries(mincArray(`likeVolume<-`(vtfce, first_file)))

test_voxels <-
  c(
    604L,
    1188L,
    488L,
    262L,
    1778L,
    979L,
    1399L,
    1421L,
    263L,
    1434L,
    562L,
    828L,
    247L,
    1556L,
    816L,
    1552L,
    833L,
    563L,
    638L,
    1180L,
    1925L,
    725L,
    1084L,
    573L,
    1443L
  )

example_values <-
  c(
    -0.000171213208313786,
    8.36697609919128e-05,
    -0.000144912590568723,
    -0.00027140487386863,
    0.000210153896640721,
    8.36697609919128e-05,
    -3.16195020479375e-05,
    5.47700425465171e-05,
    -0.000144912590568723,
    -4.47209250138795e-05,
    0.000171216272042315,
    -4.47209250138795e-05,
    5.47700425465171e-05,
    3.16225657764666e-05,
    -8.94433818920235e-05,
    8.36697609919128e-05,
    -8.94433818920235e-05,
    4.47239887424086e-05,
    -0.000137840429415466,
    9.48646336008707e-05,
    -7.74581871265381e-05,
    6.32435996886687e-05,
    -4.47209250138795e-05,
    -7.07119320172097e-05,
    -5.47751264690862e-05
  )

test_that("minc single dim TFCE works", {
  skip_on_cran()
  skip_on_travis()

  skip_if_not(has_tfce(), "TFCE from minc-stuffs is not installed")

  evalq(
    {
      tfce_vol <- mincTFCE(first_vol)
      expect_equal(example_values, tfce_vol[test_voxels])
    },
    envir = local_env
  )
})

context("Multi dim TFCE")

test_multidim <-
  cbind(v1 = first_vol, v2 = first_vol, v3 = rep(0, length(first_vol))) %>%
  as.minc %>%
  `likeVolume<-`(first_file)


test_that("minc multi dim TFCE works", {
  skip_on_cran()
  skip_on_travis()

  skip_if_not(has_tfce(), "TFCE from minc-stuffs is not installed")

  evalq(
    {
      md_tfce <- mincTFCE(test_multidim)
      expect_equal(md_tfce[test_voxels, 1], example_values)
      expect_equal(md_tfce[test_voxels, 2], example_values)
      expect_equal(md_tfce[test_voxels, 3], rep(0, length(example_values)))
      expect_equal(likeVolume(md_tfce), likeVolume(test_multidim))
    },
    envir = local_env
  )
})

context("MincLm Randomization TFCE")

## Really need a better test here
test_that("lm randomization works", {
  evalq(
    {
      skip_on_cran()
      skip_on_travis()

      skip_if_not(has_tfce(), "TFCE from minc-stuffs is not installed")

      verboseRun(lmod <- mincLm(jacobians_0.2 ~ Genotype, gf))
      # A fixed seed makes the permutations, and so the thresholds,
      # the same on each run (#212, #290)
      set.seed(1)
      verboseRun(randomization_results <- mincTFCE(lmod, R = 10))

      thresh <- thresholds(randomization_results)["1%", ]
      maxes <- apply(randomization_results$tfce, 2, max)
      # The old fixed bounds on the thresholds failed on some runs (#290).
      # Check only that the thresholds are usable values.
      expect_equal(dim(randomization_results$randomization_dist), c(10L, 2L))
      expect_true(all(is.finite(thresh)))
      expect_true(all(thresh > 0))
      expect_lt(maxes[1], thresh[1])
      expect_gt(maxes[2], thresh[2])
    },
    envir = local_env
  )
})

context("VertexTFCE numeric")

test_that("Vertex TFCE approximately matches mincTFCE", {
  skip_on_cran()
  skip_on_travis()

  skip_if_not(has_tfce(), "TFCE from minc-stuffs is not installed")

  tfce_vol <- mincTFCE(first_vol, d = .0027)
  vtfce <- vertexTFCE(
    first_vol,
    RMINC:::neighbour_list(15, 15, 15, 6),
    weights = .001
  )
  expect_false(any(abs(vtfce - tfce_vol) > 1e-4))
})

# mincPlotSliceSeries(mincArray(tfce_vol))
# mincPlotSliceSeries(mincArray(`likeVolume<-`(vtfce, first_file)))

context("mincTFCE checks the TFCE program")

# Put a fake TFCE program first on the PATH. `body` is the shell code
# that runs after `input` and `output` are set from the last two
# arguments.
with_fake_tfce <- function(body, code) {
  bin <- tempfile("fake_tfce_")
  dir.create(bin)
  writeLines(
    c(
      "#!/bin/sh",
      "for a; do input=$output; output=$a; done",
      body
    ),
    file.path(bin, "TFCE")
  )
  Sys.chmod(file.path(bin, "TFCE"), "0755")
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
  on.exit({
    Sys.setenv(PATH = old_path)
    unlink(bin, recursive = TRUE)
  })
  force(code)
}

test_that("mincTFCE stops when TFCE fails, even with an output file", {
  # pyminc can fail after it writes the data. Before, mincTFCE read the
  # output and did not report the failure.
  with_fake_tfce(
    c('cp "$input" "$output"', 'echo "fake TFCE failure" >&2', "exit 3"),
    expect_error(mincTFCE(first_vol), "exit status 3")
  )
})

test_that("mincTFCE stops when TFCE writes no output file", {
  with_fake_tfce(
    "exit 0",
    expect_error(mincTFCE(first_vol), "did not write its output file")
  )
})

test_that("mincTFCE stops when TFCE is not installed", {
  empty_bin <- tempfile("no_tfce_")
  dir.create(empty_bin)
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = empty_bin)
  on.exit({
    Sys.setenv(PATH = old_path)
    unlink(empty_bin, recursive = TRUE)
  })
  expect_error(mincTFCE(first_vol), "needs the TFCE program from minc-stuffs")
})

test_that("mincTFCE reads the output when TFCE succeeds", {
  with_fake_tfce(
    'cp "$input" "$output"',
    {
      out <- mincTFCE(first_vol)
      expect_equal(as.numeric(out), as.numeric(first_vol))
    }
  )
})
