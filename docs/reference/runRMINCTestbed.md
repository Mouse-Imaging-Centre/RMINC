# Run Testbed

Run the test bed to ensure all RMINC functions work on your system

## Usage

``` r
runRMINCTestbed(
  ...,
  dataPath = getOption("RMINC_DATA_DIR", tempdir()),
  test_q_minc = TRUE,
  method = "libcurl",
  verboseTest = FALSE
)
```

## Arguments

- ...:

  additional parameter for
  [test_dir](https://testthat.r-lib.org/reference/test_dir.html)

- dataPath:

  The directory to download and unpack the test data (unpacks in
  dataPath/rminctestdata). Default can be set with the option
  RMINC_DATA_DIR which can in turn be set with the environment variable
  RMINC_DATA_DIR. If unset a temporary directory is created.

- test_q_minc:

  Whether or not to test the cluster code, defaults to TRUE.

- method:

  Argument to pass to
  [download.file](https://rdrr.io/r/utils/download.file.html) typical
  options are `libcurl`

- verboseTest:

  Whether or not to verbosely print test output, default is to print
  simplified results

## Value

invisibly return the test results
