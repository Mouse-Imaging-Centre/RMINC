# Download Example Data

Download the example data needed to run our examples in your /tmp
directory The data can be downloaded manually from
<https://wiki.mouseimaging.ca/download/attachments/1654/rminctestdata.tar.gz>

## Usage

``` r
getRMINCTestData(
  dataPath = getOption("RMINC_DATA_DIR", tempdir()),
  method = "libcurl"
)
```

## Arguments

- dataPath:

  The directory to download and unpack the test data (unpacks in
  dataPath/rminctestdata). Default can be set with the option
  RMINC_DATA_DIR which can in turn be set with the environment variable
  RMINC_DATA_DIR. If unset a temporary directory is created.

- method:

  Argument to pass to
  [download.file](https://rdrr.io/r/utils/download.file.html) typical
  options are `libcurl` and `wget`
