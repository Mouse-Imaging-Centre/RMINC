# Writes vertex data to a file with an optional header

Writes vertex data to a file with an optional header

## Usage

``` r
writeVertex(
  vertexData,
  filename,
  headers = TRUE,
  mean.stats = NULL,
  gf = NULL,
  col.names = FALSE
)
```

## Arguments

- vertexData:

  vertex data to be written

- filename:

  full path to file where data shall be written

- headers:

  Whether or not to write header information (implies col.names=TRUE)

- mean.stats:

  mean vertex data that may also be written

- gf:

  glim matrix that can be written to the header

- col.names:

  if column names should be written

## Value

A file is generated with the vertex data and optional headers

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
writeVertex(gf$nativeRMStlink20mm,"~/RMStlink20mm.txt",FALSE,NULL,NULL)
} # }
```
