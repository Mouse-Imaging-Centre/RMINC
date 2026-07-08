# Read CIVET QC data

After running the CIVET quality control pipeline, import the results

## Usage

``` r
civet.readQC(dir, civetVersion = "1.1.12")
```

## Arguments

- dir:

  The CIVET QC directory, or vector of directories

- civetVersion:

  the version of CIVET used, currently only supports 1.1.12, 2.0.0,
  2.1.0

## Value

A table of QC results including whether or not the subjects passed
overall quality control. See
<http://www.bic.mni.mcgill.ca/ServicesSoftware/QualityControlCIVET12>
for more details.
