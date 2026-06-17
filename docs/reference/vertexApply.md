# Apply function over vertex Files

This function is used to compute an arbitrary function of every region
in a set of vertex files.

## Usage

``` r
vertexApply(
  filenames,
  fun,
  ...,
  mask = NULL,
  parallel = NULL,
  collate = simplify_masked,
  transpose = FALSE,
  column = 1
)
```

## Arguments

- filenames:

  vertex file names

- fun:

  A function to be applied to each vertex

- ...:

  additional arguments to `fun`

- mask:

  A vector of filename indicating a vertex mask (`fun` is applied to all
  vertices where mask is greater than .5)

- parallel:

  A two component vector indicating how to parallelize the computation.
  If the first element is "local" the computation will be run via the
  parallel package, otherwise it will be computed using batchtools, see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/reference/pMincApply.md)
  for details. The element should be numeric indicating the number of
  jobs to split the computation into.

- collate:

  A function to reduce the (potentially masked) list of results into a
  nice structure. Defaults to
  [simplify_masked](https://mouse-imaging-centre.github.io/RMINC/reference/simplify_masked.md)

- transpose:

  Whether to alternatively transpose the vertex matrix and apply a
  function to each subject

- column:

  Which column to treat as the input from vertex files.

## Value

The a matrix with a row of results for each vertex

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
vm <- vertexApply(gf$CIVETFILES$nativeRMStlink20mmleft, mean)
} # }
```
