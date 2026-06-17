# Apply a structure summary function across vertices

This is a wrapper around vertexApply with \`transpose\` set to \`TRUE\`
and \`fun\` wrapped in with \`tapply\` over atlas.

## Usage

``` r
vertexAtlasApply(
  filenames,
  atlas,
  fun,
  ...,
  mask = NULL,
  parallel = NULL,
  collate = simplify_masked,
  column = 1,
  atlas_column = 1
)
```

## Arguments

- filenames:

  vertex file names

- atlas:

  The atlas to use to summarize vertices.

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

- column:

  Which column to treat as the input from vertex files.

- atlas_column:

  If the atlas is a text file, which column holds the label

## Value

The a matrix with a row of results for each structure
