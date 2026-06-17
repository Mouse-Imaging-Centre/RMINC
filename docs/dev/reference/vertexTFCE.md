# Threshold Free Cluster Enhancement

Perform threshold free cluster enhancement as described in Smith and
Nichols (2008). Cluster-like structures are enhanced to allow a hybrid
cluster/voxel analysis to be performed.

## Usage

``` r
vertexTFCE(x, ...)

# S3 method for class 'numeric'
vertexTFCE(
  x,
  surface,
  E = 0.5,
  H = 2,
  nsteps = 100,
  side = c("both", "positive", "negative"),
  weights = NULL,
  ...
)

# S3 method for class 'matrix'
vertexTFCE(
  x,
  surface,
  E = 0.5,
  H = 2,
  nsteps = 100,
  side = c("both", "positive", "negative"),
  weights = NULL,
  ...
)

# S3 method for class 'vertexLm'
vertexTFCE(
  x,
  surface,
  R = 500,
  alternative = c("two.sided", "greater"),
  E = 0.5,
  H = 2,
  nsteps = 100,
  weights = NULL,
  side = c("both", "positive", "negative"),
  replace = FALSE,
  parallel = NULL,
  ...
)

# S3 method for class 'character'
vertexTFCE(
  x,
  surface,
  E = 0.5,
  H = 2,
  nsteps = 100,
  side = c("both", "positive", "negative"),
  weights = NULL,
  column = 1,
  ...
)
```

## Arguments

- x:

  A numeric vector, a filepath to a set of values, or a `matrix` object,
  or `vertexLm` object.

- ...:

  additional arguments for methods

- surface:

  Either a mesh object corresponding to the surface, an igraph graph
  object of surface created by
  [obj_to_graph](https://mouse-imaging-centre.github.io/RMINC/dev/reference/obj_to_graph.md),
  or an adjacency list (see details). For the `matrix` and
  [vertexLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/vertexLm.md)
  cases, either a single surface object may be passed and used for each
  individual, or a vector of file names

- E:

  The exponent by which to raise the extent statistic (default .5)

- H:

  The exponent by which to raise the height (default 2)

- nsteps:

  The number of steps to discretize the TFCE computation over

- side:

  Whether to consider positive and negative statistics or both (default
  both)

- weights:

  A weighting vector assigning area to vertices. The default varies by
  `surface` object type. If `surface` is an `bic_obj` the area of each
  vertex is equal to 1/3 the sum of the areas of the triangles it is a
  member of. If `surface` is an `igraph` or adjacency list the default
  is a vector of ones.

- R:

  number of randomizations to perform

- alternative:

  Whether to consider a one-sided or two-sided alternative hypothesis.
  Default "two-sided", use "greater" for a one sided test.

- replace:

  Sample with or without replacement for the randomization, defaults to
  FALSE (no replacement)

- parallel:

  A two component vector indicating how to parallelize the computation.
  If the first element is "local" the computation will be run via the
  parallel package, otherwise it will be computed using batchtools, see
  [pMincApply](https://mouse-imaging-centre.github.io/RMINC/dev/reference/pMincApply.md)
  for details. The element should be numeric indicating the number of
  jobs to split the computation into.

- column:

  Which column to treat as the input from vertex files.

## Value

The behaviour of `vertexTFCE` is to perform cluster free enhancement on
a object, in the single dimensional case, a string denoting a vertex
file or a numeric vector it returns a numeric vector with the result. In
the matrix case each column is cluster enhanced and recomposed into a
matrix. In the
[vertexLm](https://mouse-imaging-centre.github.io/RMINC/dev/reference/vertexLm.md)
case a randomization test is performed on each t-tstatistic column. The
results is a list with 3 elements

- TFCE: A matrix of the tvalue columns after randomization

- randomization_dist: An RxT matrix where R is the number of
  randomizations and T is the number of t-statistics, elements are the
  largest value obtained by the randomized TFCE

- args: Arguments passed to the internal randomzations and TFCE code

## Details

Passing an adjacency list will save some compute time but is not
recommended for general use. If an adjacency list is passed should index
starting from 0 for compatibility with c++ code. Adjacency lists of this
kind can be generated from graphs with
`` lapply(as_adj_list(graph), `-`, 1) `` using the
[as_adj_list](https://r.igraph.org/reference/as_adj_list.html) from the
igraph library.

## Methods (by class)

- `vertexTFCE(numeric)`: numeric

- `vertexTFCE(matrix)`: matrix

- `vertexTFCE(vertexLm)`: vertexLm

- `vertexTFCE(character)`: character
