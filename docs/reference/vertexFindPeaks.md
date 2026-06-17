# Vertex find peaks

Find extrema of set of values given a connectivity graph

## Usage

``` r
vertexFindPeaks(
  data_map,
  graph,
  mindist = 1,
  direction = c("both", "positive", "negative"),
  threshold = 0,
  output = c("mask", "indices")
)
```

## Arguments

- data_map:

  A vector or text file of values to search for peaks

- graph:

  An igraph object describing the connectivity, in most cases this will
  be a mesh graph produced by
  [obj_to_graph](https://mouse-imaging-centre.github.io/RMINC/reference/obj_to_graph.md)

- mindist:

  the search radius (in number of neighbours) to consider when finding a
  peak. A node is considered a peak if it is more extreme than any
  element in it's neighbourhood up to `mindist` edges away

- direction:

  Either positive, negative, or both indicating whether to consider
  maxima, minima or both.

- threshold:

  threshold above or below which to consider peaks. In the both
  direction case a two element vector may be passed indicating positive
  and negative thresholds respectively. If a single element is passed in
  the both direction case, it is treated as `c(threshold, -threshold)`.
  Defaults to 0.

- output:

  Either mask or indices, indicating whether a logical mask or a vector
  of indices is desired.

## Value

Either a logical vector indicating peaks or a vector of peak indices
determined by `output`

## Details

Note that if coordinates of peaks are desired, they can be accessed from
the parent `bic_obj` by `obj$vertex_matrix[,peaks]`
