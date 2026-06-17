# Vertex Lookup

Find the vertices closest one or more targets, potentially returning the
values for the vertices from a data map.

## Usage

``` r
vertexLookup(
  vertices,
  target,
  data_map = NULL,
  returns = c("index", "coordinates"),
  coerce = as.numeric
)
```

## Arguments

- vertices:

  A descendant of
  [mesh3d](https://dmurdoch.github.io/rgl/dev/reference/mesh3d.html),
  `bic_obj`, or matrix-like object with 3-columns, and n rows
  representing vertices.

- target:

  either a 3-element numeric vector representing x-y-z coordinates for a
  single target, or a matrix-like object as described above containing
  multiple targets.

- data_map:

  Either NULL, a vector of data about each vertex, or a file containing
  such a vector spread over multiple lines

- returns:

  Whether to return the index of each match (one per target), or the
  coordinates of the matches, the later being useful when exact matches
  aren't expected.

- coerce:

  A function to coerce the final results to a given type. Defaults to
  [as.numeric](https://rdrr.io/r/base/numeric.html), if set to NULL, no
  coercion is performed.

## Value

If a data_map is specified: a vector, typically numeric, if coerce is
set to NULL and data_map is a file, the results will be character. If
coerce is null and data_map is a vector it will return the same type as
data_map. If data_map is unspecified, it acts like
[closestVertex](https://mouse-imaging-centre.github.io/RMINC/reference/closestVertex.md)
