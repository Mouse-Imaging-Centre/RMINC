# Find Closest Vertex

Given a set of vertices (3-column matrix or data frame of x,y,z
coordinates) determine the index or coordinates of the closest vertex to
a given set of target vertices (either a 3 element vector, or x-y-z
table like structure)

## Usage

``` r
closestVertex(vertices, target, returns = c("index", "coordinates"))
```

## Arguments

- vertices:

  A matrix-like object with 3-columns, an n rows representing vertices

- target:

  either a 3-element numeric vector representing x-y-z coordinates for a
  single target, or a matrix-like object as described above containing
  multiple targets.

- returns:

  Whether to return the index of each match (one per target), or the
  coordinates of the matches, the later being useful when exact matches
  aren't expected.

## Value

In most cases, a numeric vector of match results, in the case of
multiple targets and `returns = "coordinates"` the output is a 3-column
matrix of coordinates.

## Details

Vertices are matched to the target by finding the vertex with the
minimum euclidean distance from the target.
