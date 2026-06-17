# Compute Laplace-Beltrami operator

Discrete Laplace-Beltrami operator associated with a triangle-mesh
manifold

## Usage

``` r
laplace_beltrami_operator(
  manifold = NULL,
  vertex_matrix = NULL,
  triangle_matrix = NULL
)
```

## Arguments

- manifold:

  A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
  like `bic_obj` object produced by
  [read_obj](https://mouse-imaging-centre.github.io/RMINC/dev/reference/read_obj.md)

- vertex_matrix:

  3-by-N matrix denoting the position of the N vertices defining the
  triangle-mesh manifold. Similar to argument `vertices` in
  [rgl::tmesh3d](https://dmurdoch.github.io/rgl/dev/reference/mesh3d.html).
  Required if `manifold` argument is not specified.

- triangle_matrix:

  3-by-M matrix denoting the position of the M triangles defining the
  triangle-mesh manifold. Elements of each column are indices of
  vertices defining the triangle (i.e. indices are columns of the
  vertex_matrix). Similar to argument `indices` in
  [rgl::tmesh3d](https://dmurdoch.github.io/rgl/dev/reference/mesh3d.html).
  Required if `manifold` argument is not specified.

## Value

sparse square matrix representing the discrete Laplace-Beltrami operator
for the manifold

## Details

Either supply `manifold` argument OR `vertex_matrix` and
`triangle_matrix` arguments.
