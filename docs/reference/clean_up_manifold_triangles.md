# Clean up manifold mesh

If mesh edge has 2+ attached triangles, keep the 2 triangles most
well-connected triangles

## Usage

``` r
clean_up_manifold_triangles(manifold, k = 1)
```

## Arguments

- manifold:

  A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
  like `bic_obj` object produced by
  [read_obj](https://mouse-imaging-centre.github.io/RMINC/reference/read_obj.md)

- k:

  degrees-of-separation for defining neighbours. Default is 1 (i.e.
  adjacent vertices)

## Value

A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
similar to 'manifold' argument. Problematic triangles (if they exist)
are removed from the mesh.

## Details

increasing k drastically slows down computation
