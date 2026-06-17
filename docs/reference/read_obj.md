# Read BIC .obj files

Read the BIC .obj 3D file format. This parses simple obj files to be
used with RMINC 3D plotting functions.

## Usage

``` r
read_obj(bic_obj, use_civet_triangles = FALSE)
```

## Arguments

- bic_obj:

  character file name of the obj file to read in

- use_civet_triangles:

  logical, whether or not to use the predefined triangle matrix common
  to obj files produced by CIVET 1.1.12, saves IO time when triangles
  are known in advance.

## Value

A two element list of class `bic_obj` containing a 3xV `vertex_matrix`
denoting the global coordinates of the each vertex and a 3xT
`triangle_matrix` containing triples of indices to the vertex matrix
representing individual triangles.

## Details

This parser is not robust at all and relies on a strict structure for
the .obj file at the present. It must be organized with a block of
vertices, seperated by a space from a block of colour information (which
is ignored), a space separated block of metadata and a space seperated
block of multiples of 3, followed finally by a block of triangle
membership. Only the vertex and triangle blocks are read in.
