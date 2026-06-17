# Add opacity to a mesh

Set the opacity for a brain mesh from a vector of values that correspond
to the vertices in the mesh.

## Usage

``` r
add_opacity(mesh, a_map, a_range = c(0.5, 1), a_default = 1)
```

## Arguments

- mesh:

  The brain mesh of interest

- a_map:

  The vector of values to be used to set alpha (opacity)

- a_range:

  The range of alpha values to be used in the image (after rescaling),
  must be between 0 and 1.

- a_default:

  the default alpha value for missing values in a_map

## Value

The original mesh with the alpha levels set
