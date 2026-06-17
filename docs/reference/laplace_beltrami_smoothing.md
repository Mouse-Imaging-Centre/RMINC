# Smoothing on manifold

Smoothing scalar field on a triangle-mesh manifold

## Usage

``` r
laplace_beltrami_smoothing(manifold, scalar_field, fwhm, maxiter = 1000)
```

## Arguments

- manifold:

  A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
  like `bic_obj` object produced by
  [read_obj](https://mouse-imaging-centre.github.io/RMINC/reference/read_obj.md)

- scalar_field:

  a vector whose elements are values of the field at each manifold
  vertex

- fwhm:

  the degree of smoothing. It represents the full-width-half-maximum of
  the blurring kernel if the manifold had no curvature. Regardless of
  curvature, higher the fwhm, the greater the amount of smoothing.

- maxiter:

  int specifying the maximum number of iterations for the algorithm

## Value

a vector whose elements are values of the smoothed field at each
manifold vertex

## Examples

``` r
if (FALSE) { # \dontrun{
# Load an object
manifold =
  read_obj(
    file.path("/axiom2/projects/software/cortical-thickness/",
              "MWM/c57bl6_laplacian_grid_full_surface_simplified.obj"))

# Compute the laplace beltrami operator and attach it to the manifold
#  not necessary but will make future smoothing computations on the same manifold faster
manifold$laplace_beltrami_operator = laplace_beltrami_operator(manifold)

# Generate some random uniform vertex data
init_stats = runif(ncol(manifold$vertex_matrix))

# Smooth the stats on the manifold
smooth_stats = laplace_beltrami_smoothing(manifold,init_stats,0.2)

# Plot results
plot(manifold, init_stats, colour_range = c(.5,1))
plot(manifold, smooth_stats, colour_range = c(.5,1))
} # }
```
