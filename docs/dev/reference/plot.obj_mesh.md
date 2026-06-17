# Plot an BIC obj mesh

Create a plot of BIC obj_mesh, potentially colourized by
[colour_mesh](https://mouse-imaging-centre.github.io/RMINC/dev/reference/colour_mesh.md)
and potentially including a colour bar

## Usage

``` r
# S3 method for class 'obj_mesh'
plot(x, colour_bar = TRUE, ...)
```

## Arguments

- x:

  a `obj_mesh` object

- colour_bar:

  whether or not to add a colour bar

- ...:

  additional parameters to pass to add_colour_bar

## Value

returns x invisibly
