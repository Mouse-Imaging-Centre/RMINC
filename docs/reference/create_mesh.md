# Create surface mesh for a bic_obj

Turn the bic_obj object into an rgl tmesh ready for visualization

## Usage

``` r
create_mesh(
  bic_obj,
  color = "grey",
  specular = "black",
  add_normals = TRUE,
  ...
)
```

## Arguments

- bic_obj:

  A bic_obj to be converted to a mesh

- color:

  a character indicating the colour of triangles in the mesh

- specular:

  a character indicating the colour of the specular the lighting
  produces on the object, "black" prevents the objects from looking
  shiny

- add_normals:

  logical whether to add normals to the mesh to make the surface less
  jagged

- ...:

  extra parameters passed to the material argument list of tmesh3d can
  include additional rgl.material parameters to give to your mesh

## Value

a `obj_mesh` descended from mesh3d object to be plotted alone or
subsequently colourized with
[colour_mesh](https://mouse-imaging-centre.github.io/RMINC/reference/colour_mesh.md)
