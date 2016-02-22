#' Read BIC .obj files
#' 
#' Read the BIC .obj 3D file format. This parses simple obj files to be
#' used with RMINC 3D plotting functions.
#' 
#' @param bic_obj character file name of the obj file to read in
#' @param use_civet_triangles logical, whether or not to use the predefined
#' triangle matrix common to obj files produced by CIVET 1.1.12, saves IO
#' time when triangles are known in advance.
#' @return A two element list of class \code{bic_obj} containing a 3xV \code{vertex_matrix}
#' denoting the global coordinates of the each vertex and a 3xT \code{triangle_matrix} containing
#' triples of indices to the vertex matrix representing individual triangles. 
#' @details This parser is not robust at all and relies on a strict structure for the .obj
#' file at the present. It must be organized with a block of vertices, seperated by a 
#' space from a block of colour information (which is ignored), a space separated block
#' of metadata and a space seperated block of multiples of 3, followed finally by a block
#' of triangle membership. Only the vertex and triangle blocks are read in.
#' @export
read_obj <- function(bic_obj, use_civet_triangles = FALSE) {
  
  lines <- readLines(bic_obj)
  
  general_info <- lines[1]
  
  ## Hacky parsing BIC .obj format see https://github.com/BIC-MNI/bicpl/Documentation
  ## This probably only works with output from CIVET 1.1.12
  section_ends <- 
    lines %>%
    `==`("") %>%
    which %>%
    as.list %>%
    setNames(c("vertices", "colours", "metadata", "multiples_of_three", "triangles"))
  
  parse_numbers <- 
    function(lines){
      lines %>%
        strsplit(" ") %>%
        unlist %>%
        Filter(f = function(chr) chr != "" ) %>%
        as.numeric
    } 
  
  vertices <- #Pull out vertices parse as 3xN matrix
    lines[2:section_ends$vertices] %>%  
    parse_numbers %>% 
    matrix(nrow = 3)
    
  if(!use_civet_triangles){
    polygons <- #Pull out triangles, parse as 3xN matrix
      lines[-(1:(section_ends$multiples_of_three))] %>% 
      parse_numbers %>%
      `+`(1) %>%
      matrix(nrow = 3)
  } else {
    polygons <- civet_polygons
  }  
  
  structure(list(vertex_matrix = vertices, triangle_matrix = polygons), class = "bic_obj")
}

add_mesh <-
  function(mesh){
    if(!is.null(mesh$id)) 
      
    shade3d(mesh)
    rgl_id <- rgl.ids()$id %>% tail(1)
    mesh$id <- rgl_id
  }


#' Create surface mesh for a bic_obj
#'
#' Turn the bic_obj object into an rgl tmesh ready for visualization
#' 
#' @param bic_obj A bic_obj to be converted to a mesh
#' @param color a character indicating the colour of triangles in the mesh
#' @param specular a character indicating the colour of the specular the lighting produces
#' on the object, "black" prevents the objects from looking shiny
#' @param add_normals logical whether to add normals to the mesh to make the surface less jagged
#' @param ... extra parameters passed to the material argument list of \link{tmesh3d} can include
#' additional \link{rgl.material} parameters to give to your mesh
#' 
#' @return a \link{obj_mesh} descended from \link{mesh3d} object to be 
#' plotted alone or subsequently colourized with \link{colour_mesh}
#' @export
create_mesh <- 
  function(bic_obj, color = "grey", specular = "black", add_normals = TRUE, ...){
    mesh <- with(bic_obj, tmesh3d(vertex_matrix, triangle_matrix, 
                                  material  = list(specular = specular,
                                                   color = color,
                                                   ...),
                                  homogeneous = FALSE))
    
    if(add_normals) mesh <- mesh %>% addNormals()
    
    class(mesh) <- c("obj_mesh", class(mesh))
    
    mesh
  }

#' Colourize a mesh
#' 
#' Add colour information to your mesh, either from a vertex atlas like AAL
#' or from a statistic/measurement map like those produced by CIVET
#' 
#' @param mesh \link{mesh3d} object ideally produced by \link{create_mesh}
#' @param colour_map either a vector with a label/measure/statistic for every vertex
#' or a character file path pointing to a file with such a vector in a rowwise format.
#' @param colour_range a two element numeric vector indicating the min and max values of 
#' allowable labels/measures/statistics to be includedon the surface
#' @param colour_default The colour given to vertices excluded by colour_range
#' @param reverse Whether to have a positive and negative colour scale (not yet implemented)
#' @param an \code{obj_mesh} object descended from \link{mesh3d}, with added colour information
#' and an additional \code{legend} element to be used in building a colour bar
#' @export  
colour_mesh <- function(mesh, 
                        colour_map,
                        colour_range = NULL,
                        colour_default = "grey",
                        reverse = NULL,
                        palette = heat.colors(255)){
  
  #Check colour_map is a numeric vector of colours per vertex or file name
  #Of a file containing such a vector spread over lines
  if(is.character(colour_map)){
    colour_map <- readLines(colour_map) %>% as.numeric
  } else if(!is.numeric(colour_map)){
    stop("Colour_map must either be a vector or a text file path")
  }
  
  if(is.null(colour_range)) colour_range <- range(colour_map)
  colour_map[!between(colour_map, colour_range[1], colour_range[2])] <- NA
  
  colour_depth <- length(palette)
  
  colour_indices <- 
    floor(
      (colour_map - min(colour_map, na.rm = TRUE)) / 
        diff(range(colour_map, na.rm = TRUE)) * (colour_depth - 1)) + 1
  
  #Internally in tmesh3d, the vertex matrix is expanded
  #Such that the vertices for each triangle appear sequentially
  #in informal groups of three, the colours need to be triplicated
  #and reordered to match, this is acheived by using the polygon matrix 
  #as an index to the colours
  colours <- palette[colour_indices][mesh$it]
  
  mesh$legend <- list(colour_range = colour_range, 
                      palette = palette)
  
  mesh$material$color <- colours
  
  class(mesh) <- c("obj_mesh", class(mesh))
  
  mesh
}


#' Plot a BIC obj
#' 
#' Create a basic plot of a BIC obj with a uniformly coloured mesh in the
#' current rgl device, opening a new one if necessary.
#' 
#' @param x A \code{bic_obj} probably created by \link{read_obj}
#' @param ... additional arguments to \link{create_mesh} including but not limited
#' to colour, specular, and add_normals
#' @return Invisibly returns the created mesh
#' @export
plot.bic_obj <- 
  function(x, ...){
   mesh <- create_mesh(x, ...) 
   mesh %>% shade3d()
   invisible(mesh)
  }

#' Plot an BIC obj mesh
#' 
#' Create a plot of BIC obj_mesh, potentially colourized by \link{colour_mesh}
#' and potentially including a colour bar
#' 
#' @param x a \link{obj_mesh} object
#' @param colour_bar whether or not to add a colour bar
#' @param additional parameters to pass to add_colour_bar
#' @return returns x invisibly
#' @export
plot.obj_mesh <-
  function(x, colour_bar = TRUE, ...){
    x %>% shade3d
    if(colour_bar) x %>% add_colour_bar
    
    invisible(x)
  }

#' Plot a BIC obj with a colour map
#' 
#' Create a colour mapped bic_obj and plot it in the current rgl
#' device, opening a new one if necessary.
#' 
#' @param bic_obj A \code{bic_obj} probably created by \link{read_obj}
#' @param colour_map A numeric vector equal in length to the number of vertices
#' in the \code{bic_obj} or the path to a text file with one line per vertex with
#' colour information.
#' @param A palette, AKA look-up-table, providing a linear colour scale for the colours in
#' \code{colour_map}
#' @param ... additional arguments to \link{create_mesh} including but not limited
#' to colour, specular, and add_normals
#' @param invisibly return the colourized mesh object
#' @export
plot_bic_obj_overlay <- 
  function(bic_obj, 
           colour_map, 
           palette = heat.colors(255),
           colour_bar = TRUE,
           ...){
    
    mesh <-
      bic_obj %>%
      create_mesh(...)  %>%
      colour_mesh(colour_map)
      
    mesh %>% shade3d(override = FALSE)
    if(colour_bar) mesh %>% add_colour_bar
    
    invisible(mesh)
  }

#' Add a colour bar for a mesh
#' 
#' Add a colour bar that corresponds to the colours in a given colourized mesh
#' to the current rgl device
#' 
#' @param mesh A \code{obj_mesh} object created with \link{colour_mesh}
#' @param title The label to give the colour bar
#' @param ... additional plot parameters to be passed down to the \link{image}
#' function.
#' @param column_widths define the horizonatal segmentation of the plot device,
#' values are scaled to sum to 1. Defaults to a 15% region on the left to hold
#' the colour bar
#' @param which_col The column index of where to draw the colour bar, defaults
#' to the last panel, useful for creating a non-standard colour bar location.
#' @return invisible NULL
#' @export
add_colour_bar <- function(mesh,
                           title = "",
                           ...,
                           column_widths = c(.75, .25), 
                           which_col = length(column_widths)){
  if(is.null(mesh$legend)) stop("Your mesh has not colour information")
  with(mesh$legend, {
    bgplot3d({
      par(mar = c(4,8,4,2))
      
      stat_range <- seq(colour_range[1], colour_range[2], 
                        length.out = length(palette))
      
      layout(matrix(1:length(column_widths), nrow = 1), widths = column_widths)
      replicate(which_col - 1, plot.new())
      
      image(y = stat_range, 
            z =  matrix(stat_range, nrow = 1), 
            col = palette,
            xaxt = "n",
            ylab = title,
            useRaster = TRUE,
            ...)
    })
  })
  
  invisible(NULL)
}

#' Brain Montage
#' 
#' Plot the left and right hemispheres of brain colourized by two colour maps
#' with a colour bar in the middle
#' 
#' @param left_obj  A \code{bic_obj} probably created by \link{read_obj} with the left hemisphere
#' of a subject's brain
#' @param right_obj  A \code{bic_obj} probably created by \link{read_obj} with the right hemisphere
#' of a subject's brain
#' @param left_map a colour map to apply to the left hemisphere see \link{colour_mesh} for details
#' @param right_map a colour map to apply to the right hemisphere see \link{colour_mesh} for details
#' @param output Either NULL or a file path to write the snapshot.
#' @param ... additonal parameters to be passed to \link{create_mesh}
#' @param add_normals Whether or not to add normals to the surface objects, see \link{create_mesh} for
#' details
#' @param close_on_output Whether or not to close the output after taking a snapshot, defaults to
#' TRUE
#' @details
#' This function is designed to do a simple 6 angle plot for statistic maps of the left and right
#' hemispheres of subject's brain. It defaults to leaving the rgl device open so that you can
#' take a snapshot after tweaking the angles, or changing the colour bar with \link{add_colour_bar}.
#' Its other mode is to take a snapshot after it has finished adding the 3d objects to the scene,
#' in this mode, whether or not to keep the window open can be configured with close_on_output.
#' @return invisible NULL
#' @export
obj_montage <- function(left_obj, 
                        right_obj, 
                        left_map, 
                        right_map,
                        output = NULL,
                        ...,
                        add_normals = TRUE,
                        colour_title = "",
                        close_on_output = TRUE){
  left_mesh <- 
    create_mesh(left_obj, add_normals = TRUE, ...) %>% 
    colour_mesh(left_map)
  
  right_mesh <- 
    create_mesh(right_obj, add_normals = TRUE, ...) %>% 
    colour_mesh(right_map)
  
  open3d(windowRect = c(100,100, 900, 900))
  #This isn't my fault I swear, see the ?bgplot3d examples, weird bugs happen
  #without waiting, e.g. brains will only draw when there is an inactive rgl device
  Sys.sleep(.25)
  par3d(viewport = c(0,0,800,800))

  parent_scene <- currentSubscene3d()
  
  subscenes <- mfrow3d(3,2)
  left_or_right <- rep(c("left", "right"), 3)
  
  mapply(function(subscene, view_matrix, hemisphere){
    useSubscene3d(subscene)
    if(hemisphere == "left"){
      left_mesh %>% shade3d(override = FALSE)
    } else {
      right_mesh %>% shade3d(override = FALSE)
    }
    
    par3d(userMatrix = view_matrix)
  }, subscenes, hemisphere_viewpoints, left_or_right)

  useSubscene3d(parent_scene)
  
  left_mesh %>% add_colour_bar(title = colour_title, column_widths = c(6,2,6), which_col = 2, cex.axis = 2, cex.lab = 2)
  
  if(!is.null(output)) snapshot3d(output)
  if(!is.null(output) && close_on_output) rgl.close()
  
  invisible(NULL)
}

