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


#' Create surface mesh for a bic_obj
#'
#' Turn the bic_obj object into an rgl tmesh ready for visualization
#' 
#' @param bic_obj A bic_obj to be converted to a mesh
#' @param color a character indicating the colour of triangles in the mesh
#' @param specular a character indicating the colour of the specular the lighting produces
#' on the object, "black" prevents the objects from looking shiny
#' @param add_normals logical whether to add normals to the mesh to make the surface less jagged
#' @param ... extra parameters passed to the material argument list of \link[rgl]{tmesh3d} can include
#' additional \link[rgl]{rgl.material} parameters to give to your mesh
#' 
#' @return a \code{obj_mesh} descended from \link[rgl]{mesh3d} object to be
#' plotted alone or subsequently colourized with \link{colour_mesh}
#' @export
create_mesh <- 
  function(bic_obj, color = "grey", specular = "black", add_normals = TRUE, ...){
    mesh <- with(bic_obj, rgl::tmesh3d(vertex_matrix, triangle_matrix, 
                                       material  = list(specular = specular,
                                                        color = color,
                                                        ...),
                                       homogeneous = FALSE))
    
    if(add_normals) mesh <- mesh %>% rgl::addNormals()
    
    class(mesh) <- c("obj_mesh", class(mesh))
    
    mesh
  }

#' Colourize a mesh
#' 
#' Add colour information to your mesh, either from a vertex atlas like AAL
#' or from a statistic/measurement map like those produced by CIVET
#' 
#' @param mesh \link[rgl]{mesh3d} object ideally produced by \link{create_mesh}
#' @param colour_map either a vector with a label/measure/statistic for every vertex
#' or a character file path pointing to a file with such a vector in a rowwise format.
#' @param colour_range a two element numeric vector indicating the min and max values of 
#' allowable labels/measures/statistics to be includedon the surface
#' @param colour_default The colour given to vertices excluded by colour_range
#' @param reverse Whether to have a positive and negative colour scale (not yet implemented)
#' @param palette A palette, AKA look-up-table, providing a linear colour scale for the colours in
#' \code{colour_map}
#' @return an \code{obj_mesh} object descended from \link[rgl]{mesh3d}, with added colour information
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
      (colour_map - colour_range[1]) / 
        diff(colour_range) * (colour_depth - 1)) + 1
  
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
#' Create a basic plot of a BIC obj in the current rgl device, opening a 
#' new one if necessary. If colour_map is supplied, an overlay is added to the mesh
#' 
#' @param x A \code{bic_obj} probably created by \link{read_obj}
#' @param colour_map A numeric vector equal in length to the number of vertices
#' in the \code{bic_obj} or the path to a text file with one line per vertex with
#' colour information.
#' @param colour_range a two element numeric vector indicating the min and max values of 
#' allowable labels/measures/statistics to be includedon the surface
#' @param colour_default The colour given to vertices excluded by colour_range
#' @param reverse Whether to have a positive and negative colour scale (not yet implemented)
#' @param palette A palette, AKA look-up-table, providing a linear colour scale for the colours in
#' \code{colour_map}
#' @param colour_bar whether to draw a colour bar
#' @param add whether or not to add this object to the current rgl device (if possible)
#' defaults to opening a new device
#' @param ... additional arguments to \link{create_mesh} including but not limited
#' to colour, specular, and add_normals
#' @return invisibly returns the mesh object
#' @method plot bic_obj
#' @export
plot.bic_obj <- 
  function(x, 
           colour_map = NULL, 
           colour_range = NULL,
           colour_default = "grey",
           reverse = NULL,
           palette = heat.colors(255),
           colour_bar = TRUE,
           add = FALSE,
           ...){
    
    if(!add) rgl::open3d()
    
    rgl::.check3d()
    window_dimensions <- rgl::par3d("windowRect")
    
    #Small windows prevent legend plotting
    if(diff(window_dimensions[c(1,3)]) < 400 && 
       diff(window_dimensions[c(2,4)]) < 400 ){
      rgl::par3d(windowRect = c(100,100, 900, 900))
      #This isn't my fault I swear, see the ?bgplot3d examples, weird bugs happen
      #without waiting, e.g. brains will only draw when there is an inactive rgl device
      Sys.sleep(.25)
      rgl::par3d(viewport = c(0,0,800,800)) 
    }
    
    mesh <-
      x %>%
      create_mesh(...)
    
    if(!is.null(colour_map))
      mesh <- mesh %>% colour_mesh(colour_map, 
                                   colour_range, 
                                   colour_default, 
                                   reverse, 
                                   palette)
    
    mesh %>% rgl::shade3d(override = FALSE)
    if(colour_bar && !is.null(colour_map)) mesh %>% add_colour_bar
    
    invisible(mesh)
  }

#' Plot an BIC obj mesh
#' 
#' Create a plot of BIC obj_mesh, potentially colourized by \link{colour_mesh}
#' and potentially including a colour bar
#' 
#' @param x a \code{obj_mesh} object
#' @param colour_bar whether or not to add a colour bar
#' @param ... additional parameters to pass to add_colour_bar
#' @return returns x invisibly
#' @method plot obj_mesh
#' @export
plot.obj_mesh <-
  function(x, colour_bar = TRUE, ...){
    x %>% rgl::shade3d
    if(colour_bar) x %>% add_colour_bar(...)
    
    invisible(x)
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
  if(is.null(mesh$legend)) stop("Your mesh has no colour information")
  with(mesh$legend, {
    rgl::bgplot3d({
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
#' @param colour_title legend title for the colour bar if requested
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
  
  rgl::open3d(windowRect = c(100,100, 900, 900))
  #This isn't my fault I swear, see the ?bgplot3d examples, weird bugs happen
  #without waiting, e.g. brains will only draw when there is an inactive rgl device
  Sys.sleep(.25)
  rgl::par3d(viewport = c(0,0,800,800))

  parent_scene <- rgl::currentSubscene3d()
  
  subscenes <- rgl::mfrow3d(3,2)
  left_or_right <- rep(c("left", "right"), 3)
  
  mapply(function(subscene, view_matrix, hemisphere){
    rgl::useSubscene3d(subscene)
    if(hemisphere == "left"){
      left_mesh %>% rgl::shade3d(override = FALSE)
    } else {
      right_mesh %>% rgl::shade3d(override = FALSE)
    }
    
    rgl::par3d(userMatrix = view_matrix)
  }, subscenes, hemisphere_viewpoints, left_or_right)

  rgl::useSubscene3d(parent_scene)
  
  left_mesh %>% add_colour_bar(title = colour_title, 
                               column_widths = c(6,2,6), 
                               which_col = 2, 
                               cex.axis = 2, 
                               cex.lab = 2)
  
  if(!is.null(output)) rgl::snapshot3d(output)
  if(!is.null(output) && close_on_output) rgl::rgl.close()
  
  invisible(NULL)
}

#' Find Closest Vertex
#' 
#' Given a set of vertices (3-column matrix or data frame of x,y,z coordinates)
#' determine the index or coordinates of the closest vertex to a given set of
#' target vertices (either a 3 element vector, or x-y-z table like structure)
#' 
#' @param vertices A matrix-like object with 3-columns, an n rows representing vertices
#' @param target either a 3-element numeric vector representing x-y-z coordinates for
#' a single target, or a matrix-like object as described above containing multiple targets.
#' @param returns Whether to return the index of each match (one per target), or the coordinates
#' of the matches, the later being useful when exact matches aren't expected.
#' @return In most cases, a numeric vector of match results, in the case of multiple targets
#' and \code{returns = "coordinates"} the output is a 3-column matrix of coordinates.
#' 
#' @details Vertices are matched to the target by finding the vertex with the minimum euclidean
#' distance from the target.
#' @export
closestVertex <-
  function(vertices, target, returns = c("index", "coordinates")){
    if(!is.data.frame(target) & !is.matrix(target))
      target <- matrix(target, ncol = 3)
    
    return_type <- match.arg(returns)
    
    results <- apply(target, 1, single_closest_vertex, vertices = vertices, returns = return_type)
    
    results <- t(results)
    if(nrow(results) == 1) dim(results) <- NULL
    
    results
  }

single_closest_vertex <-
  function(vertices, target, returns){
    if(ncol(vertices) == 3) vertices <- cbind(vertices, 1)
    if(length(target) == 3) target <- c(target, 1)
    
    distances <- colSums((t(vertices) - as.numeric(target))^2)
    closest <- which(distances == min(distances))
    
    if(returns == "coordinates") return(vertices[closest, -4])
    
    return(closest)
  }

#' Select Vertices From a Surface
#' 
#' Interactively select one or more vertices from an \code{rgl} surface.
#' 
#' @param object An rgl id for the object from which to choose vertices. Defaults to
#' the first id. To check the ids available for selection call \link[rgl]{rgl.ids}()
#' unfortunately there is no elegant way at the present for determining which object is which.
#' @param tolerance The tolerance, in fractions of the visible window, for finding points.
#' See details for more.
#' @param multiples Whether to stop after the first selection, or to allow multiples (explicitly stopped
#' by pressing escape). 
#' @param indicate Whether to draw indicator points to highlight the clicked points. Useful for ensuring
#' point click accuracy.
#' 
#' @details 
#' The vertex selection algorithm has two parts, selecting the vertices near the clicked point 
#' and determining which points are closest to the observer. \cr
#' To determine the vertices near the click point
#' \itemize{
#' \item{Determine the x-y coordinates of the click in window space}
#' \item{Convert those coordinates to a rectangle defined by 
#' (x - tolerance, y - tolerance, x + tolerance, y + tolerance)}
#' \item{Convert the rectangle coordinates from window to user coordinates}
#' \item{Filter the vertices of the object within the rectangle}
#' }
#' To determine the vertex closest to observer
#' \itemize{
#' \item{Determine the x-y click location}
#' \item{Append z=0 to the coordinates to create a 3d window coordinate, since the
#' rgl observe looks in the positive z direction, a window z coordinate of zero puts the point
#' as far from the surface as possible while remaining in the window coordinates}
#' \item{Convert from window to user coordinates}
#' \item{Calculate the vertex from the set identified above that is closest to those coordinates}
#' }
#' This algorithm is not perfect, and can yeild spurious coordinates for irregular topologies.
#' For more accurate vertex selection, use a high magnification (controlled with the scroll wheel), 
#' the higher the magnification the more accurate the vertex selection becomes. Additionally, keep
#' indicate = TRUE, this will place indicator points on the identified vertices, they will allow you
#' to ensure your coordinates are accurate, and can always be removed with \link[rgl]{pop3d}()
#' @export
vertexSelect <-
  function(object = first(rgl::rgl.ids()$id),
           tolerance = 0.01,
           multiples = FALSE,
           indicate = TRUE){
    
    first = TRUE
    selected_vertices <- matrix(NA, nrow = 0, ncol = 3)
    
    while(first || multiples){
      first <- FALSE
      selected_vertex <- select_one_vertex(object, tolerance, indicate)
      if(is.null(selected_vertex)) break
      
      selected_vertices <- rbind(selected_vertices, selected_vertex)
    }
    rownames(selected_vertices) <- NULL
    colnames(selected_vertices) <- c("x", "y", "z")
    
    if(indicate){
      bbox <- rgl::par3d("bbox")
      bbox_smallest_dimension <- min(abs(bbox[4:6] - bbox[1:3]))
      
      rgl::spheres3d(selected_vertices[,1], 
                selected_vertices[,2], 
                selected_vertices[,3], 
                alpha = .2,
                radius = bbox_smallest_dimension / 100,
                specular = "black")
    }
    
    if(nrow(selected_vertices) == 1){
      dim(selected_vertices) <- NULL
      names(selected_vertices) <- c("x", "y", "z")
    }
    
    return(selected_vertices)
  }

select_one_vertex <- 
  function(object, tolerance, indicate, ...){
    
    selection <- as.list(environment(rgl::select3d(...)))
    if(is.null(selection$rect)) return(invisible(NULL))
    
    rect_bounds <- selection$rect
    click_location <- rect_bounds[1:2]
    
    if(rect_bounds[1] != rect_bounds[3] || rect_bounds[2] != rect_bounds[4])
      warning("Click-and-drag selection does not work with vertexSelect, using position of button press",
              call. = FALSE)
    
    inflated_bounds <- rep(click_location, 2) - c(1,1,-1,-1) * tolerance
    
    object_vertices <- rgl::rgl.attrib(object, "vertices")
    window_coords <- rgl::rgl.user2window(object_vertices, projection = selection$proj)
    
    selected_vertices <-
      which(window_coords[,1] > inflated_bounds[1] &
              window_coords[,2] > inflated_bounds[2] &
              window_coords[,1] < inflated_bounds[3] &
              window_coords[,2] < inflated_bounds[4])
    
    if(length(selected_vertices) == 0) return(NULL)
    
    selected_coords <-
      object_vertices[selected_vertices,]
    
    unique_coords <-
      selected_coords[!duplicated(selected_coords, MARGIN = 1),]
    
    target_location <- rgl::rgl.window2user(click_location[1], click_location[2], 0, 
                                       projection = selection$proj)
    
    closest_vertex <- closestVertex(unique_coords, target_location, returns = "coordinates")
    
    closest_vertex
  }

#' Vertex Lookup
#' 
#' Find the vertices closest one or more targets, potentially returning the values
#' for the vertices from a data map.
#' 
#' @param vertices A descendent of \link[rgl]{mesh3d}, \code{bic_obj}, or matrix-like object with 3-columns, 
#' and n rows representing vertices.
#' @param target either a 3-element numeric vector representing x-y-z coordinates for
#' a single target, or a matrix-like object as described above containing multiple targets.
#' @param data_map Either NULL, a vector of data about each vertex, or a file containing such
#' a vector spread over multiple lines
#' @param returns Whether to return the index of each match (one per target), or the coordinates
#' of the matches, the later being useful when exact matches aren't expected.
#' @param coerce A function to coerce the final results to a given type. Defaults to \link{as.numeric},
#' if set to NULL, no coersion is performed.
#' @return If a data_map is specified: a vector, typically numeric, if coerce is set to NULL 
#' and data_map is a file, the results will be character. If coerce is null and data_map is a vector
#' it will return the same type as data_map. If data_map is unspecified, it acts like \link{closestVertex}
#' @export
vertexLookup <- 
  function(vertices, target, data_map = NULL,
           returns = c("index", "coordinates"),
           coerce = as.numeric){
    if(vertices %>% is("mesh3d")) vertices <- t(vertices$vb)
    if(vertices %>% is("bic_obj")) vertices <- t(vertices$vertex_matrix)
    
    if(is.null(data_map)){
      return_type <- match.arg(returns)
      closest <- closestVertex(vertices, target, returns = return_type)
      return(closest)
    }
    
    closest <- closestVertex(vertices, target)
    if(is.character(data_map) && length(data_map) == 1)
      data_map <- readLines(data_map)
    
    map_value <- data_map[closest]
    
    if(is.null(coerce)) return(map_value)
    
    coersion_function <- match.fun(coerce)
    coersion_function(map_value)
  }

#' Read a BIC-obj line file
#' 
#' Parse the BIC obj format for when the object contains
#' lines instead of a mesh. 
#' 
#' @param lines_obj Path to the object file of interest
#' @return \code{bic_lines} object, which is a list of matrices, each 
#' matrix coresponds to one line in the object. The matrices are 3xN matrices 
#' of world coordinates.
#' @export
read_line_obj <-
  function(line_obj){
    lines <- readLines(line_obj)
    
    general_info <- lines[1]
    
    ## Hacky parsing BIC .obj format see https://github.com/BIC-MNI/bicpl/Documentation
    ## This probably only works with output from CIVET 1.1.12
    section_ends <- 
      lines %>%
      `==`("") %>%
      which %>%
      as.list %>%
      setNames(c("vertices", "colours", "line_ends"))
    
    parse_numbers <- 
      function(lines){
        lines %>%
          strsplit(" ") %>%
          unlist %>%
          Filter(f = function(chr) chr != "" ) %>%
          as.numeric
      }
    
    vertices <- 
      lines[2:section_ends$vertices] %>%
      parse_numbers %>%
      matrix(nrow = 3)
    
    line_ends <-
      lines[(section_ends$colours + 1):(section_ends$line_ends - 1)] %>%
      parse_numbers
    
    line_frame <-
      data_frame(end = line_ends, start = lag(line_ends, default = 0))
    
    lines_list <- 
      mapply(
        function(start, end){
          vertices[,(start+1):end]
        }, 
        start = line_frame$start, end = line_frame$end,
        SIMPLIFY = FALSE
      )
    
    class(lines_list) <- c("bic_lines", "list")
    attr(lines_list, "coord_system") <- "world" 
    
    return(lines_list)
  }

#' Convert Lines to Voxel Coordinates
#' 
#' Convert a \code{bic_lines} object to world coordinates for plotting
#' 
#' @param line_obj The \code{bic_lines} object of interest
#' @param minc_file The reference file for computing voxel coordinates
#' @return A \code{bic_lines} object in voxel coordinates
#' @export
line_obj_to_voxel <-
  function(line_obj,
           minc_file){
    
    stopifnot(inherits(line_obj, "bic_lines"), attr(line_obj, "coord_system") == "world")
    
    lines_list <-
      lapply(line_obj, function(line){
        mincConvertWorldMatrix(minc_file, line, nearest_voxel = FALSE)
      })
    
    class(lines_list) <- c("bic_lines", "list")
    attr(lines_list, "coord_system") <- "voxel" 
    
    return(lines_list)
  }

#' Plot A bic_lines object
#' 
#' Add lines corresponding to the coordinates in a bic_lines
#' object to a figure. Only draws the projection of the lines
#' on a single dimension, no regard is given for whether the
#' lines are near the slice of interest.
#' 
#' @param x an \code{bic_lines} object
#' @param dimension which axis to display the lines on
#' @param ... additional parameters to pass to \link{segments}
#' @return NULL invisibly
#' @export
plot.bic_lines <-
  function(x, dimension = 2, ...){
    stopifnot(inherits(x, "bic_lines"))
    
    lapply(x,
           function(line){
             line <- line[-dimension,]
             line_frame <- 
               line %>%
               t %>%
               as.data.frame %>%
               setNames(c("y0", "x0")) %>% #X and Y are transposed
               mutate(                     #in mincImage so follow suit
                 x1 = lag(x0),
                 y1 = lag(y0)
               ) %>%
               with(segments(x0, y0, x1, y1, ...))
             
             NULL
           })
    
    return(invisible(NULL))
  }






