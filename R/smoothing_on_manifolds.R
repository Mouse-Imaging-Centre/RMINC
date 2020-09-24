#___ Functions to clean up manifold ___

# find neighbouring vertices
find_neigh_verts = function(idx,manifold,clean=F) {
        triangle_matrix = manifold$triangle_matrix
        ret = as.vector(triangle_matrix[,apply(triangle_matrix, 2, function(i) any(i==idx))])
        ret = unique(ret)
        if (clean) {
                ret = sort(ret)
                ret[ret != idx]
        }
        ret
}

# find neighbouring vertices with less than n-degrees of seperation
find_n_neigh_verts = function(idx,n,manifold,clean=F) {
        ret = find_neigh_verts(idx,manifold, clean=F)
        while (TRUE) {
                n = n - 1
                if ( n == 0 ) {break}
                ret = unique(unlist(lapply(ret, function(i) find_neigh_verts(i,manifold, clean=F))))
        }
        if (clean) {
                ret = sort(ret)
                ret[ret != idx]
        }
        ret
}

# find number of neighbouring vertices with less than n-degrees of seperation
find_num_n_neigh_vert = function(idx,n,manifold) {
        length(find_n_neigh_verts(idx,n,manifold))
}

# find all edges in manifold triangle mesh 
#  If symm is true, it will repeat symmetric edges 
#   i.e: edge a-b and edge b-a are counted
find_all_edges = function(manifold,symm=T) {
        triangle_matrix = manifold$triangle_matrix
        cot_matrix = manifold$cot_matrix
        edges = do.call('rbind',lapply(1:3, function(z) {
                       edges = t(triangle_matrix[-z,])
                       colnames(edges) = c('i','j')
                       edges
        }))
        if (symm) {
                edges = rbind(edges,edges[,c(2,1)])
        }
        edges
}

# If an edge has more than 2 attached triangle,
#  Remove all but 2 triangles
#  Pick triangles with vertices that have the most neighbours (at most k degrees-of-seperation)
filter_problematic_triangles = function(edge,manifold,k=1) {
        triangle_matrix = manifold$triangle_matrix
        bool1 = apply(triangle_matrix == edge[1],2,any)
        bool2 = apply(triangle_matrix == edge[2],2,any)
        bool = bool1 & bool2

        if ( sum(bool1 & bool2) < 3 ) {stop('critical error')}

        sub_tri_mat = manifold$triangle_matrix[ , bool ]
        test_verts = apply(sub_tri_mat, 2, function(x) x[!x %in% edge] )

        # find number of 1st degree neighbors for each ambiguous vertex.
        # keep triangles with the ambiguous vertex with the two highest 1st degree neighbors
        num_neigh = unlist(lapply(test_verts, function(i) find_num_n_neigh_vert(i,k,manifold)))
        which_verts = order(num_neigh,decreasing=T)[1:2]

        idxs = which(bool)
        idxs[-which_verts]
}

#' Clean up manifold mesh
#'
#' If mesh edge has 2+ attached triangles, keep the 2 triangles most well-connected triangles
#'
#' @param manifold A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
#' like \code{bic_obj} object produced by \link{read_obj}
#' @param k degrees-of-seperation for defining neighbours. Default is 1 (i.e. adjacent vertices)
#' @return A list of length 2 with names 'vertex_matrix' and 'triangle_matrix' 
#' similar to 'manifold' argument. Problematic triangles (if they exist) are removed from the mesh.
#' @details increasing k drastically slows down computation
#' @export
clean_up_manifold_triangles = function(manifold,k=1) {
   # clean up triangle mesh
   # for object to be a proper mesh, all edges must not have more than 2 attached triangles
   # thus, each symmetric edge must not be counter more than twice
   edges = find_all_edges(manifold)
   edges = edges[duplicated(edges),]
   problematic_edges = unique(edges[duplicated(edges),])

   if (length(problematic_edges) != 0) {
      cat ('\n',
           'There are some edges in the manifold mesh that are attached to','\n',
           'more than two triangles. These will be removed for computing','\n',
           'the laplace-beltrami operator for this manifold',
           '\n')
      problematic_triangle_idxs = unlist(lapply(
             1:nrow(problematic_edges) ,
             function(i) filter_problematic_triangles(problematic_edges[i,],manifold,k=k)))

      manifold$triangle_matrix = manifold$triangle_matrix[,-problematic_triangle_idxs]
   }
   return(manifold)
   # TODO: this function is a bit of a bottleneck. It can probably be made faster.
}

#___ Helper functions for computing the discrete laplace-beltrami operator ___
# given 2 points, compute the distance-squared between them
position_to_edgelengths_sq = function(posmat) {
        sum((posmat[,1] - posmat[,2])^2)
}

# given 3 points defining a triangle, compute all pairwise distance-squared
triangle_to_edgelengths_sq = function(vertex_matrix) {
        sapply(1:3, function(z) {
                       position_to_edgelengths_sq(vertex_matrix[,-z])
         })
}

# for a manifold defined by a triangle mesh, compute all pairwise distance-squared 
#  for all triangles
find_edgelengths_sq = function(manifold) {
        vertex_matrix   = manifold$vertex_matrix
        triangle_matrix = manifold$triangle_matrix
        apply(triangle_matrix,2,function(i) triangle_to_edgelengths_sq(vertex_matrix[,i]))
}

# Given distance-squared for the sides of triangle, compute cotangent of all interior angles
#  Note: each angle is paired with the edge opposite to it
edgelen_sq_to_cot_angles = function(edgelen_sq) {
        sapply(1:3, function(z) {
                       x = 0.5*(sum(edgelen_sq[-z]) - edgelen_sq[z])/sqrt(prod(edgelen_sq[-z]))
                       x/sqrt(1-x^2)
         })
}

# for a manifold defined by a triangle mesh, compute the cotangent of all angles
#  Note: you have to run \code{find_edgelengths_sq} function first
find_cot_angles = function(manifold) {
        vertex_matrix   = manifold$vertex_matrix
        triangle_matrix = manifold$triangle_matrix
        edgelen_sq = manifold$edgelen_sq
        apply(edgelen_sq,2,edgelen_sq_to_cot_angles)
}

# Given two vectors connected tail-to-tail defining 2 edges of a triangle, 
#  compute the area of the triangle
cross_product_triangle_area = function(u,v) {
        xprodvec = c(
          u[2]*v[3] - u[3]*v[2],
          u[3]*v[1] - u[1]*v[3],
          u[1]*v[2] - u[2]*v[1]
          )
        0.5 * sqrt(sum(xprodvec^2))
}

# for a manifold defined by a triangle mesh, compute the area of triangles
find_triangle_areas = function(manifold) {
        triangle_matrix = manifold$triangle_matrix
        vertex_matrix   = manifold$vertex_matrix

        traingle_area = apply(triangle_matrix,2,function(x) {
                               y = vertex_matrix[,x]
                               cross_product_triangle_area( y[,2] - y[,1] , y[,3] - y[,1] )
          })
        traingle_area
}

# for a manifold defined by a triangle mesh, compute the vertex area of all vertices
#  Note: for any vertex, the vertex area is the sum of the areas of the 
#   triangles connected to it, divided by 3.
find_vertex_areas = function(manifold) {
        triangle_matrix = manifold$triangle_matrix
        traingle_area = manifold$traingle_area
        vertex_matrix = manifold$vertex_matrix

        vertex_area = rep(NA,ncol(vertex_matrix))

        tmpsum = aggregate(
                  rep(traingle_area, each=3),
                  list(as.vector(triangle_matrix)),
                  sum)
        colnames(tmpsum) = c('vertex','area')
        tmpsum[,'area'] = tmpsum[,'area']/3
        vertex_area[tmpsum[,'vertex']] = tmpsum[,'area']

        vertex_area
}

#___ Functions for computing the discrete laplace-beltrami operator ___
# The discrete laplace-beltrami operator for a mesh defined my N vertices is
#  a sparse N-by-N matrix.
# If points in the mesh are indexed by i and j, 
#  the value of a non-diagonal element in the matrix is:
#    0                if i and j are not connected
#    0.5*(ca + cb)/Ai if i and j are connected
#                     where Ai is the vertex area of point i
#                           ca and cb are the cotangent of the two angles opposite edge ij
# To find the value of the diagonal elements, 
#  take the negative sum of of all non-diagonal elements in that row

# Function to compute the nonzero elements of the laplace-beltrami operator
find_laplace_beltrami_nonzero_elements = function(manifold) {

        triangle_matrix = manifold$triangle_matrix
        cot_matrix = manifold$cot_matrix
        traingle_area = manifold$traingle_area
        vertex_area = manifold$vertex_area

        edges = do.call('rbind',lapply(1:3, function(z) {
                       edges = t(triangle_matrix[-z,])
                       colnames(edges) = c('i','j')
                       edges
        }))
        edges = rbind(edges,edges[,c(2,1)])

        cots = unlist(lapply(1:3, function(z) {
                       cot_matrix[z,]
        }))
        cots = c(cots,cots)

        tmpordr = order(edges[,'i'],edges[,'j'])
        edges = edges[tmpordr,]
        cots  = cots[tmpordr]

        # check if there are any edges attached to more than 2 trianges
        tst = edges %>%
                (function(x) {x[duplicated(x),]}) %>%
                (function(x) {x[duplicated(x),]}) %>%
                nrow
        if (tst != 0) {stop('some problematic edges remain')}

        # find edges with two attached triangles
        bool1 = duplicated(edges)
        bool2 = rev(duplicated(edges[rev(1:nrow(edges)),]))

        if ( any(apply(edges[bool1,] != edges[bool2,],1,any)) ) {
                stop('critical error')
        }

        # find edges with one triangle attached
        single_bool = ((!bool1) & (!bool2))

        # for these triangles, use reflective boundary condition
        # (caution: area scaling may be off but shouldn't change results much)
        ret_single_edges = edges[single_bool, ]
        ret_single_cot   = cots[single_bool]

        # for edges with two attached triangles, take mean of cot
        ret_double_edges = edges[bool1,]
        ret_double_cot   = 0.5*(cots[bool1] + cots[bool2])

        # bind off-diagonal elements
        ret_odiag_edges = rbind(ret_double_edges,ret_single_edges)
        ret_odiag_cot   = c(ret_double_cot,ret_single_cot)

        # symmetric elements are the negative sum of off-diagonals
        tmpres = aggregate(ret_odiag_cot,by = list(ret_odiag_edges[,1]), FUN = sum)
        ret_diag_edges = cbind(tmpres[,1],tmpres[,1]) ; colnames(ret_diag_edges) = c('i','j')
        ret_diag_cot   = -tmpres[,2]

        # bind all elements
        ret_edges = rbind(ret_odiag_edges,ret_diag_edges)
        ret_cot   = c(ret_odiag_cot,ret_diag_cot)

        # scale elements with inverse vertex area
        ret_cot = ret_cot/vertex_area[ret_edges[,'i']]

        list(edges = ret_edges, elements = ret_cot)
}

#' Compute Laplace-Beltrami operator
#'
#' Discrete Laplace-Beltrami operator associated with a triangle-mesh manifold
#'
#' @param manifold A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
#' like \code{bic_obj} object produced by \link{read_obj}
#' @param vertex_matrix 3-by-N matrix denoting the position of the N vertices 
#' defining the triangle-mesh manifold. Similar to argument \code{vertices} in
#' \link[rgl:mesh3d]{rgl::tmesh3d}. Required if \code{manifold} argument is not specified.
#' @param triangle_matrix 3-by-M matrix denoting the position of the M triangles
#' defining the triangle-mesh manifold. Elements of each column are indices of vertices 
#' defining the triangle (i.e. indices are columns of the vertex_matrix). Similar to
#' argument \code{indices} in \link[rgl:mesh3d]{rgl::tmesh3d}. 
#' Required if \code{manifold} argument is not specified.
#' @return sparse square matrix representing the discrete Laplace-Beltrami operator
#' for the manifold
#' @details Either supply \code{manifold} argument OR \code{vertex_matrix} and \code{triangle_matrix}
#' arguments. 
#' @export
laplace_beltrami_operator = function(
                                     manifold = NULL ,
                                     vertex_matrix = NULL ,
                                     triangle_matrix = NULL ) {
        defense_func = function(i) {
                c1 = function() {
                        cat('\n',
                            'You must either supply the manifold argument','\n',
                            'OR BOTH vertex_matrix and triangle_matrix arguments',
                        '\n')
                }
                c2 = function() {
                        cat('\n',
                           'manifold argument must be a list with names','\n',
                           '"vertex_matrix" and "triangle_matrix"',
                        '\n')
                }
                c3 = function() {
                        cat('\n',
                   'Let N be number of vertices and M be number of triangles','\n',
                   'vertex_matrix:  denoting the position of every vertex (3 by N).','\n',
                   '   Equivilant to "vertices" argument in the rgl package function','\n',
                   'triangle_matrix: denoting the triangle mesh connecting the vetrices (3 by M)','\n',
                   '   Equivilant to "indices" argument in the rgl package function',
                        '\n')
                }
                if (i == 1) {c1();c2();c3()}
                if (i == 2) {c2();c3()}
                if (i == 3) {c3()}
        }

        if (is.null(manifold) & (is.null(vertex_matrix) | is.null(triangle_matrix))) {
                defense_func(1); cat('\n')
                stop('Error in input arguments')
        }
        if (!is.null(manifold)) {
                if (!all(c('vertex_matrix','triangle_matrix') %in% names(manifold))) {
                        defense_func(2); cat('\n')
                        stop('Error in input arguments')
                }
        } else {
                manifold = list(
                                vertex_matrix = vertex_matrix,
                                triangle_matrix = triangle_matrix
                                )
        }
        if (is.null(dim(manifold$vertex_matrix)) | is.null(dim(manifold$triangle_matrix)) ) {
                defense_func(3); cat('\n')
                stop('Error in input arguments')
        }
        if ( (nrow(manifold$vertex_matrix) != 3) | (nrow(manifold$triangle_matrix) != 3) ) {
                defense_func(3); cat('\n')
                stop('Error in input arguments')
        }

        # check if there are any edges connected to more than 2 triangles and remove them
        manifold = clean_up_manifold_triangles(manifold)

        # calculate length of all edges (squared). 
        # Result (edgelen_sq) has the same shape as triangle_matrix 
        #   In triangle_matrix, each element is a vertex
        #   In edgelen_sq     , each element is the edge length square opposite the vertex
        manifold$edgelen_sq = find_edgelengths_sq(manifold)

        # calculate cotangent of all angles in the mesh
        # Result (cot_matrix) has the same shape as triangle_matrix
        #   In triangle_matrix, each element is a vertex
        #   In cot_matrix     , each element of the cotangent of the vertex angle
        manifold$cot_matrix = find_cot_angles(manifold)

        # Find area of all triangles
        manifold$traingle_area = find_triangle_areas(manifold)

        # Find vertex area of all vertices
        #  The vertex area is the sum of triangle areas attached to the vertex divided by 3
        manifold$vertex_area = find_vertex_areas(manifold)

        # calculate the indices and values of the non-zero elements of the 
        #  laplace-beltrami operator
        manifold$laplace_beltrami_elements = find_laplace_beltrami_nonzero_elements(manifold)

        n = ncol(manifold$vertex_matrix)
        # define laplace_beltrami as sparse matrix
        #  i and j index the vertices
        #  element at i and j are non-zero if: 
        #   1) i and j are connected by an edge
        #   2) i=j
        laplace_beltrami_operator =
                sparseMatrix(
                    i=manifold$laplace_beltrami_elements$edges[,1],
                    j=manifold$laplace_beltrami_elements$edges[,2],
                    x=manifold$laplace_beltrami_elements$elements,dims=list(n,n))

        return(laplace_beltrami_operator)
}

#___Functions for smoothing simulations___

# Function to compute optimum time-step-size for the simulation
#    If step size is too small, convergence is slow
#    If step size is too large, 
#       variance of scalar field values increases 
#       the field has values outside the range of values in the initial conditions
#    This is nonsensical for the heat equation because
#      scalar field should become MORE homogeneous with time (reduced variance)
#      and should never have values outside the minimum and maximum values of the initial conditions
# Select the optimum time-step for the simulation by 
#    minimizing the variance at each iteration
#    penalizing time-steps where the scalar field values are outside the 
#     range of values in the initial conditions
# Note: time is returned as log-10
tstep_size_select = function(cmat, u, mnu = NULL , mxu = NULL) {
        n = ncol(cmat)
        if (is.null(mnu)) {mnu = min(u)-0.001}
        if (is.null(mxu)) {mxu = max(u)+0.001}
        ru = diff(range(u))
        tmpfunc = function(x) {
                y = ((cmat*(10^x) + Diagonal(n)) %*% u)[,1]
                v1 = var(y)
                v2 = 10*as.numeric(max(y)>mxu)+1
                v3 = 10*as.numeric(min(y)<mnu)+1
                v1 * v2 * v3
        }
        optimize(tmpfunc, c(-10,10))$minimum
}

# There is an easy relationship between simulation time and FWHM if:
#   diffusion constant is 0.5 (easy assumption because it just requires rescaling time units) 
#    AND
#   space is euclidean (bad assumption)
# Nevertheless, this relationship helps put 
#  simulation time in intuitive units of FWHM (if space is flat)
tsim_to_fwhm = function(tsim) {
        sqrt(8 * log(2) * tsim)
}
fwhm_to_tsim = function(fwhm) {
        ((fwhm)^2)/(8*log(2))
}

#' Smoothing on manifold
#'
#' Smoothing scalar field on a triangle-mesh manifold
#'
#' @param manifold A list of length 2 with names 'vertex_matrix' and 'triangle_matrix'
#' like \code{bic_obj} object produced by \link{read_obj}
#' @param scalar_field a vector whose elements are values of the field at each manifold vertex
#' @param fwhm the degree of smoothing. It represents the full-width-half-maximum 
#' of the blurring kernel if the manifold had no curvature. Regardless of curvature,
#' higher the fwhm, the greater the amount of smoothing.
#' @param maxiter int specifying the maximum number of iterations for the algorithm
#' @return a vector whose elements are values of the smoothed field at each manifold vertex
#' @examples
#' \dontrun{
#' # Load an object
#' manifold =
#'   read_obj(
#'     file.path("/axiom2/projects/software/cortical-thickness/",
#'               "MWM/c57bl6_laplacian_grid_full_surface_simplified.obj"))
#'
#' # Compute the laplace beltrami operator and attach it to the manifold
#' #  not necessary but will make future smoothing computations on the same manifold faster
#' manifold$laplace_beltrami_operator = laplace_beltrami_operator(manifold)
#'
#' # Generate some random uniform vertex data
#' init_stats = runif(ncol(manifold$vertex_matrix))
#'
#' # Smooth the stats on the manifold
#' smooth_stats = laplace_beltrami_smoothing(manifold,init_stats,0.2)
#' 
#' # Plot results
#' plot(manifold, init_stats, colour_range = c(.5,1))
#' plot(manifold, smooth_stats, colour_range = c(.5,1))
#'}
#' @export
laplace_beltrami_smoothing = function( manifold, scalar_field, fwhm , maxiter = 1000) {
        if (!'laplace_beltrami_operator' %in% names(manifold)) {
                lbmat = laplace_beltrami_operator(manifold)
        } else {
                lbmat = manifold$laplace_beltrami_operator
        }

        # number of vertices
        n = ncol( lbmat )

        if (n != length(scalar_field)) {
           stop(
                'number of vertices [',
                n,
                '] should equal number of elements in scalar_field [',
                length(scalar_field),
                ']'
                )
        }
        # convert fwhm into simulation time
        total_sim_time = fwhm_to_tsim(fwhm)

        # define some values to help optimize stepsize for convergence
        min_scal_val = min(scalar_field)-1E-10
        max_scal_val = max(scalar_field)+1E-10

        t_current = 0 ; brk = F
        for (i in 1:maxiter) {
                if ((i %% 10) == 0) {
                 cat(
                     'Iter:',i,
                     'Current FWHM:',tsim_to_fwhm(t_current),
                     'Target FWHM:',fwhm,
                     '\n')
                }
                tstep_log10 = tstep_size_select(lbmat,scalar_field, min_scal_val, max_scal_val)
                tstep = 10^tstep_log10
                t_remaining = (total_sim_time - t_current )

                if ( t_remaining < tstep ) {
                        tstep = t_remaining
                        brk = T
                }

                t_operator = lbmat*tstep + Diagonal(n)
                scalar_field = (t_operator %*% scalar_field)[,1]
                t_current = tstep + t_current

                if (brk) {
                        cat('simulation completed in',i,'iterations','\n')
                        break
                }
        }
        if (!brk) {
                cat('simulation DID NOT complete in',i,'iterations','\n')
                cat('Max iterations reached','\n')
        }
        return(scalar_field)
}


