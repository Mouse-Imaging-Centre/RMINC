create_parallel_mask <-
  function(sample_file, mask = NULL, n, temp_dir = getwd(), prefix = "pMincMask", tinyMask = FALSE){
    
    sample_volume <- mincGetVolume(sample_file)
    mask_file <- "unset"
    
    temp_dir <- path.expand(temp_dir)
    if(!file.exists(temp_dir)) dir.create(temp_dir)
    
    mask_file <- tempfile(prefix, tmpdir = temp_dir, fileext = ".mnc")
    
    if(is.null(mask)){
      if(n %% 1 != 0) stop("the number of jobs must be an integer")
      
      nVoxels <- length(sample_volume)
      mask_values <- groupingVector(nVoxels, n)
    } else {
      mask_values <- mincGetVolume(mask)
      if(tinyMask) mask_values[mask_values > 1.5] <- 0
      
      nVoxels <- sum(mask_values > .5)
      
      mask_values[mask_values > .5] <- 
        groupingVector(nVoxels, n)
      
    }
    
    mincWriteVolume(mask_values, mask_file, like.filename = sample_file)
    
    return(mask_file)
  }

tenacious_remove_registry <- 
  function(reg){
    tryCatch(removeRegistry(reg = reg)
             , error = 
               function(e){
                 killJobs(reg = reg, findNotDone(reg = reg))
                 removeRegistry(reg = reg)
               })
  }

new_file <- function(name, dir = getwd()){
  files <- grep(name, list.files(dir), value = TRUE)

  i <- 1
  while(sprintf("%s%03d", name, i) %in% files)
    i <- i + 1

  sprintf("%s%03d", name, i)
}

#' Parallel MincApply
#' 
#' Apply an arbitrary R function across a collection of minc files, distributing
#' the computation to multiple cores or workers on a grid computing environment
#' 
#' @param filenames Paths to the minc files to be applied accross
#' @param fun The function to apply
#' @param ... Additional arguments to fun through \link{mincApplyRCPP} see 
#' notes for a warnings
#' @param mask The path to a mask for the minc files
#' @param tinyMask whether to use a small subset of voxels to test the computation
#' @param batches The number of jobs to break the computation into, ignored for
#' snowfall/local mode
#' @param slab_sizes A 3 element vector indicating how large a chunk of data to read from each minc file 
#' at a time defaults to one slice along the first dimension.
#' @param method A deprecated argument formerly used to configure how to parallelize the jobs
#' now this is handled with \code{conf_file}
#' @param local boolean whether to run the jobs locally (with the parallel package) or with batchtools
#' @param cores defaults to 1 or  
#' \code{max(getOption("mc.cores"), parallel::detectCores() - 1)} if running locally
#' see \link{qMincApply} for details.
#' @param resources A list of resources to request from the queueing system
#' common examples including vmem, walltime, nodes, and modules see
#' \code{system.file("parallel/pbs_script.tmpl", package = "RMINC")} and
#' \code{system.file("parallel/sge_script.tmpl", package = "RMINC")} for
#' more details
#' @param packages Character vector of packages to load for all jobs
#' @param vmem The number of gigabytes of memory to request for each batched
#' job. It is a compatibility argument and will overload \code{vmem} 
#' set in the resource list (if it is defined)
#' @param walltime The amount of walltime to request for each batched job.
#' It is a compatibility argument and will overwrite \code{walltime}
#' set in the resource list (if it is defined)
#' @param workers The number of workers to use. It is a compatibility option
#' and will overwrite \code{batches} if it is supplied.
#' @param temp_dir A path to a temporary directory to hold the job registry
#' created when using a true queuing system and for writing temporary mask
#' files. This must be a location read/writable by all nodes when using a
#' true queuing system (so /tmp will not work).
#' @param cleanup Whether to clean up registry after a queue parallelization job
#' @param collate A function to be applied to collapse the results of the
#' the pMincApply. Defaults to \link{simplify2minc}.
#' @inheritParams qMincApply
#' @return The results of applying \code{fun} to each voxel accross \code{filenames}
#' after collation with \code{collate}
#' @details This is a convenience wrapper for two underlying functions \link{qMincApply}
#' and \link{mcMincApply} for queueing and multicore processing respectively. Each of
#' these functions divides all of the voxels that are masked by \code{mask} into
#' \code{batches}. Batches are processed in parallel, with calling
#' \link{mincApplyRCPP} to process the voxels in the batch. Arguments passed in through
#' \code{...} will be bound by \link{mincApplyRCPP} before \code{fun}, so be wary
#' of potential partial matches. When in doubt, partially apply your function before
#' hand, and do not rely on positional matching.
#' @export
pMincApply <-
  function(filenames, fun, ...,
           mask = NULL,
           tinyMask = FALSE,
           batches = 4,
           slab_sizes = NULL,
           method = NULL,
           local = FALSE,
           cores = NULL,
           resources = NULL,
           packages = NULL,
           vmem = NULL,
           walltime = NULL,
           workers = batches,
           temp_dir = getwd(),
           cleanup = TRUE,
           collate = simplify2minc,
           conf_file = getOption("RMINC_BATCH_CONF"),
           registry_name = new_file("pMincApply_registry"),
           registry_dir = getwd()
  ){
    
    if(!is.null(method))
      warning("Use of method is deprecated, please use a `conf_file` or"
              , "the local argument")
    
    enoughAvailableFileDescriptors(length(filenames))

    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(filenames[1])
      slab_sizes[1] <- 1
    }
    
    if(local){
      if(is.null(cores))
        cores <- max(getOption("mc.cores"), parallel::detectCores() - 1)
      
      results <- mcMincApply(filenames, fun, ...,
                             slab_sizes = slab_sizes,
                             mask = mask, 
                             cores = cores,
                             temp_dir = temp_dir, collate = collate)
    } else {
      if(is.null(cores))
        cores <- 1
      
      results <- qMincApply(filenames, fun, ..., mask = mask, 
                            slab_sizes = slab_sizes,
                            tinyMask = tinyMask, batches = batches, 
                            resources = resources, packages = packages,
                            clobber = TRUE, collate = collate,
                            cleanup = cleanup,
                            conf_file = conf_file,
                            registry_dir = registry_dir,
                            registy_name = registry_name) 
    }

  
    results
  }

#' Local multicore mincApply
#' 
#' Split a minc apply job into batches and process it locally
#' using a fork cluster generated by the parallel package.
#' 
#' @param filenames Paths to the minc files to apply accross
#' @param fun An arbitrary R function to be applied
#' @param ... Additional arguments to pass to fun, see details for a warning
#' @param mask The mask used to select voxels to apply to
#' @param tinyMask Shrink the mask for testing
#' @param slab_sizes A 3 element vector indicating large a chunk of data to read from each minc file 
#' at a time defaults to one slice along the first dimension.
#' @param temp_dir A directory to hold mask files used in the job batching
#' @param cores the number of cores to use, defaults to the option
#' \code{mc.cores} or one less than the number of detected cores if \code{mc.cores} 
#' is unset.
#' @param return_raw An internal use argument that prevents the resulting object
#' from being reordered and expanded.
#' @param cleanup Whether to delete temporary parallelization masks 
#' @param mask_vals values of the mask over which to parallelize, defaults
#' to subdividing all masked voxels into the specified number of batches
#' @param collate A function to collate the list into another object type
#' @export
mcMincApply <-
  function(filenames, fun, ...,
           mask = NULL,
           tinyMask = FALSE,
           slab_sizes = NULL,
           temp_dir = getwd(),
           cores = getOption("mc.cores", parallel::detectCores() - 1),
           return_raw = FALSE,
           cleanup = TRUE,
           mask_vals = NULL,
           collate = simplify2minc
           ){
    
    filenames <- as.character(filenames)
    enoughAvailableFileDescriptors(length(filenames))
    
    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(filenames[1])
      slab_sizes[1] <- 1
    }
    
    sample_file <- filenames[1]
    sample_volume <- mincGetVolume(sample_file)
    mask_file <- "unset"
    
    if(is.null(mask_vals)){
      mask_file <- create_parallel_mask(sample_file = sample_file
                                        , mask = mask
                                        , n = cores
                                        , tinyMask = tinyMask)
      on.exit(try(unlink(mask_file)))
      mask_vals <- 1:cores
    } else {
      mask_file <- mask
    }
    
    dot_args <- list(...)
    mincApplyArguments <- 
      c(list(filenames = filenames,
             fun = match.fun(fun),
             slab_sizes = slab_sizes,
             mask = mask_file,
             collate = identity),
        dot_args)
    #override important argument
    mincApplyArguments$filter_masked <- TRUE
    mincApplyArguments$return_indices <- TRUE
    
    results <- 
      failing_mcmapply(mincApplyRCPP,
                     maskval = mask_vals,
                     MoreArgs = mincApplyArguments,
                     SIMPLIFY = FALSE,
                     mc.cores = cores)
    
    inds <- unlist(lapply(results, function(el) el$inds))
    vals <- unlist(lapply(results, function(el) el$vals), recursive = FALSE)
    
    if(return_raw){
      results <- list(inds = inds, vals = vals)
      results <- setMincAttributes(results, list(filenames = filenames,
                                                 likeVolume = filenames[1],
                                                 mask = mask))
      return(results)
    }
    
    result_order <-
      order(inds)
    
    results <- vals[result_order]
    
    if(!is.null(mask)){
      expanded_results <- rep(list(getOption("RMINC_MASKED_VALUE")), length(sample_volume))
      expanded_results[inds[result_order] + 1] <- results ##add one to inds to convert c++ to R
      results <- expanded_results
    }
    
    collation_function <- match.fun(collate)
    results <- collation_function(results)
    results <- setMincAttributes(results, list(filenames = filenames,
                                               likeVolume = filenames[1],
                                               mask = mask))
    
    return(results)
  }


#' True cluster mincApply
#' 
#' Split a minc apply job into batches and process it either locally
#' or a true grid computing setup. Endeavours to provide an abstract
#' and customizable interface for job scheduling based on the batchtools
#' package. Basic steps of the apply is to 
#' \itemize{
#' \item{create a registry with \link{qMincRegistry} where jobs are coordinated and 
#' results are deposited}
#' \item{map a function over batches of voxels in a collection of minc volumes with
#' \link{qMincMap}, generating submission scripts for the queue scheduler and 
#' submitting the jobs}
#' \item{Collect the results from each batch with \code{qMincReduce}, reorganizing the voxel 
#' results as necessary to reproduce the original order, and then collating the results into a
#' usable object}
#' }
#' Interfaces are provided to perform all three steps at once, either through 
#' \link{qMincApply} or the more general \link{pMincApply}. By default
#' \link{qMincApply} will wait for the cluster to finish processing all jobs
#' although the jobs can be submitted and the R session closed while still
#' maintaining the ability to access results when the jobs have finished. 
#' 
#' @param filenames Paths to the minc files to apply accross
#' @param fun An arbitrary R function to be applied
#' @param ... extra arguments to pass down through \code{qMincMap} to \code{mcMincApply}
#' to \code{mincApplyRCPP} and finally to \code{fun}, there is a chance arguments 
#' here will be trapped by one of the functions on this chain, when in doubt partially
#' apply \code{fun} to its arguments before hand and do not use positional arguments, they
#' are almost certainly not going to work as expected. 
#' @param mask The mask used to select voxels to apply to
#' @param tinyMask Shrink the mask for testing
#' @param slab_sizes A 3 element vector indicating large a chunk of data to read from each minc file 
#' at a time defaults to one slice along the first dimension.
#' @param batches The number of batches to divide the job into, this is ignored for
#' multicore jobs, with the number of batches set to the number of cores.
#' @param cores the number of cores to parallelize across for each worker, defaults to 1
#' but higher numbers may be useful for batchtools multicore or systems like SciNet that do
#' not allocate single core jobs.
#' @param resources The resources to request for each job, overrides the \code{default.resources}
#' specified in the configuration list. 
#' @param packages packages to be loaded for each job in a registry
#' @param temp_dir A directory to store files needed for the parallelization
#' and job management
#' @param wait Whether to wait for your results or return a registry object
#' to be checked on later
#' @param cleanup Whether to empty the registry after a successful run defaults
#' to true
#' @param clobber Whether to overwrite an existing registry at \code{registry_dir} 
#' @param collate A function to collate the returned list into another object type.
#' @param ignore_incompletes Whether to reduce the results with \code{qMincReduce}
#' even if all jobs are not complete.
#' @param registry_dir where batchtools should create the registry
#' @param registry_name a name for the registry
#' @param registry a pre-existing job registry
#' @param conf_file A batchtools config file, defaults to \code{option("RMINC_BATCH_CONF")}
#' @details RMINC's batching facilities are inherited with little modification from
#' the batchtools package, mostly just providing handy wrappers to handle registry
#' creation, batching, submission, and reduction. The abstractions provided are very leaky
#' and it is worth learning about batchtools to handle more complex situations. Formerly one
#' could set the parallelization method from this function, this has been removed. 
#' Controlling how and where to execute the parallel jobs is now handled by 
#' the conf_file argument.
#' @return 
#' \itemize{
#' \item{If \code{qMincApply} is called with \code{wait = TRUE} or if \code{qMincReduce}
#' is called, the results are returned after collation with \code{collate}
#' }
#' \item{If \code{qMincApply} is called with \code{wait = FALSE} or if \code{qMincRegistry} or
#' \code{qMincMap} are called  a batchtools registry is returned that can be used to 
#' query job states, kill jobs, and collected results
#' }
#' }
#' @seealso \url{https://www.jstatsoft.org/article/view/v064i11} \link{pMincApply} \link{mcMincApply}
#' @export
qMincApply <- 
  function(filenames, fun, ..., 
           mask=NULL, batches=4, tinyMask=FALSE,
           slab_sizes = NULL,
           resources = NULL,
           packages = c("RMINC"),
           registry_dir = getwd(),
           registry_name = "qMincApply_registry",
           temp_dir = getwd(),
           cores = 1, #max(getOption("mc.cores"), parallel::detectCores() - 1),
           wait = TRUE,
           cleanup = TRUE,
           clobber = FALSE,
           collate = simplify2minc,
           conf_file = getOption("RMINC_BATCH_CONF")) {

    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(filenames[1])
      slab_sizes[1] <- 1
    }

    qMinc_registry <-
      qMincRegistry(registry_name = new_file(registry_name, registry_dir),
                    registry_dir = registry_dir,
                    packages = packages,
                    clobber = clobber,
                    resources = resources,
                    conf_file = conf_file)
    
    on.exit({
      if(cleanup && wait){
        try(lapply(list.files(temp_dir, pattern = "pMincMask.*\\.mnc"), unlink))
        tenacious_remove_registry(qMinc_registry)
      } 
    })
    
    qMincMap(qMinc_registry,
             filenames, 
             fun = match.fun(fun), 
             ..., 
             slab_sizes = slab_sizes,
             batches = batches,
             cores = cores,
             mask = mask, 
             tinyMask = tinyMask,
             temp_dir = temp_dir)
    
    if(wait){
      qMinc_results <- qMincReduce(qMinc_registry, wait = TRUE, collate = collate)
      return(qMinc_results)
    }
    
    return(qMinc_registry)
  }


#' @describeIn qMincApply registry 
#' @export
qMincRegistry <- function(registry_name = "qMincApply_registry",
                          packages = c("RMINC"),
                          registry_dir = getwd(),
                          clobber = FALSE,
                          resources = NULL,
                          conf_file = getOption("RMINC_BATCH_CONF")){
  
  if(! "RMINC" %in% packages)
    packages <- c("RMINC", packages)
  
  if(clobber)
    try(removeRegistry(reg = loadRegistry(registry_dir)), silent = TRUE)
  
  qMinc_registry <-
    makeRegistry(registry_name,
                 registry_dir,
                 packages = packages,
                 conf.file = conf_file)

  if(!is.null(resources))
    qMinc_registry$default.resources <- resources
  
  return(qMinc_registry)
}

#' @describeIn qMincApply map
#' @export
qMincMap <- 
  function(registry, filenames, fun, ..., mask = NULL, 
           slab_sizes = NULL,
           batches = 4, tinyMask = FALSE, temp_dir = getwd(),
           resources = NULL,
           cores = 1){
    
    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(filenames[1])
      slab_sizes[1] <- 1
    }
    
    new_mask_file <- create_parallel_mask(sample_file = filenames[1]
                                          , mask = mask
                                          , n = batches * cores
                                          , tinyMask = tinyMask)
    
    #Group the indices
    mask_indices <- 
      lapply(1:batches, function(i){
        seq((i -1)*cores + 1,
            i * cores)
      })
    
    #Create a list of all additional args to pass to mcMincApply
    dot_args <- list(...)
    mincApplyArguments <- 
      c(list(filenames = filenames, 
             fun = match.fun(fun),
             slab_sizes = slab_sizes,
             cores = cores,
             mask = new_mask_file,
             tinyMask = tinyMask,
             temp_dir = temp_dir,
             collate = identity),
        dot_args)
    #Override these if passed through ...
    mincApplyArguments$return_raw <- TRUE

    
    ids <-
      batchMap(reg = registry,
               mcMincApply, 
               mask_vals = mask_indices,
               more.args = mincApplyArguments)
    
    submitJobs(ids, reg = registry)
    
    return(registry)
  }

#' @describeIn qMincApply reduce
#' @export
qMincReduce <- 
  function(registry, ignore_incompletes = FALSE, wait = FALSE, collate = simplify2minc){
    
    if(wait)
      waitForJobs(reg = registry)
    
    if((!ignore_incompletes) && length(findNotDone(registry) != 0))
      stop("Some jobs have not terminated, use `ignore_incompletes` to reduce anyway, or set `wait`")
    
    results <- reduceResultsList(reg = registry, use.names = FALSE)
    result_attributes <- mincAttributes(results[[1]])
    
    result_indices <- unlist(lapply(results, function(el) el$inds))
    result_order <- order(result_indices)
    results <- unlist(lapply(results, function(el) el$vals), recursive = FALSE)[result_order]
    
    if(!is.null(result_attributes$mask)){
      total_voxels <- prod(minc.dimensions.sizes(result_attributes$mask))
      expanded_results <- rep(list(getOption("RMINC_MASKED_VALUE")), total_voxels)
      ##add one to indices to convert from c++ to R
      expanded_results[result_indices[result_order] + 1] <- results
      results <- expanded_results
    }
    
    collation_function <- match.fun(collate)
    results <- collation_function(results)
    
    if(!is.null(result_attributes))
      results <- setMincAttributes(results, result_attributes)
    
    return(results)
  }

groupingVector <- function(len, groups){
  group_len <- ceiling(len / groups)
  n_overfill <- (group_len * groups) - len
  groups <- rep(1:groups, each = group_len)
  if(n_overfill != 0)
    groups <- groups[-(1:n_overfill)]
  
  return(groups)
}

### Legacy code for pMincApply 
# pMincApply <- 
#   function(filenames, function.string,
#            mask=NULL, workers=4, tinyMask=FALSE, 
#            method="snowfall",global="",packages="", 
#            modules="",vmem="8",walltime="01:00:00") {
#   
#   REDUCE = TRUE; # For now this option is not exposed
# 
#   # if no mask exists use the entire volume
#   if (is.null(mask)) {
#     maskV = mincGetVolume(filenames[1])
#     nVoxels = length(maskV)
#     maskV[maskV >= min(maskV)] <- as.integer(cut(seq_len(nVoxels), workers)) 
#   }
#   else {
#     maskV <- mincGetVolume(mask)
#     # optionally make the mask a fraction of the original size - for testing
#     if (tinyMask!=FALSE) {
#       maskV[maskV>1.5] <- 0
#     }
#     nVoxels <- sum(maskV>0.5)
#     maskV[maskV>0.5] <- as.integer(cut(seq_len(nVoxels), workers)) 
#   }
#   
#   # Saving to /tmp does not always work...
#   maskFilename <- paste("pmincApplyTmpMask-", Sys.getpid(), ".mnc", sep="")
#   
#   #If the current working directory isn't writeable, 
#   #write to a tempdir instead
#   if(file.access(getwd(), 2) != 0) maskFilename <- file.path(tempdir(), maskFilename)
#   
#   mincWriteVolume(maskV, 
#                   maskFilename, 
#                   clobber=TRUE) 
#   
#   
#   # create the packageList that will be used for the snowfall and sge options
#   # if packages contains multiple libraries, the test (packages == "") 
#   # will return as many TRUE/FALSE as the length of the vector. So to test
#   # for "", first test that the length of the packages vector is 1
#   if(length(packages) < 2) {
#     if(packages == "") {
#       packageList = c("RMINC")
#     }
#     else {
#       packageList = c(packages,"RMINC")
#     }
#   }
#   else {
#     packageList = c(packages,"RMINC")
#   }
#   
#   pout <- list()
#   
#   if (method == "local") {
#     stop("Lovely code ... that generates inconsistent results because something somewhere is not thread safe ...")
#     
#     if(!(requireNamespace("doMC", quietly = TRUE) & requireNamespace("foreach", quietly = TRUE))) 
#       stop("One or both of doMC and foreach is missing, please install these packages")
#     
#     registerDoMC(workers)
#     
#     # run the job spread across each core
#     pout <- foreach(i=1:workers) %dopar% { mincApply(filenames, function.string,
#                                                    mask=maskFilename, maskval=i) }
#     #cat("length: ", length(pout), "\n")
#   }
# 
#   # The pbs options use mpirun and snow to parallelize mincApply over multiple cores.
#   # It is currently configured to send all the jobs to one node, and parallelize over the
#   # cores at that node. Also note the amount of virtual memory requested is set at 8g.
#  
#   # It operates as follows:
#   # 1) Save global variables to disk
#   # 2) Write out a .R file that will execute mpi operations once submitted to the cluster
#   # 3) Write out a .sh file that will be submitted to the cluster
#   # 4) Submit the .sh file to the cluster
#   # 5) Wait for the jobs to finish
#   # 6) Read the output from disk
# 
#   else if (method == "pbs") {
#  
#    if(is.null(getOption("MAX_NODES"))) {
# 		 ppn = 8
#                  #maximize node usage
# 	 	 nodes = ceiling(workers/ppn)
# 		 workers = nodes*ppn-1
#    }
#    else {
# 		 ppn = getOption("MAX_NODES")
# 		 nodes = ceiling(workers/ppn)
#     }
# 
#     # 1) Save variables which will be referenced on the cluster to disk (including user specified global variables)
#     rCommand = 'save(\'maskV\',\'filenames\',\'workers\',\'REDUCE\',\'function.string\',\'maskFilename\','
#     for (nVar in 1:length(global)) {
# 	if(global[nVar] != "") {
# 		rCommand = paste(rCommand,'\'',global[nVar],'\',',sep="")
#         }
#     }
#     rCommand = paste(rCommand,'file=\'mpi-rminc-var\')',sep="")
#     
#     eval(parse(text=rCommand))
# 
#     # 2) Write out an R file to disk. This file will be executed via mpirun on the cluster
#     fileConn <- file("mpi-rminc.R",open='w')
# 
#     # Load R Packages
#     
# 
#     # snow is used to coordinate operations, but Rmpi could be used for greater control
#     packageList = c(packageList,"snow")
# 
#     for (nPackage in 1:length(packageList)) {     
#     	writeLines(paste("library( ",packageList[nPackage],")",sep = ""),fileConn)
#     }
# 
#     # Load the variables we saved in step 1
#     writeLines("load(\"mpi-rminc-var\")",fileConn)
# 
#     # Create the snow cluster
#     writeLines("cl <- makeCluster(workers, type = \"MPI\")",fileConn) 
# 
#     # Create a wrapper function that we can easily run with clusterApply
#     writeLines("wrapper <- function(i) { return(mincApply(filenames, function.string, mask = maskFilename, maskval = i, reduce = REDUCE))}",fileConn)
# 
#     # Export all neccessary variables to each slave
#     writeLines("clusterExport(cl,c('filenames','REDUCE','function.string','maskFilename','maskV'))",fileConn)
# 
#     # Main mpi exection
#     writeLines("clusterOut <- clusterApply(cl,1:workers,wrapper)",fileConn)
# 
#     # At this point we are done.
#     writeLines("stopCluster(cl)",fileConn)
# 
#     # Test data at one voxel to determine how many ouytputs
#     writeLines(" x <- mincGetVoxel(filenames, 0,0,0)",fileConn)
#     writeLines("test <- eval(function.string)",fileConn)
#     writeLines("if (length(test) > 1) {",fileConn)
#     writeLines("output <- matrix(0, nrow=length(maskV), ncol=length(test))",fileConn)
#     writeLines("class(output) <- class(clusterOut[[1]])",fileConn)
#     writeLines("attr(output, \"likeVolume\") <- attr(clusterOut[[1]], \"likeVolume\")",fileConn)
#     writeLines("} else {",fileConn)
#     writeLines("output <- maskV",fileConn)
#     writeLines("}",fileConn)
#     writeLines("for(i in 1:workers) {",fileConn)
#     writeLines("if (length(test)>1) {",fileConn)
#     writeLines("if(REDUCE == TRUE)",fileConn)	
#     writeLines("output[maskV == i,] <- clusterOut[[i]]",fileConn)
#     writeLines("else",fileConn)
#     writeLines("output[maskV == i,] <- clusterOut[[i]][maskV == i, ]",fileConn)
#     writeLines("}",fileConn)
#     writeLines("else {",fileConn)
#     writeLines("output[maskV==i] <- clusterOut[[i]]",fileConn)
#     writeLines("}",fileConn)
#     writeLines("}",fileConn)
# 
#     # Write the output (R data) to disk, so the initiating R session can access the data
#     writeLines("save('output',file='mpi-rminc-out')",fileConn)
# 
#     close(fileConn)
# 
#     # 3) Write out a .sh file which will be what is submitted to the cluster
#     qfileConn <- file("q-mpi-rminc.sh",open='w')
#     writeLines("#!/bin/bash -x",qfileConn)
#    
#     # Errors and Output are written to the current directory, but this could be set to a user defined spot.
#     writeLines("#PBS -e ./",qfileConn)
#     writeLines("#PBS -o ./",qfileConn)
#     writeLines("#PBS -N pMincApply",qfileConn) 
# 
#     # Allocate nodes and ppn
#     writeLines(paste("#PBS -l nodes=",as.character(nodes),":ppn=",as.character(ppn),sep=""),qfileConn)
# 
#     # Allocate walltime and vmem
#     writeLines(paste("#PBS -l vmem=",vmem,"g,walltime=",walltime,sep=""),qfileConn)
# 
#     # Load modules
#     for (nModule in 1:length(modules)) {     
#         if(modules[nModule] != "") {
#     		writeLines(paste("module load ",modules[nModule],sep = ""),qfileConn)
#        }
#     }
# 
#     # Define temp directory
#     writeLines(paste("export TMPDIR=",getOption("TMPDIR"),sep=""),qfileConn)
# 
#    # For Testing
#     writeLines(paste("cp -R ~/Software/RMINC/rminctestdata /tmp/",sep=""),qfileConn)
# 
#     # Move to working directory
#     writeLines(paste("cd ",getOption("WORKDIR"),sep=""),qfileConn)
# 
#     # Neccessary to initiate mpi operations.
#     writeLines("mpirun -np 1 R CMD BATCH mpi-rminc.R",qfileConn)
# 
#     close(qfileConn)
# 
#     # 4) Submit the job
#     result <- system("qsub q-mpi-rminc.sh",intern=TRUE)
#     ptm = proc.time()
#     
#     # 5) Wait for Job to finish
#     status <- system(paste("qstat ",result," | grep C"),intern=TRUE) 
#     # 1 Indicates not completed
#     while(length(status) == 0) {
#        flush.console()
#        status <- system(paste("qstat ",result," | grep R"),intern=TRUE)
#        runTime = proc.time()-ptm
#        if(length(status) == 0) {
# 	  cat(paste("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bQueued:  ",sprintf("%3.3f seconds",(runTime[3]))))
#        } else {
# 	  cat(paste("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRunning: ",sprintf("%3.3f seconds",(runTime[3]))))
#        }
#        Sys.sleep(2)
#        status <- system(paste("qstat ",result," | grep C"),intern=TRUE)
#       
#     }
#     print(paste("Completed: ",as.character(runTime[3])," Seconds"))
#     # 6) Read output from disk
#     load("mpi-rminc-out")
#     
#     # Clean up intermediate files
#     system(paste("rm",maskFilename))
#     system("rm q-mpi-rminc.sh")
#     system("rm mpi-rminc.R")
#     system("rm mpi-rminc-var")
#     system("rm mpi-rminc-out")
# 
#     return(output)
# 
#   }
#   
#   else if (method == "snowfall") {
#     
#     if(!requireNamespace("snowfall", quietly = TRUE))
#       stop("The snowfall package is required to run code with snowfall parallelism")
#     
#     sfInit(parallel=TRUE, cpus=workers)
# 
#     for (nPackage in 1:length(packageList)) {
#       sfLibrary(packageList[nPackage],character.only=TRUE) 
#     }
# 		
#     sfExport(list = global) 
# 
#     wrapper <- function(i) {
#       cat( "Current index: ", i, "\n" ) 
#       return(mincApply(filenames, function.string, mask=maskFilename,
#                       maskval=i, reduce=REDUCE))
#     }
#     # use all workers in the current cluster if # of workers not specified
#     if (is.null(workers)) {
#       workers <- length(sfSocketHosts())
#     }
#     
#     #sink("/dev/null");
#     pout <- sfLapply(1:workers, wrapper)
#     
#     sfStop();
#   }
#   else {
#     stop("unknown execution method")
#   }
#   
#   # Need to get one voxel, x, to test number of values returned from function.string
#   x <- mincGetVoxel(filenames, 0,0,0)
#   test <- eval(function.string) 
#   
#   # recombine the output into a single volume
#   if (length(test) > 1) {
#     output <- matrix(0, nrow=length(maskV), ncol=length(test))
#     class(output) <- class(pout[[1]])
#     attr(output, "likeVolume") <- attr(pout[[1]], "likeVolume")
#   }
#   else {
#     output <- maskV
#   }
#   
#   for(i in 1:workers) {
#     if (length(test)>1) {
#       if(REDUCE == TRUE)	
#       	output[maskV == i,] <- pout[[i]]
#       else
#       	output[maskV == i,] <- pout[[i]][maskV == i, ] 
#     }
#     else {
#       output[maskV==i] <- pout[[i]]
#     }
#   }
#   unlink(maskFilename)
#   return(output)
# }
