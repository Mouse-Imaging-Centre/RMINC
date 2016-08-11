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
#' @param slab_sizes A 3 element vector indicating large a chunk of data to read from each minc file 
#' at a time defaults to one slice along the first dimension.
#' @param method The parallelization method, local and snowfall perform 
#' computations on multiple cores (snowfall is an alias to local for backcompatibility)
#' pbs refers to a torque queueing system, and sge refers to the
#' sge queueing system.
#' @param cores defaults to 1 in the case of PBS or SGE jobs, and 
#' \code{max(getOption("mc.cores"), parallel::detectCores() - 1)} for local/snowfall
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
           method = c("local", "snowfall", "pbs", "sge", "none"),
           cores = getOption("mc.cores", parallel::detectCores() - 1),
           resources = list(),
           packages = NULL,
           vmem = NULL,
           walltime = NULL,
           workers = batches,
           temp_dir = tempdir(),
           cleanup = TRUE,
           collate = simplify2minc
  ){
    
    enoughAvailableFileDescriptors(length(filenames))
    method <- match.arg(method)
    
    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(filenames[1])
      slab_sizes[1] <- 1
    }
    
    if(method == "local" || method == "snowfall")
      results <- mcMincApply(filenames, fun, ...,
                             slab_sizes = slab_sizes,
                             mask = mask, tinyMask = tinyMask, 
                             cores = cores,
                             temp_dir = temp_dir, collate = collate)
    
    #Initialize queue resources from arguments if not supplied
    if(!is.null(vmem)) resources$vmem <- vmem
    if(!is.null(walltime)) resources$walltime <- walltime
    
    if(method == "sge")
      results <- qMincApply(filenames, fun, ..., mask = mask, 
                            slab_sizes = slab_sizes,
                            tinyMask = tinyMask, batches = batches, 
                            resources = resources, packages = packages,
                            clobber = TRUE, queue = "sge", collate = collate,
                            cleaup = cleanup)
    
    if(method == "pbs")
      results <- qMincApply(filenames, fun, ..., mask = mask,
                            slab_sizes = slab_sizes,
                            tinyMask = tinyMask, batches = batches, 
                            resources = resources, packages = packages,
                            clobber = TRUE, queue = "pbs", collate = collate,
                            cleanup = cleanup)
    
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
           temp_dir = tempdir(),
           cores = getOption("mc.cores", parallel::detectCores() - 1),
           return_raw = FALSE,
           cleanup = TRUE,
           mask_vals = NULL,
           collate = simplify2minc){
    
    filenames <- as.character(filenames)
    enoughAvailableFileDescriptors(length(filenames))
    
    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(lmod$fr[1,1])
      slab_sizes[1] <- 1
    }
    
    sample_file <- filenames[1]
    sample_volume <- mincGetVolume(sample_file)
    mask_file <- "unset"
    
    if(is.null(mask_vals)){
      temp_dir <- path.expand(temp_dir)
      if(!file.exists(temp_dir)) dir.create(temp_dir)
      
      mask_file <- tempfile("pMincMask", tmpdir = temp_dir, fileext = ".mnc")
      
      on.exit({
        if(cleanup) unlink(mask_file)
      })
      
      if(is.null(mask)){
        if(cores %% 1 != 0) stop("the number of cores must be an integer")
        
        nVoxels <- length(sample_volume)
        mask_values <- groupingVector(nVoxels, cores)
      } else {
        mask_values <- mincGetVolume(mask)
        if(tinyMask) mask_values[mask_values > 1.5] <- 0
        
        nVoxels <- sum(mask_values > .5)
        
        mask_values[mask_values > .5] <- 
          groupingVector(nVoxels, cores)
        
      }
      
      mincWriteVolume(mask_values, mask_file, like.filename = sample_file)
    }
    
    if(mask_file == "unset"){
      mask_file <- mask
      on.exit() #clear the on.exit for extra safety
    }
    
    if(is.null(mask_vals))
      mask_vals <- 1:cores
    
    dot_args <- list(...)
    mincApplyArguments <- 
      c(list(filenames = filenames), #list wrapping ensures c() works
        fun = match.fun(fun),
        dot_args,
        slab_sizes = slab_sizes,
        mask = mask_file,
        collate = identity)
    #override important argument
    mincApplyArguments$filter_masked <- TRUE
    mincApplyArguments$return_indices <- TRUE
    
    results <- 
      parallel::mcmapply(mincApplyRCPP,
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
#' and customizable interface for job scheduling based on the BatchJobs
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
#' @param parallel_method Where to run the batches, defaults to \code{"none"} which results
#' in using the currently loaded configuration for BatchJobs. The other options
#' are \code{"sge"} for Sun Grid Engine, \code{"pbs"} for Torque, \code{"multicore"} for
#' local (not recommended, use \link{mcMincApply}), \code{"interactive"} for sequential local
#' processing, and \code{"custom"} in which case a cluster.functions object must be
#' passed in through \code{cluster_functions} (also not recommend, for more control 
#' use the .BatchJobs.R configuration system, see details).
#' @param cluster_functions Custom cluster functions to use if parallel_method = "custom"
#' @param batches The number of batches to divide the job into, this is ignored for
#' multicore jobs, with the number of batches set to the number of cores.
#' @param cores the number of cores to use, in all cases except multicore this
#' defaults to 1, with a reminder message if method "none" or "custom" are used,
#' for multicore it defaults to \code{max(getOption("mc.cores"), parallel::detectCores() - 1)}
#' @param resources The resources to request for each job, overrides the default.resources
#' specified in .BatchJobs.R config files (see details). These include things like vmem, 
#' walltime, and nodes.
#' @param packages packages to be loaded for each job in a registry
#' @param temp_dir A directory to store files needed for the parallelization
#' and job management
#' @param registry_name The name to give your BatchJobs registry
#' @param registry_dir The path for your BatchJobs registry, defaults to
#' \code{temp_dir/registry_name}
#' @param wait Whether to wait for your results or return a registry object
#' to be checked on later
#' @param cleanup Whether to empty the registry after a successful run defaults
#' to true
#' @param clobber Whether to overwrite an existing registry at \code{registry_dir} 
#' @param collate A function to collate the returned list into another object type.
#' @param ignore_incompletes Whether to reduce the results with \code{qMincReduce}
#' even if all jobs are not complete.
#' @param registry A pre-made BatchJobs registry
#' @details RMINC's batching facilities are inherited with little modification from
#' the BatchJobs package, mostly just providing handy wrappers to handle registry
#' creation, batching, submission, and reduction. The abstractions provided are very leaky
#' and it is worth learning about BatchJobs to handle more complex situations. This being
#' said a high degree of flexibility is already available. RMINC honours the hierarchy
#' of configuration files <BatchJobs>/.BatchJobs.R < ~/.BatchJobs.R < getwd()/.BatchJobs.R
#' these files are the place set sensible defaults for your usage. The template scripts
#' provided for SGE and PBS are located in \code{system.file("parallel/", package = "RMINC")}
#' as well as some example .BatchJobs.R files.  
#' @return 
#' \itemize{
#' \item{If \code{qMincApply} is called with \code{wait = TRUE} or if \code{qMincReduce}
#' is called, the results are returned after collation with \code{collate}
#' }
#' \item{If \code{qMincApply} is called with \code{wait = FALSE} or if \code{qMincRegistry} or
#' \code{qMincMap} are called  a BatchJobs registry is returned that can be used to 
#' query job states, kill jobs, and collected results
#' }
#' }
#' @seealso \url{https://www.jstatsoft.org/article/view/v064i11} \link{pMincApply} \link{mcMincApply}
#' @export
qMincApply <- 
  function(filenames, fun, ..., 
           mask=NULL, batches=4, tinyMask=FALSE,
           slab_sizes = NULL,
           parallel_method = c("none", "sge", "pbs", "multicore", 
                               "interactive", "custom"),
           cluster_functions = NULL,
           resources = list(),
           packages = c("RMINC"),
           temp_dir = tempdir(),
           registry_name = "qMincApply_registry",
           registry_dir = file.path(temp_dir, registry_name),
           cores = NULL, #max(getOption("mc.cores"), parallel::detectCores() - 1),
           wait = TRUE,
           cleanup = TRUE,
           clobber = FALSE,
           collate = simplify2minc) {
    
    parallel_method <- match.arg(parallel_method)
    
    if(parallel_method %in% c("sge", "pbs") &&
       missing(temp_dir))
      stop("When using a sge or pbs systems using the default temp_dir is unlikely to work.",
           "Manually specify temp_dir = tempdir() to try anyway")
    
    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(lmod$fr[1,1])
      slab_sizes[1] <- 1
    }
    
    if(is.null(cores)){
      cores <-
        switch(parallel_method,
               sge = 1,
               pbs = 1,
               multicore = max(getOption("mc.cores"), parallel::detectCores() - 1),
               interactive = 1,
               {message("Defaulting to 1 core for parallel method ", parallel_method); 1})
    }
    
    qMinc_registry <-
      qMincRegistry(registry_name = registry_name,
                    registry_dir = registry_dir,
                    parallel_method = parallel_method,
                    cluster_functions = cluster_functions,
                    packages = packages,
                    cores = cores,
                    clobber = clobber)
    
    on.exit({
      if(cleanup && wait){
        removeRegistry(qMinc_registry, ask = "no")
        lapply(list.files(temp_dir, pattern = "pMincMask.*\\.mnc"), unlink)
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
             temp_dir = temp_dir,
             resources = resources,
             parallel_method = parallel_method)
    
    if(wait){
      qMinc_results <- qMincReduce(qMinc_registry, wait = TRUE, collate = collate)
      return(qMinc_results)
    }
    
    return(qMinc_registry)
  }


#' @describeIn qMincApply registry 
#' @export
qMincRegistry <- function(registry_name = "qMincApply_registry",
                          parallel_method = c("none", "sge", "pbs", 
                                              "multicore", "interactive", "custom"),
                          packages = c("RMINC"),
                          cluster_functions = NULL,
                          registry_dir = file.path(tempdir(), registry_name),
                          cores = NULL,
                          clobber = FALSE
){
  
  config <- getConfig()
  
  parallel_method <- match.arg(parallel_method)
  script_directory <- system.file("parallel", package = "RMINC")
  
  if(is.null(cores)){
    cores <-
      switch(parallel_method,
             sge = 1,
             pbs = 1,
             multicore = max(getOption("mc.cores"), parallel::detectCores() - 1),
             interactive = 1,
             {message("Defaulting to 1 core for parallel method ", parallel_method); 1})
  }
  
  config$cluster.functions <-
    switch(parallel_method,
           "none" = config$cluster.functions,
           "pbs" = makeClusterFunctionsTorque(file.path(script_directory, "pbs_script.tmpl")),
           "sge" = makeClusterFunctionsSGE(file.path(script_directory, "sge_script.tmpl")),
           "multicore" = makeClusterFunctionsMulticore(max.jobs = cores),
           "interactive" = makeClusterFunctionsInteractive(),
           "custom" = cluster_functions)
  
  setConfig(config)
  
  if(! "RMINC" %in% packages)
    packages <- c("RMINC", packages)
  
  if(clobber)
    try(removeRegistry(loadRegistry(registry_dir), ask = "no"), silent = TRUE)
  
  qMinc_registry <-
    makeRegistry(registry_name,
                 registry_dir,
                 packages = packages)
  
  return(qMinc_registry)
}

#' @describeIn qMincApply map
#' @export
qMincMap <- 
  function(registry, filenames, fun, ..., mask = NULL, 
           slab_sizes = NULL,
           parallel_method = c("none", "sge", "pbs", 
                               "multicore", "interactive", "custom"),
           batches = 4, tinyMask = FALSE, temp_dir = tempdir(),
           resources = list(),
           cores = NULL){
    
    parallel_method <- match.arg(parallel_method)
    
    if(is.null(slab_sizes)){
      slab_sizes <- minc.dimensions.sizes(lmod$fr[1,1])
      slab_sizes[1] <- 1
    }
    
    #Determine default number of cores
    if(is.null(cores)){
      cores <-
        switch(parallel_method,
               sge = 1,
               pbs = 1,
               multicore = max(getOption("mc.cores"), parallel::detectCores() - 1),
               interactive = 1,
               {message("Defaulting to 1 core for parallel method ", parallel_method); 1})
    }
    
    # If multicore batches is the number of cores
    if(parallel_method == "multicore"){
      batches <- cores
      cores <- 1 
    }
    
    #Read in a likefile
    sample_file <- filenames[1]
    
    temp_dir <- path.expand(temp_dir)
    if(!file.exists(temp_dir)) dir.create(temp_dir)
    
    sample_volume <- mincGetVolume(sample_file)
    new_mask_file <- tempfile("pMincMask", tmpdir = temp_dir, fileext = ".mnc")
    
    
    if(is.null(mask)){
      if(batches %% 1 != 0 ||
         cores %% 1 != 0) 
        stop("the number of batches and cores must be integers")
      
      nVoxels <- length(sample_volume)
      mask_values <- groupingVector(nVoxels, batches * cores)
    } else {
      mask_values <- mincGetVolume(mask)
      if(tinyMask) mask_values[mask_values > 1.5] <- 0
      
      nVoxels <- sum(mask_values > .5)
      mask_values[mask_values > .5] <- 
        groupingVector(nVoxels, batches * cores)
    }
    
    mincWriteVolume(mask_values, new_mask_file, like.filename = sample_file)
    
    #Group the indices
    mask_indices <- 
      lapply(1:batches, function(i){
        seq((i -1)*cores + 1,
            i * cores)
      })
    
    #Create a list of all additional args to pass to mcMincApply
    dot_args <- list(...)
    mincApplyArguments <- 
      c(list(filenames = filenames), #list wrapping ensures c() works
        fun = match.fun(fun),
        dot_args,
        slab_sizes = slab_sizes,
        cores = cores,
        mask = new_mask_file,
        tinyMask = tinyMask,
        temp_dir = temp_dir,
        collate = identity)
    #Override these if passed through ...
    mincApplyArguments$return_raw <- TRUE
    
    batchMap(registry,
             mcMincApply, 
             mask_vals = mask_indices,
             more.args = mincApplyArguments)
    
    submitJobs(registry, resources = resources)
    
    return(registry)
  }

#' @describeIn qMincApply reduce
#' @export
qMincReduce <- 
  function(registry, ignore_incompletes = FALSE, wait = FALSE, collate = simplify2minc){
    
    if(wait)
      waitForJobs(registry)
    
    if((!ignore_incompletes) && length(findNotTerminated(registry) != 0))
      stop("Some jobs have not terminated, use `ignore_incompletes` to reduce anyway, or set `wait`")
    
    results <- loadResults(registry, use.names = FALSE)
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
