anatGetFile <- function(filename, 
                        atlas, 
                        method = "jacobians", 
                        defs = getOption("RMINC_LABEL_DEFINITIONS"), 
                        dropLabels=FALSE, 
                        side = "both" ) {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  out <- NULL
  tmpfile <- tempfile(pattern="RMINC-", fileext=".txt")
  if (method == "jacobians") {
    system(paste("label_volumes_from_jacobians", atlas, filename, ">", tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "labels") {
    # filename here should be a set of labels unique to this brain
    system(paste("volumes_from_labels_only", filename, tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "means") {
    system(paste("compute_values_across_segmentation", "-m",
                 filename, atlas, tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "sums") {
    system(paste("compute_counts_for_labels",
                  atlas, filename, ">", tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "slow_sums") {
    system(paste("compute_values_across_segmentation", "-s",
                 filename, atlas, tmpfile, sep=" "))
    out <- read.csv(tmpfile, header=FALSE)
  }
  else if (method == "text") {
    # values are already extracted and stored in a text file
    out <- read.table(filename, header=FALSE)
  }
  else {
    # unrecognized option...
    stop("Unrecognized option used for \"method\" (anatGetFile/anatGetAll). Available options are: jacobians, labels, means, sums, text.")
  }
  #cat("FILENAME:", filename, "\n")
  if (dropLabels == TRUE) {
    labels <- read.csv(defs)
    usedlabels <- 0
    if (side == "right") {
      usedlabels <- labels$right.label
    }
    else if (side == "left") {
      usedlabels <- labels$left.label
    }
    else if (side == "both") {
      usedlabels <- c(labels$left.label, labels$right.label)
    }
    out <- out[out$V1 %in% usedlabels,]
  }
  on.exit(unlink(tmpfile))
  return(out)
}

anatRenameRows <- function(anat, defs=getOption("RMINC_LABEL_DEFINITIONS")) {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  defs <- read.csv(defs)
  rn <- rownames(anat)
  on <- as.character(rn)
  for (i in 1:length(rn)) {
    # if structure exists in both hemisphere with same number
    if (rn[i] %in% defs$right.label & rn[i] %in% defs$left.label) {
      index <- match(rn[i], defs$right.label)
      on[i] <- paste(defs$Structure[index])
      #cat("Match", rn[i], "to both:", on[i], "\n")
    }
    # left side only
    else if (rn[i] %in% defs$left.label) {
      index <- match(rn[i], defs$left.label)
      on[i] <- paste("left", defs$Structure[index])
      #cat("Match", rn[i], "to left:", on[i], "\n")
    }
    else if (rn[i] %in% defs$right.label) {
      index <- match(rn[i], defs$right.label)
      on[i] <- paste("right", defs$Structure[index])
      #cat("Match", rn[i], "to right:", on[i], "\n")
    }
    else {
      #cat("Something weird", rn[i], "\n")
    }
  }
  rownames(anat) <- on
  # make it an anat class
  class(anat) <- c("anatUnilateral", "anatMatrix", "matrix")
  # make sure we can map label numbers for later use
  attr(anat, "anatIDs") <- rn
  return(anat)
}

#' Faster AnatGet
#' 
#' Computes volumes, means, sums, and similar values across a
#' segmented atlas
#' @inheritParams anatGetAll_old
#' @param strict check if any files differ in step sizes
#' @param parallel how many processors to run on (default=single processor).
#' Specified as a two element vector, with the first element corresponding to
#' the type of parallelization, and the second to the number
#' of processors to use. For local running set the first element to "local" or "snowfall"
#' for back-compatibility, anything else will be run with batchtools see \link{pMincApply}
#' and \link{configureMincParallel} for details.
#' Leaving this argument NULL runs sequentially.
#' @param conf_file A batchtools configuration file defaulting to \code{getOption("RMINC_BATCH_CONF")}

#' @details 
#' anatGetAll needs a set of files along with an atlas and a set of
#' atlas definitions. In the end it will produce one value per label
#' in the atlas for each of the input files. How that value is
#' computed depends on the methods argument:
#' \itemize{
#'   \item{jacobians -}{ Each file contains log jacobians, and the volume for
#'   each atlas label is computed by multiplying the jacobian with
#'   the voxel volume at each voxel.
#'   }
#'   \item{labels -}{ Each file contains integer labels (i.e. same as the atlas).
#'   The volume is computed by counting the number of voxels with
#'   each label and multiplying by the voxel volume.
#'   }
#'   \item{means -}{ Each file contains an arbitrary number and the mean of all
#'   voxels inside each label is computed.
#'   }
#'   \item{sums -}{ Each file contains an aribtrary number and the sum of all
#'   voxels inside each label is computed.
#'   }
#'   \item{text -}{ Each file is a comma separated values text file and is simply
#'   read in.
#'   }
#' }
#' @return A matrix with ncols equal to the number of labels in the atlas and
#' nrows equal to the number of files.
#' @export
anatGetAll <-
  function(filenames, atlas = NULL, 
           defs = getOption("RMINC_LABEL_DEFINITIONS"), 
           method = c("jacobians", "labels", "sums", "means", "text"), 
           side = c("both", "left", "right"),
           parallel = NULL, 
           strict = TRUE
         , conf_file = getOption("RMINC_BATCH_CONF")){
    
    method <- match.arg(method)
    
    if(method != "labels" && is.null(atlas))
      stop("An atlas is required for all methods other than \"labels\" ")
      
    mincFileCheck(filenames)
    
    ## Handle dispatch of worker functions
    compute_summary <-
      function(filenames){
        switch(method
               , labels =
                   return(
                     label_counts(filenames = filenames
                                , defs = defs
                                , side = side
                                , strict = strict))
               , test = 
                   return(
                     Reduce(function(acc, file){ cbind(acc, read.csv(file, header = FALSE)) }
                            , filenames
                            , NULL))
               , #default
                   return(
                     anatomy_summary(filenames = filenames
                                     , atlas = atlas
                                     , method = method
                                     , defs = defs
                                     , side = side
                                     , strict= strict)))
        }
    
    ## Special matrix reducer to ensure rows are bound by index
    reduce_matrices <- 
      function(mat_list){
        mat_list %>%
          lapply(function(m) tibble::rownames_to_column(as.data.frame(m))) %>%
          Reduce(function(df, acc) full_join(df, acc, by = "rowname")
                 , ., data_frame(rowname = character())) %>%
          select_(~ -rowname) %>%
          as.matrix
      }
    
    ## Handle parallelism choices
    if (is.null(parallel)) {
      out <- compute_summary(filenames)
    }
    else {
      n_groups <- as.numeric(parallel[2])
      groups <- split(seq_along(filenames), groupingVector(length(filenames), n_groups))
      # a vector with two elements: the methods followed by the # of workers
      if (parallel[1] %in% c("local", "snowfall")) {
        out <- 
          failing_mclapply(groups, function(group){
            compute_summary(filenames[group])
          }, mc.cores = n_groups) %>%
          reduce_matrices %>%
          `colnames<-`(filenames) %>%
          setNA(0)
      }
      else {
        reg <- qMincRegistry(new_file("anatGet_registry")
                           , conf_file = conf_file)
        on.exit( tenacious_remove_registry(reg) )
        
        suppressWarnings( #Warning suppression for large env for function (>10mb)
          ids <- batchMap(reg = reg, function(group){
            compute_summary(filenames[group])
          }, group = groups)
        )
        
        submitJobs(ids, reg = reg)
        waitForJobs(reg = reg)
        
        out <-
          reduceResultsList(reg = reg) %>%
          reduce_matrices %>%
          `colnames<-`(filenames) %>%
          setNA(0)
      } 
    }
    
    missing_labels <- abs(rowSums(out)) == 0
    
    ## Handle creating the label frame
    if(!is.null(defs)){
      label_frame <- create_labels_frame(defs, side = side)
    } else if(!is.null(atlas)) {
      warning("No definitions provided, using indices from atlas \n",
              "to set a default label set options(RMINC_LABEL_DEFINITIONS),",
              " or set it as an environment variable")
      atlas_vol <- mincGetVolume(atlas) %>% round %>% as.integer
      uniq_labels <- unique(atlas_vol)
      label_frame <-
        data_frame(Structure = as.character(uniq_labels)
                   , label   = uniq_labels) 
    } else {
      warning("No definitions provided, using observed labels \n",
              "to set a default label set options(RMINC_LABEL_DEFINITIONS),",
              " or set it as an environment variable")
      uniq_labels <- which(!missing_labels)
      label_frame <-
        data_frame(Structure = as.character(uniq_labels)
                   , label   = uniq_labels) 
    }
    
    ## Dress up the results and do label error checking
    create_anat_results(out, label_frame, filenames = filenames, atlas = atlas)
  }

# Convert an RMINC style atlas label set to a nice data.frame
create_labels_frame <-
  function(defs, side = c("both", "left", "right"), hierarchy = NULL){
    side <- match.arg(side)
    
    # if the definitions are given, check to see that we can read the 
    # specified file
    if(file.access(as.character(defs), 4) == -1){
      stop("The specified label definitions can not be read: "
           , defs
           , "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
    }
    
    label_defs <- 
      read.csv(defs, stringsAsFactors = FALSE)
    
    if(!is.null(hierarchy) && hierarchy %in% names(label_defs)){
      label_defs <- 
        select(label_defs
               , !!! c("Structure"
                       , "left.label"
                       , "right.label"
                       , hierarchy)) %>% ##new tidyeval
        rename(hierarchy = !!hierarchy)
    } else {
      label_defs <- select(label_defs, !!! c("Structure", "left.label", "right.label"))
    }
    
    label_defs <- 
      label_defs %>%
      mutate(both_sides = .data$right.label == .data$left.label) %>%
      gather("hemisphere", "label", c("right.label", "left.label")) %>%
      mutate_(Structure =
                ~ ifelse(both_sides
                         , Structure
                         , ifelse(hemisphere == "right.label"
                                  , paste0("right ", Structure)
                                  , paste0("left ", Structure))))
    
    if("hierarchy" %in% names(label_defs))
      label_defs <-
      label_defs %>%
      mutate_(hierarchy = 
                ~ with(.
                       , case_when(is.na(hierarchy) | hierarchy == "" ~ ""
                                   , both_sides ~ hierarchy
                                   , hemisphere == "right.label" ~ paste0("right ", hierarchy)
                                   , hemisphere == "left.label" ~ paste0("left ", hierarchy))))
    
    label_defs <-
      switch(side
             , left  = filter_(label_defs, ~ both_sides | hemisphere == "left.label")
             , right = filter_(label_defs, ~ both_sides | hemisphere == "right.label")
             , both  = label_defs)
    
    label_defs <-
      label_defs %>%
      select_(~ -c(hemisphere, both_sides)) %>%
      filter_(~ ! duplicated(label))
    
    label_defs
  }

# create the results matrix from the c++ results and a label frame 
create_anat_results <-
  function(results, label_frame, filenames = NULL, atlas = NULL){
    
    missing_labels <- abs(rowSums(results)) == 0
    if(nrow(results) == 1) stop("No structures found, check your atlas or label volumes")
    
    results <- 
      as.data.frame(results) %>%
      slice(-1) %>% #removes the zero label
      mutate(indices = 1:n()) %>%
      filter_(~ ! missing_labels[-1]) %>%
      left_join(label_frame, by = c("indices" = "label"))
    
    extra_structures <- results %>% filter_(~ is.na(Structure)) %>% .$indices
    if(length(extra_structures) != 0)
      message("Extra Structures  found in files but not in labels: "
              , paste0(extra_structures, collapse = ", "))
    
    missing_structures <- 
      with(label_frame, Structure[! label %in% results$indices])
    if(length(missing_structures) != 0)
      message("Missing Structures found in labels but not in any files: "
              , paste0(missing_structures, collapse = ", "))
    
    results <- 
      results %>%
      filter_(~ !is.na(Structure))
    
    label_inds <- results$indices
    structures <- results$Structure
    
    results <-
      results %>%
      select_(~ -c(indices,Structure)) %>% 
      t %>%
      `colnames<-`(structures) %>%
      `rownames<-`(NULL)
    
    sapply(colnames(results), function(cn){
      missings <- which(results[,cn] == 0)
      
      if(length(missings) != 0)
        message("Files ", paste0(missings, sep = ", "), " are missing label ", cn)
    })
    
    attr(results, "anatIDs") <- label_inds
    if(!is.null(atlas))     attr(results, "atlas") <- atlas
    if(!is.null(filenames)) attr(results, "input") <- filenames
    class(results) <- c("anatUnilateral", "anatMatrix", "matrix")
    
    results
  }

#' @export
print.anatMatrix <- function(x, n = min(6, nrow(x)), width = min(6, ncol(x)), ...){
  cat("anatMatrix, showing ", n, " of ", nrow(x), " rows, and ", width, " of ", ncol(x), "columns\nprint more by supplying n and width to the print function\n")
  sub_x <- x[seq_len(n), seq_len(width)]
  
  print(sub_x)
  invisible(x)
}

#' @export
print.anatModel <- function(x, n = min(6, nrow(x)), width = min(6, ncol(x)), ...){
  cat("anatModel, showing ", n, " of ", nrow(x), " rows, and ", width, " of ", ncol(x), "columns\nprint more by supplying n and width to the print function\n")
  sub_x <- x[seq_len(n), seq_len(width)]
  attributes(sub_x) <- attributes(sub_x) %>% { .[names(.) %in% c("dim", "dimnames")]}
  
  print.default(sub_x)
  invisible(x)
}

#' @export
`[.anatModel` <- function(x, i, j, drop = TRUE){
  orig_x <- x
  mdrop <- missing(drop)
  n_args <- nargs() - !mdrop
  model_attrs <- attributes(x) %>% { .[! names(.) %in% c(attributes_to_skip, "dimnames", "class", "df", "stat-type")]}
  class(x) <- "matrix"
  
  if(n_args == 2){
    x <- x[i]
  } else {
    x <- x[i,j,drop = drop] 
  }
  
  if(!is.null(dim(x))){
    if(missing(j))
      j <- seq_len(ncol(x))
    
    if(is.character(j))
      j <- match(j, colnames(orig_x))
    
    if( is.logical(j))
      j <- which(j)
    
    class(x) <- class(orig_x)
    attributes(x) <- c(attributes(x), model_attrs)
    
    ## Deal with stats and dfs
    orig_stats <- attr(orig_x, "stat-type")
    orig_dfs <- attr(orig_x, "df")
    
    knownStats <- c("t", "F", "u", "chisq", "tlmer")
    if(length(orig_stats) == 1) orig_stats <- rep(orig_stats, ncol(orig_x))
    stat_columns <- which(orig_stats %in% knownStats)
    selected_stat_columns <- match(j, stat_columns, 0)
    
    if(length(orig_dfs) == 1) orig_dfs <- rep(orig_dfs, length(stat_columns))

    attr(x, "stat-type") <- orig_stats[j]
    attr(x, "df") <- orig_dfs[selected_stat_columns]
  }
  
  x
}

#' @export
rbind.anatModel <- function(..., deparse.level = 1){
  args <- list(...)
  attrs <- lapply(args, function(model) attributes(model) %>% { .[! names(.) %in% "call"]})
  lapply(attrs, function(att){
    if(! identical(att, attrs[[1]])) stop("Models aren't compatible, differ in some attributes")
  })
  
  bound_models <- do.call(rbind, lapply(args, `class<-`, "matrix"))
  final_attrs <- attrs[[1]] %>% { .[! names(.) %in% c("dimnames", "dim")]}
  attributes(bound_models) <- c(attributes(bound_models), final_attrs)
  rownames(bound_models) <- unlist(lapply(args, rownames))
  colnames(bound_models) <- colnames(args[[1]])

  bound_models
}

#' Create Anatomy data.frame
#' 
#' Convert an anatomy frame to data.frame for use with tidy tools
#' 
#' @param x An \code{anatMatrix} object produced by \link{anatGetAll}
#' @param ... additional parameters for methods
#' @return A data.frame augmented with addition attributes from anat
#' @details Note that structure names are munged such that each structure name has non-standard
#' characters replaced by underscores
#' @export
as.data.frame.anatMatrix <- 
  function(x, ...){
    minc_attrs <- mincAttributes(x)
    minc_attrs <- minc_attrs[! names(minc_attrs) %in% c("dimnames", "class") ]
    
    x %>%
      unclass %>%
      as.data.frame %>%
      setNames(fix_names(names(.))) %>%
      setMincAttributes(minc_attrs)
  }

anatomy_summary <- 
  function(filenames, atlas, 
           defs = getOption("RMINC_LABEL_DEFINITIONS"), 
           method = c("jacobians", "labels", "sums", "means"), 
           side = c("both", "left", "right"),
           strict = TRUE){
    
    method <- match.arg(method)
    side <- match.arg(side)
    atlas_vol <- mincGetVolume(atlas) %>% round %>% as.integer
    
    check_step_sizes(filenames, atlas, strict)

    cpp_res <- anat_summary(filenames, atlas_vol, method)
    cpp_res$counts <- cpp_res$counts
    cpp_res$values <- cpp_res$values
    
    results <-
      switch(method
             , jacobians = cpp_res$value
             , means     = cpp_res$value / cpp_res$counts
             , sums      = cpp_res$value)
    
    results[!is.finite(results)] <- 0
    
    results
  }

label_counts <- 
  function(filenames,
           defs = getOption("RMINC_LABEL_DEFINITIONS"), 
           method = c("jacobians", "labels", "sums", "means"), 
           side = c("both", "left", "right"),
           strict = TRUE){
    
    method <- match.arg(method)
    side <- match.arg(side)

    check_step_sizes(filenames, strict = strict)
    
    cpp_res <- count_labels(filenames)
    missing_labels <- abs(rowSums(cpp_res)) == 0
    results <- cpp_res
    
    results
  }

#' Summarize anatGetAll results by hierarchy
#' 
#' Take the results of \link{anatGetAll} and summarize a grouping label, typically the structures
#' group in the anatomical hierarchy. 
#' 
#' @param anat The results of a call to \link{anatGetAll}
#' @param summarize_by either a data frame with grouping information or a path to label definitions. 
#' If a data frame is passed it must have two columns: "label" containing the structure
#' name (see \code{colnames(anat)} to check the relevant labels) and "grouping" containg the name 
#' of the group each structure belongs to.
#' @param defs A text file containing the label definitions if \code{summarize_by} is a string
#' @param discard_missing logical controlling how to handle structures with no ("") group information.
#' If TRUE filter these structures, if FALSE give each structure a group label to match their structure
#' name.
#' @return A matrix with columns representing groups from the hierarchy and rows representing individuals
#' with values equal to the sum of the individual members of each group.
#' @export 
anatSummarize <-
  function(anat
           , summarize_by = "hierarchy"
           , defs = getOption("RMINC_LABEL_DEFINITIONS")
           , discard_missing = FALSE){
    
    if(is.character(summarize_by) && length(summarize_by == 1)){
      summarize_by <- 
        create_labels_frame(defs, hierarchy = summarize_by) %>%
        select_(~ -label) %>%
        rename_(label = "Structure", group = "hierarchy")
    }
    
    if(!discard_missing){
      summarize_by <-
        summarize_by %>%
        mutate_(group = ~ifelse(is.na(group) | group == "", label, group))
    } else {
      summarize_by <- summarize_by %>% filter_(~ group != "")
    }
    
    anat %>%
      as.data.frame.matrix %>%
      (tibble::rownames_to_column) %>%
      gather_("label", "value", setdiff(colnames(anat), "rowname")) %>%
      inner_join(summarize_by, by = "label") %>%
      group_by_("group", "rowname") %>%
      summarize_(value = ~ sum(value)) %>%
      spread_("group", "value") %>%
      arrange_(~ as.numeric(rowname)) %>%
      select_(~ -rowname) %>%
      as.matrix
  }


###########################################################################################
#' @description Computes volumes, means, sums, and similar values across a
#' segmented atlas
#' @title Get values given a set of files and an atlas
#' @param filenames A vector of filenames (strings) which contain the
#' information to be extracted at every structure in the atlas.
#' @param atlas A single filename containing the atlas definitions. This MINC
#' volume has to be of the same sampling (sizes and dimension
#' order) as the filenames specified in the first argument and
#' use a separate integer for each atlas label.
#' @param method A string specifying the way information is to be computed at
#' every voxel. See the details section for the possible options
#' and what they mean.
#' @param defs A string pointing to the filename containing the label
#' definitions. Used to map the integers in the atlas to a
#' proper name for the structure and contains additional
#' information for laterality of each structure. See \link{voxel_atlas_defs}
#' for details.
#' @param dropLabels Whether to return a value for every structure in the defs
#' or just for the ones actually contained in each file.
#' @param side Three choices, "right", "left", and "both" (the default)
#' which specify what labels to obtain.
#' @details 
#' anatGetAll needs a set of files along with an atlas and a set of
#' atlas definitions. In the end it will produce one value per label
#' in the atlas for each of the input files. How that value is
#' computed depends on the methods argument:
#' \itemize{
#'   \item{jacobians -}{ Each file contains log jacobians, and the volume for
#'   each atlas label is computed by multiplying the jacobian with
#'   the voxel volume at each voxel.
#'   }
#'   \item{labels -}{ Each file contains integer labels (i.e. same as the atlas).
#'   The volume is computed by counting the number of voxels with
#'   each label and multiplying by the voxel volume.
#'   }
#'   \item{means -}{ Each file contains an arbitrary number and the mean of all
#'   voxels inside each label is computed.
#'   }
#'   \item{sums -}{ Each file contains an aribtrary number and the sum of all
#'   voxels inside each label is computed.
#'   }
#'   \item{text -}{ Each file is a comma separated values text file and is simply
#'   read in.
#'   }
#' }
#' @return A matrix with ncols equal to the number of labels in the atlas and
#' nrows equal to the number of files.

#' @seealso anatLm,anatCombineStructures
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' filenames <- read.csv("/tmp/rminctestdata/filenames.csv")
#' volumes <- anatGetAll(filenames=filenames$absolute_jacobian, 
#'                       atlas="/tmp/rminctestdata/test_segmentation.mnc", 
#'                       method="jacobians",
#'                       defs="/tmp/rminctestdata/test_defs.csv")
#'}
#'@export
anatGetAll_old <- function(filenames, atlas, method="jacobians", 
                       defs = getOption("RMINC_LABEL_DEFINITIONS"), 
                       dropLabels = TRUE, 
                       side="both") {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  # Get output dimensions from full set of label definitions
  labeldefs <- read.csv(defs) 
  labels <- c(labeldefs$right.label, labeldefs$left.label)
  labels.sorted <- sort(labels)
  usedlabels <- labels.sorted[!duplicated(labels.sorted)]
  
  output <- matrix(nrow = length(usedlabels), ncol = length(filenames))
  rownames(output) <- usedlabels

  for (i in 1:length(filenames)){
    vol <- anatGetFile(filenames[i], atlas, method, defs, dropLabels, side)

    ### some quick error checking ###
    # get the list of all the labels that was found in the current file
    labels_found <- vol$V1
    # find labels that occurred in the file, but not in the mapping
    diff_between_file_and_mapping <- setdiff(labels_found, usedlabels)
    if(length(diff_between_file_and_mapping) > 0) {
      print(paste("WARNING: Labels found in the inputfile, but not in the mapping: ", diff_between_file_and_mapping))
    }
    # find labels that occurred in the mapping, but not in the file
    diff_between_mapping_and_file <- setdiff(usedlabels,labels_found)
    if(length(diff_between_mapping_and_file) > 0) {
      print(paste("WARNING: Labels found in the mapping, but not in the file: ", diff_between_mapping_and_file))
    }
    ### end of error checking ###

    for (j in 1:length(usedlabels)){
       if(usedlabels[j] == vol[j,1]){
          output[j,i] <- vol[j,2]
       }
       else {
          #If labels don't exist for a particular volume, set their value to zero
          toadd <- list(usedlabels[j],0)
          dim = nrow(vol)
          vol <- rbind(vol[1:j-1,], toadd, vol[j:dim,]) 
          output[j,i] <- 0
       }			
    }
  }
 
  attr(output, "atlas") <- atlas
  attr(output, "input") <- filenames
  if (! is.null(defs)) {
    output <- anatRenameRows(output, defs)
  }
  return(t(output))
}
###########################################################################################
#' @description Combines left and right labels from volumes obtained from anatGetAll call
#' @name anatCombineStructures
#' @title Combine left and right volumes
#' @param vols Matrix output from call to anatGetAll
#' @param method A string specifying the way information was computed at
#' every voxel ("jacobians","labels","means","sums")
#' @param defs A string pointing to the filename containing the label
#' definitions. Used to map the integers in the atlas to a
#' proper name for the structure and contains additional
#' information for laterality of each structure.
#' @details anatCombineStructures collapses left and right volume information into one
#' measure. If "jacobians","sums",or "labels" is selected then the sum of the left and right is produced, otherwise
#' the mean is produced.
#' @return A matrix with ncols equal to the number of collapsed labels

#' @seealso anatLm,anatGetAll
#' @examples
#' \dontrun{
#' getRMINCTestData()
#' filenames <- read.csv("/tmp/rminctestdata/filenames.csv")
#' volumes <- anatGetAll(filenames=filenames$absolute_jacobian, 
#'                       atlas="/tmp/rminctestdata/test_segmentation.mnc", 
#'                       method="jacobians",
#'                       defs="/tmp/rminctestdata/test_defs.csv")
#' volumes_combined <- 
#'      anatCombineStructures(vols=volumes, 
#'                            method="jacobians",
#'                            defs="/tmp/rminctestdata/test_defs.csv")
#' }
#' @export
anatCombineStructures <- function(vols, method = "jacobians", 
                                  defs = getOption("RMINC_LABEL_DEFINITIONS")) {
  if(defs == ""){
    stop("No label definitions specified. Either use the defs argument, or use the environment variable $RMINC_LABEL_DEFINITIONS.")    
  }
  # if the definitions are given, check to see that we can read the 
  # specified file
  if(file.access(as.character(defs), 4) == -1){
    stop("The specified label definitions can not be read: ", defs, "\nUse the defs argument or the $RMINC_LABEL_DEFINITIONS variable to change.")
  }
  labels <- read.csv(defs)
  combined.labels <- matrix(nrow=nrow(vols), ncol=nrow(labels))
  labelNumbers <- attr(vols, "anatIDs")
  for (i in 1:nrow(labels)) {
    if (labels$right.label[i] == labels$left.label[i]) {
      combined.labels[,i] <- vols[, labelNumbers == as.character(labels$right.label[i])]
    }
    else {
      if (method == "jacobians" || method == "labels" || method == "sums") {
        combined.labels[,i] <-
          vols[, labelNumbers == as.character(labels$right.label[i])] + vols[, labelNumbers == as.character(labels$left.label[i])]
      }
      else if (method == "means"){
        combined.labels[,i] <-
          (vols[, labelNumbers == as.character(labels$right.label[i])] + vols[, labelNumbers == as.character(labels$left.label[i])]) /2
      }
        
    }
  }
  colnames(combined.labels) <- labels$Structure
  class(combined.labels) <- c("anatCombined", "anatMatrix", "matrix")
  attr(combined.labels, "atlas") <- attr(vols, "atlas")
  attr(combined.labels, "definitions") <- defs
  return(combined.labels)
}
###########################################################################################
#' @description This function is used to compute an arbitrary function of every region in an anat structure.
#' @name anatApply
#' @title Apply function over anat structure
#' @param vols anatomy volumes
#' @param grouping grouping with which to perform operations
#' @param method The function which to apply [default mean]
#' @return  out: The output will be a single vector containing as many
#'          elements as there are regions in the input variable by the number of groupings
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- anatApply(gf$lobeThickness,gf$Primary.Diagnosis)
#' }
#' @export
anatApply <- function(vols, grouping, method=mean) {
  ngroups <- length(levels(grouping))
  output <- matrix(nrow=ncol(vols), ncol=ngroups)

  for (i in 1:ncol(vols)) {
    output[i,] <- tapply(vols[,i], grouping, method)
  }
  colnames(output) <- levels(grouping)
  rownames(output) <- colnames(vols)
  return(output)
}
###########################################################################################
#' Calculates statistics and coefficients for linear model of specified anat structure
#' @param formula a model formula. The RHS of the formula may contain one term with a matrix. If
#' so only the + operator may be used, and only two terms may appear on the RHS
#' @param data a data.frame containing variables in formula 
#' @param anat an array of atlas labels vs subject data
#' @param subset rows to be used, by default all are used
#' @param weights An optional set of weights to use the regression, must be one per subject
#' @return Returns an object containing the R-Squared,value,coefficients,F 
#' and t statistcs that can be passed directly into anatFDR. Additionally
#' has the attributes for model,stat type and degrees of freedom.
#' @seealso mincLm,anatLm,anatFDR 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' rmincLm= anatLm(~ Sex,gf,gf$lobeThickness)
#' } 
#' @export
anatLm <- function(formula, data, anat, subset=NULL, weights = NULL) {
  
  #INITIALIZATION
  matrixFound = FALSE
  mmatrix =  matrix()
  matrixName = '';

  # Build model.frame
  m <- m_orig <- match.call()
  mf <- match.call(expand.dots=FALSE)
  # Add a row to keep track of what subsetting does
  mf$rowcount <- seq(1, nrow(data))
  m <- match(c("formula", "data", "subset", "rowcount"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

   if(length(grep("\\$",formula[[2]])) > 0) {
	stop("$ Not Permitted in Formula")  
  }

  # note - the formula parsing appears to be implemented again 
  # in parseLmFormula; should reconcile.
  # Only 1 Term on the RHS
  if(length(formula[[2]]) == 1) {
    
    term <- data[[as.character(formula[[2]])]]
    
    if (is.matrix(term)) {
      # Save term name for later
      rows = c('Intercept',formula[[2]])
      matrixName = formula[[2]]
      matrixFound = TRUE
      data.matrix.left <- t(anat[mf[,"(rowcount)"],])
      #data.matrix.left <- vertexTable(filenames)
      data.matrix.right <- t(term)
    }  
  }
  # Multiple Terms on RHS
  else {
    for (nTerm in 2:length(formula[[2]])){
      rCommand = paste("term <- data$",formula[[2]][[nTerm]],sep="")
      if(!all(as.character(formula[[2]][[nTerm]]) %in% names(data))) {
        next
      }
      eval(parse(text=rCommand))	
      if (is.matrix(term)) {
        matrixName = formula[[2]][[nTerm]]
        matrixFound = TRUE
        #data.matrix.left <- vertexTable(filenames)
        data.matrix.left <- t(anat[mf[,"(rowcount)"],])
        data.matrix.right <- t(term)
        
      }  else  {
        tmpFormula = formula
        rCommand = paste("formula <-",formula[[1]],"~",formula[[2]][[nTerm]],sep="")
        eval(parse(text=rCommand))	
        mmatrix <- model.matrix(formula, mf)	
        formula = tmpFormula	
      }
    }
    rows = colnames(mmatrix)
    rows = append(rows,matrixName)
  }	
  
  
  
  # Call subroutine based on whether matrix was found
  if(!matrixFound) {    	
    mmatrix <- model.matrix(formula, mf)	
    data.matrix.left <- t(anat[mf[,"(rowcount)"],])
    rows = colnames(mmatrix) 
    data.matrix.right = matrix()
  }
  
  if(!is.null(weights)){
    if(length(weights) != nrow(mmatrix)) stop("Incorrect number of weights supplied, need one per observation")
    result <- .Call("vertex_wlm_loop",data.matrix.left,data.matrix.right,mmatrix,weights,PACKAGE="RMINC") 
  } else {
    result <- .Call("vertex_lm_loop",data.matrix.left,data.matrix.right,mmatrix,PACKAGE="RMINC") 
  }
  
  rownames(result) <- colnames(anat)
  
  # the order of return values is:
  #
  # f-statistic
  # r-squared
  # betas
  # t-stats
  #
  betaNames                   <- paste('beta-',rows, sep='')
  tnames                      <- paste('tvalue-',rows, sep='')
  colnames(result)            <- c("F-statistic", "R-squared", betaNames, tnames, "logLik")
  class(result)               <- c("anatModel", "matrix")
  attr(result, "atlas")       <- attr(anat, "atlas")
  attr(result, "definitions") <- attr(anat, "definitions")
  attr(result, "model")       <- as.matrix(mmatrix)
  attr(result, "data")        <- data 
  attr(result, "call")        <- m_orig
  attr(result, "stat-type")   <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2), "logLik")
  
  Fdf1 <- ncol(attr(result, "model")) -1
  Fdf2 <- nrow(attr(result, "model")) - ncol(attr(result, "model"))

  # degrees of freedom are needed for the fstat and tstats only
  dflist <- vector("list", (ncol(result)-2)/2 + 1)
  dflist[[1]] <- c(Fdf1, Fdf2)
  dflist[2:length(dflist)] <- Fdf2
  attr(result, "df") <- dflist
  
  # run the garbage collector...
  gcout <- gc()
  
  return(result)
}

#' Anatomy Linear Mixed Effects Model
#' 
#' Fit a linear mixed effects model for each structure in the results of
#' \link{anatGetAll}.
#' 
#' @param anat a subject by label matrix of anatomical 
#' summaries typically produced by \link{anatGetAll}
#' @inheritParams mincLmer
#' @details \code{anatLmer}, like its relative \link{mincLmer} provides an interface to running 
#' linear mixed effects models at every vertex. Unlike standard linear models testing hypotheses 
#' in linear mixed effects models is more difficult, since the denominator degrees of freedom are 
#' more difficult to  determine. RMINC provides estimating degrees of freedom using the
#' \code{\link{anatLmerEstimateDF}} function. For the most likely models - longitudinal
#' models with a separate intercept or separate intercept and slope per subject - this
#' approximation is likely correct. Be careful in using this approximations if
#' using more complicated random effects structures.
#' @export 
anatLmer <-
  function(formula, data, anat, REML = TRUE
           , control = lmerControl(), verbose = FALSE
           , start = NULL, parallel = NULL, safely = FALSE
           , summary_type = "fixef"
           , weights = NULL){
    
    mc <- mcout <- match.call()
    original_data <- data
    
    # Add a dummy response column
    mc$formula <- update(formula,  RMINC_DUMMY_LHS ~ .)
    data$RMINC_DUMMY_LHS <- rnorm(nrow(data))
    mc[["data"]] <- as.symbol("data")
    
    mc[[1]] <- quote(lme4::lFormula)
    
    # remove lme4 unknown arguments, since lmer does not know about them and keeping them
    # generates obscure warning messages
    mc <- mc[!names(mc) %in% c("anat", "subset", "parallel", "safely", "summary_type")]
    lmod <- eval(mc, environment())
    
    # code ripped from lme4:::mkLmerDevFun
    rho <- new.env(parent = parent.env(environment()))
    rho$pp <- do.call(merPredD$new, c(lmod$reTrms[c("Zt", "theta", 
                                                    "Lambdat", "Lind")],
                                      n = nrow(lmod$X), list(X = lmod$X)))
    REMLpass <- if (REML) 
      ncol(lmod$X)
    else 0L
    
    
    mincLmerList <- list(lmod, mcout, control, start, verbose, rho, REMLpass)
    
    summary_fun <- summary_type
    if(is.character(summary_type) && length(summary_type) == 1)
      summary_fun <- switch(summary_type
                          , fixef = fixef_summary
                          , ranef = ranef_summary
                          , both = effect_summary
                          , anova = anova_summary
                          , stop("invalid summary type specified"))

    if(!is.function(summary_fun))
      stop("summary_type must be a string specifying a summary, or a function")
    
    mincLmerOptimizeAndExtractSafely <-
      function(x, mincLmerList, summary_fun){
        tryCatch(mincLmerOptimizeAndExtract(x, mincLmerList, summary_fun),
                 error = function(e){warning(e); return(NA)})
      }
    
    optimizer_fun <-
      `if`(safely, mincLmerOptimizeAndExtractSafely, mincLmerOptimizeAndExtract)
    
    out <- 
      matrixApply(t(anat)
                  , optimizer_fun
                  , mincLmerList = mincLmerList
                  , summary_fun = summary_fun
                  , parallel = parallel)
    
    out[is.infinite(out)] <- 0            #zero out infinite values produced by vcov
    
    #termnames <- colnames(lmod$X)
    #betaNames <- paste("beta-", termnames, sep="")
    #tnames <- paste("tvalue-", termnames, sep="")
    #colnames(out) <- c(betaNames, tnames, "logLik", "converged")
    rownames(out) <- colnames(anat)
    
    # generate some random numbers for a single fit in order to extract some extra info
    mmod <- mincLmerOptimize(rnorm(length(lmod$fr[,1])), mincLmerList)
    
    res_cols <- colnames(out)
    attr(out, "stat-type") <- ## Handle all possible output types
      check_stat_type(res_cols, summary_type)
    
    # get the DF for future logLik ratio tests; code from lme4:::npar.merMod
    attr(out, "logLikDF") <- length(mmod@beta) + length(mmod@theta) + mmod@devcomp[["dims"]][["useSc"]]
    attr(out, "REML") <- REML
    attr(out, "mincLmerList") <- mincLmerList
    attr(out, "atlas") <- attr(anat, "atlas")
    attr(out, "definitions") <- attr(anat, "definitions")
    attr(out, "data") <- data
    attr(out, "anat") <- anat
    
    class(out) <- c("anatLmer", "anatModel", "matrix")
    
    return(out)
  }

#' Estimate Degrees of freedom for an anatLmer model
#' 
#' Estimate the degrees of freedom for an \code{anatLmerModel} object
#' produced by \link{anatLmer}. See \link{anatLmer} for more details
#' on degrees of freedom estimation for linear mixed effects models
#' @param model an \code{anatLmerModel}
#' @param n number of structures to use for DF estimation
#' @param verbose Whether or not to print progress
#' @return the same model, now with degrees of freedom set
#' @export
anatLmerEstimateDF <- function(model, n=50, verbose = FALSE) {
  # set the DF based on the Satterthwaite approximation
    
  mincLmerList <- attr(model, "mincLmerList")
  anat <- attr(model, "anat")
  original_data <- attr(model, "data")

  ## It seems LmerTest cannot compute the deviance function for mincLmers
  ## in the current version, instead extract the model components from
  ## the mincLmerList and re-fit the lmers directly at each structure,
  ## Slower but yeilds the correct result
  lmod <- mincLmerList[[1]]
  
  # estimated DF depends on the input data. Rather than estimate separately at every structure,
  # instead select a small number of structures and estimate DF for those structures, then keep the
  # min
  nstructures <- min(n, nrow(model))
  rstructures <- sample(seq_len(nrow(model)), nstructures)
  dfs <- matrix(nrow = nstructures, ncol = sum(attr(model, "stat-type") %in% "tlmer"))

  for(i in seq_len(nstructures)){
    if(verbose) cat("Iteration ", i, " of ", n)
    
    rand_ind <- rstructures[i]
    structureData <- anat[,rand_ind]
    
    original_data$RMINC_DUMMY_LHS <- structureData

    ## Work around for slowness in recent lme4, fixed in upstream lme4
    ## thanks to https://github.com/lme4/lme4/issues/410#issuecomment-311092416
    model_env <- list2env(original_data)
    form <- lmod$formula
    environment(form) <- model_env

    mmod <-
      lmerTest::lmer(form, REML = lmod$REML,
                     start = mincLmerList[[4]], control = mincLmerList[[3]],
                     verbose = mincLmerList[[5]])
    
    dfs[i,] <- 
      suppressMessages(
        tryCatch(lmerTest::summary(mmod)$coefficients[,"df"]
               , error = function(e){ 
                 warning("Unable to estimate DFs for structure "
                       , colnames(anat)[rand_ind]
                       , call. = FALSE)
                 NA}))
  }
  
  df <- apply(dfs, 2, median, na.rm = TRUE)
  cat("Mean df: ", apply(dfs, 2, mean, na.rm = TRUE), "\n")
  cat("Median df: ", apply(dfs, 2, median, na.rm = TRUE), "\n")
  cat("Min df: ", apply(dfs, 2, min, na.rm = TRUE), "\n")
  cat("Max df: ", apply(dfs, 2, max, na.rm = TRUE), "\n")
  cat("Sd df: ", apply(dfs, 2, sd, na.rm = TRUE), "\n")
  
  attr(model, "df") <- df
  
  return(model)
}

#' Performs ANOVA on each region specified 
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param anat  an array of atlas labels vs subject data
#' @param subset rows to be used, by default all are used
#' @return Returns an array with the F-statistic for each model specified 
#' by formula with the following attributes
#' \itemize{
#' \item{model}{ design matrix}
#' \item{stat-type}{ type of statistic used}
#' \item{df}{ degrees of freedom of each statistic}
#' } 
#' @seealso mincAnova,vertexAnova 
#' @examples
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET",
#'                            TRUE,"1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' rmincAnova = anatAnova(~ Sex,gf,gf$lobeThickness);
#' } 
#' @export
anatAnova <- function(formula, data=NULL, anat=NULL, subset=NULL) {
  # Create Model
  m  <- match.call()
  mf <- match.call(expand.dots=FALSE)
  # in order to keep track of which subjects we need (when subsetting),
  # we will store the row numbers
  mf$rowcount <- seq(1, nrow(data))
  m  <- match(c("formula", "data", "subset", "rowcount"), names(mf), 0)  
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mmatrix <- model.matrix(formula, mf)
  
  # Get the data from the anatomy matrix using the same
  # subset as specified in the formula (using the rownumbers)
  anatmatrix <- t(anat[mf[,"(rowcount)"],])
  
  # same stats as used for the vertices
  result <- .Call("vertex_anova_loop", anatmatrix, mmatrix,attr(mmatrix, "assign"), PACKAGE="RMINC");
  
  rownames(result) <- colnames(anat)
  # unlike in the anatLm function, the column names here should be the terms of the model
  colnames(result) <- attr(terms(formula), "term.labels")
  class(result) <- c("anatModel", "matrix")
  attr(result, "atlas") <- attr(anat, "atlas")
  attr(result, "definitions") <- attr(anat, "definitions")
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "stat-type") <-  rep("F", ncol(result))
  
  # get the structure in order to get the degrees of freedom
  firstStructure <- anatmatrix[1,] 
  l <- lm.fit(mmatrix, firstStructure)
  asgn <- l$assign[l$qr$pivot]
  dfr <- df.residual(l)
  dfs <- c(unlist(lapply(split(asgn, asgn), length)))
  dflist <- vector("list", ncol(result))
  for (i in 1:ncol(result)) {
    dflist[[i]] <- c(dfs[[i+1]], dfr)
  }
  attr(result, "df") <- dflist
  
  return(result)
}

#' Create a Statistic Volume
#' 
#' Convert an anatomy volume into a statistic volume
#' with the results of an \link{anatLm}
#' 
#' @param anat A minc volume object with a \code{labels} attribute
#' @param filename the filename for the new minc volume
#' @param column which column of the \link{anatLm} results to use
#' @return Invisibly returns the created volume
#' @export
anatCreateVolume <- function(anat, filename, column=1) {
  labels <- read.csv(attr(anat, "definitions"))
  volume <- mincGetVolume(attr(anat, "atlas"))
  newvolume <- volume
  for (i in 1:nrow(labels)) {
    newvolume[volume < labels$right.label[i] + 0.5 &
              volume > labels$right.label[i] - 0.5] <-
                anat[labels$Structure[i], column]
    newvolume[volume < labels$left.label[i] + 0.5 &
              volume > labels$left.label[i] - 0.5] <-
                anat[labels$Structure[i], column]
  }
  mincWriteVolume(newvolume, filename, attr(anat, "atlas"))
  
  return(invisible(newvolume))
}

#' anatSummaries
#' 
#' These functions are used to compute the mean, standard deviation,
#'    		sum, or variance of every region in an anat structure.
#' @param anat anat structure.
#' @return  out: The output will be a single vector containing as many
#'          elements as there are regions in the input variable. 
#' @examples 
#' \dontrun{
#' getRMINCTestData() 
#' gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
#' gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
#' gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
#' vm <- anatMean(gf$lobeThickness)
#' }
#' @name anatSummaries
NULL

#' @describeIn anatSummaries mean
#' @export 
anatMean <- function(anat) {
   return(rowMeans(t(anat)))
}

#' @describeIn anatSummaries sum
#' @export
anatSum <- function(anat) {
   return(rowSums(t(anat)))
}

#' @describeIn anatSummaries variance
#' @export
anatVar <- function(anat) {
   return(apply(t(anat),1,var))
}

#' @describeIn anatSummaries standard deviation
#' @export
anatSd <- function(anat) {
   return(apply(t(anat),1,sd))
}
