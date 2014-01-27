# compute a linear model over every structure


anatGetFile <- function(filename, atlas, method="jacobians", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv", dropLabels=FALSE, side="both" ) {
  out <- NULL
  if (method == "jacobians") {
    system(paste("label_volumes_from_jacobians", atlas, filename, "> tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "labels") {
    # filename here should be a set of labels unique to this brain
    system(paste("volumes_from_labels_only.py", filename, "tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "means") {
    system(paste("compute_values_across_segmentation", "-m",
                 filename, atlas, "tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "sums") {
    system(paste("compute_counts_for_labels",
                  atlas, filename, "> tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "slow_sums") {
    system(paste("compute_values_across_segmentation", "-s",
                 filename, atlas, "tmp.txt", sep=" "))
    out <- read.csv("tmp.txt", header=FALSE)
  }
  else if (method == "text") {
    # values are already extracted and stored in a text file
    out <- read.table(filename, header=FALSE)
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
  return(out)
}

anatRenameRows <- function(anat, defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv") {
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

# both unilateral and bilateral matrices can be printed the same way
print.anatMatrix <- function(x, ...) {
  print.table(x)
}

anatGetAll <- function(filenames, atlas, method="jacobians", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv", dropLabels=TRUE, side="both") {
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
  if (! is.null(defs)) {
    output <- anatRenameRows(output, defs)
  }
  return(t(output))
}

anatCombineStructures <- function(vols, method="jacobians", defs="/projects/mice/jlerch/cortex-label/c57_brain_atlas_labels.csv") {
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
#' @param formula a model formula
#' @param data a data.frame containing variables in formula 
#' @param anat an array of atlas labels vs subject data
#' @param subset rows to be used, by default all are used
#' @return Returns an object containing the coefficients,F 
#' and t statistcs that can be passed directly into anatFDR.
#' @seealso mincLm,anatLm,anatFDR 
#' @examples 
#' gf = read.csv("~/SubjectTable.csv") 
#' civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
#' gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
#' result = anatLm(~Primary.Diagnosis,gf,gf$lobeVolume)
#' anatFDR(result)

###########################################################################################  
anatLm <- function(formula, data, anat, subset=NULL) {
  
  #INITIALIZATION
  matrixFound = FALSE
  mmatrix =  matrix()

  # Build model.frame
  m <- match.call()
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


  # Only 1 Term on the RHS
  if(length(formula[[2]]) == 1) {
	  rCommand = paste("term <- data$",formula[[2]],sep="")
	  eval(parse(text=rCommand))

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
		  eval(parse(text=rCommand))	
		  if (is.matrix(term)) {
                        matrixName = formula[[2]][[nTerm]]
			matrixFound = TRUE
			#data.matrix.left <- vertexTable(filenames)
			data.matrix.left <- t(anat[mf[,"(rowcount)"],])
			data.matrix.right <- t(term)

			}
		  else  {
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
  if(matrixFound) {
	   result <- .Call("vertex_lm_loop_file",data.matrix.left,data.matrix.right,mmatrix,PACKAGE="RMINC") 
	}
  else  {      	
	mmatrix <- model.matrix(formula, mf)	
	data.matrix.left <- t(anat[mf[,"(rowcount)"],])
	rows = colnames(mmatrix)
	result <- .Call("vertex_lm_loop", data.matrix.left, mmatrix, PACKAGE="RMINC");		 
        }

  rownames(result) <- colnames(anat)

  # the order of return values is:
  #
  # f-statistic
  # r-squared
  # betas
  # t-stats
  #
  betaNames = paste('beta-',rows, sep='')
  tnames = paste('tvalue-',rows, sep='')
  colnames(result) <- c("F-statistic", "R-squared", betaNames, tnames)
  class(result) <- c("anatModel", "matrix")
  attr(result, "atlas") <- attr(anat, "atlas")
  attr(result, "definitions") <- attr(anat, "definitions")
  attr(result, "model") <- as.matrix(mmatrix)
  attr(result, "stat-type") <- c("F", "R-squared", rep("beta",(ncol(result)-2)/2), rep("t",(ncol(result)-2)/2))
  
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
}

anatFDR <- function(buffer, method="FDR") {
  vertexFDR(buffer, method)
}
