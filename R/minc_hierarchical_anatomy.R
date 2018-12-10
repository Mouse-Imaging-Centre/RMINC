#' Converts a column from a hierarchical tree into a volume for viewing or saving
#'
#' @param anatTree The anatomical tree
#' @param labelVolume The volume containing the label definitions
#' @param column String indicating which column to turn into a volume
#'
#' @return The volume with the values from the anatomical tree
#' @export
#'
#' @examples \dontrun{
#' labelVol <- mincArray(mincGetVolume("some-labels.mnc"))
#' hLm <- hanatLm(~Sex, gfBasic, vols)
#' statsvol <- hanatToVolume(hLm, labelVol, "F.statistic")
#' mincPlotSliceSeries(anatVol, statsvol, anatLow = 700, anatHigh = 1400, 
#'   low=1, high=10, symmetric = F, begin=50, end=-50)
#' }
hanatToVolume <- function(anatTree, labelVolume, column) {
  out <- array(0, dim=dim(labelVolume))
  labels <- c()
  values <- c()
  nodeList <- Traverse(anatTree, filterFun = isLeaf)
  for (i in 1:length(nodeList)) {
    nlabels <- length(nodeList[[i]]$label)
    labels <- c(labels, nodeList[[i]]$label)
    values <- c(values, rep(nodeList[[i]][[column]], nlabels))
  }
  if (is.numeric(values)) {
    out <- array(0, dim=dim(labelVolume))
    replaceValues(labelVolume, out, labels, values)
  }
  else if (is.character(values)) {
    labels <- c(labels, 0)
    values <- c(values, NA)
    out <- array("", dim=dim(labelVolume))
    replaceColours(labelVolume, out, labels, values)
  }
  return(out)
}

#' linear model applied to every item in an anatomical tree
#' 
#' See \link{anatLm} for more detail
#' 
#' @param formula The model formula
#' @param data The data frame containing the right side of the formula
#' @param anatTree The anatomical tree
#'   
#' @return The anatomical tree with the model information added
#' @export
#' 
#' @details The volumes inside the anatomical hierarchy and the data must be of
#'   the same length and ordering.
#'   
#' @examples \dontrun{
#' vols <- addVolumesToHierarchy(hdefs, allvols)
#' hLm <- hanatLm(~Sex, gf, vols)
#' }
hanatLm <- function(formula, data, anatTree) {
  # create a copy of the anatTree
  modelTree <- Clone(anatTree)
  # get a simple matrix of volumes
  volTable <- modelTree$Get("volumes")
  # run standard anatLm
  modelled <- anatLm(formula, data, volTable)
  # and return the model table back to the tree
  modelToTree(modelled, modelTree)
  
  return(modelTree)
}

#' ANOVA applied to every item in an anatomical tree
#' 
#' See \link{anatAnova} for more detail
#' 
#' @param formula The model formula
#' @param data The data frame containing the right side of the formula
#' @param anatTree The anatomical tree
#'   
#' @return The anatomical tree with the model information added
#' @export
#' 
#' @details The volumes inside the anatomical hierarchy and the data must be of
#'   the same length and ordering.
#'   
#' @examples \dontrun{
#' vols <- addVolumesToHierarchy(hdefs, allvols)
#' hLm <- hanatAnova(~Sex, gf, vols)
#' }
hanatAnova <- function(formula, data, anatTree) {
  # create a copy of the anatTree
  modelTree <- Clone(anatTree)
  # get a simple matrix of volumes
  volTable <- modelTree$Get("volumes")
  # run standard anatAnova
  modelled <- anatAnova(formula, data, volTable)
  # and return the model table back to the tree
  modelToTree(modelled, modelTree)
  
  return(modelTree)
}

#' linear mixed effects model applied to every item in an anatomical tree
#' 
#' See \link{anatLmer} for more detail
#'
#' @param formula The model formula
#' @param data The data frame containing the right side of the formula
#' @param anatTree The anatomical tree
#' @param ... Extra options (REML, control of parallel execution, etc.) passed on to anatLmer
#'
#' @return The anatomical tree with the model information added
#' @export
#' 
#' @details The volumes inside the anatomical hierarchy and the data must be of
#'   the same length and ordering.
#'
#' @examples \dontrun{
#' vols <- addVolumesToHierarchy(hdefs, allvols)
#' hLm <- hanatLmer(~Sex + (1|ID), gf, vols)
#' }
hanatLmer <- function(formula, data, anatTree, ...) {
  # create a copy of the anatTree
  modelTree <- Clone(anatTree)
  # get a simple matrix of volumes
  volTable <- modelTree$Get("volumes")
  # run standard anatLmer
  modelled <- anatLmer(formula, data, volTable, ...)
  # and return the model table back to the tree
  modelToTree(modelled, modelTree)
  
  return(modelTree)
}

modelToTree <- function(modelled, modelTree) {
  nodeList <- Traverse(modelTree)
  for (i in 1:length(nodeList)) { 
    for (n in colnames(modelled)) {
      nn <- gsub('-', '.', n) # the $ accessor does not play nicely with dashes
      nodeList[[i]][[nn]] <- modelled[i,n]
    }
    nodeList[[i]]$model <- modelled[i,]
  }
  # copy attributes
  modelTree <- copyExtraAttributes(modelled, modelTree)
  attr(modelTree, "modelClass") <- attr(modelled, "class")
  attr(modelTree, "modelDimnames") <- attr(modelled, "dimnames")
}

copyExtraAttributes <- function(source, target) {
  sourceAttributes <- attributes(source)
  for (n in names(sourceAttributes)) {
    if (n %in% c("dim", "dimnames", "class")) {
      attr(target, paste0("original", n)) <- attr(source, n)
    }
    else {
      attr(target, n) <- attr(source, n)
    }
  }
  return(target)
}

#' Compute the False Discovery Rate for an anatomical hierarchy
#' 
#' @param buffer What to compute FDR on
#'   
#' @return the anatomical hierarchy with FDR added (see details)
#' @export
#' 
#' @details Uses the Benjamini & Yekutieli (2001) algorithm for computing FDR to
#'   account for dependence (which is obviously present in a hierarchy where 
#'   parent values are directly derived from their children).
#'   
#'   Stores the values using the same name as the incoming statistic but with a
#'   q prepended. I.e. tvalue.Sexmale becomes qtvalue.Sexmale. Note that the original, 
#'   input, is modified with these new additions as well.
#' @examples \dontrun{
#' vols <- addVolumesToHierarchy(hdefs, allvols)
#' hLm <- hanatLm(~Sex, gfBasic, vols)
#' hLm <- hanatFDR(hLm)
#' thresholds(hLm)
#' }
hanatFDR <- function(buffer) {
  modelData <- buffer$Get("model")
  if (is.null(dim(modelData))) {
    dim(modelData) <- c(length(modelData), 1)
  }
  else {
    modelData <- t(modelData)
  }
  modelData <- copyExtraAttributes(buffer, modelData)
  dimnames(modelData) <- attr(buffer, "modelDimnames")
  class(modelData) <- attr(buffer, "modelClass")
  aFDR <- anatFDR(modelData, method="BY")
  
  nodeList <- Traverse(buffer)
  for (i in 1:length(nodeList)) { 
    for (n in colnames(aFDR)) {
      nn <- paste0("q", gsub('-', '.', n)) 
      nodeList[[i]][[nn]] <- aFDR[i,n]
    }
    nodeList[[i]]$aFDR <- aFDR[i,]
  }
  attr(buffer, "thresholds") <- attr(aFDR, "thresholds")
  class(buffer) <- c(class(buffer), "mincQvals")
  return(buffer)
}

#' Estimates degrees of freedom
#' 
#' see \link{anatLmerEstimateDF}
#'
#' @param buffer input hierarchical tree
#' @param n number of structures to use for DF estimation
#'
#' @return input tree with df added
#' @export
#' 
#' @details The volumes inside the anatomical hierarchy and the data must be of
#'   the same length and ordering.
#'
#' @examples \dontrun{
#' vols <- addVolumesToHierarchy(hdefs, allvols)
#' hLm <- hanatLmer(~Sex + (1|ID), gf, vols)
#' hLm <- hanatLmerEstimateDF(hLm)
#' hLm <- hanatFDR(hLm)
#' thresholds(hLm)
#' }
hanatLmerEstimateDF <- function(buffer, n=50) {
  modelData <- buffer$Get("model")
  if (is.null(dim(modelData))) {
    dim(modelData) <- c(length(modelData), 1)
  }
  else {
    modelData <- t(modelData)
  }
  modelData <- copyExtraAttributes(buffer, modelData)
  dimnames(modelData) <- attr(buffer, "modelDimnames")
  class(modelData) <- attr(buffer, "modelClass")
  modelData <- anatLmerEstimateDF(modelData, n=n)
  attr(buffer, "df") <- attr(modelData, "df")
  return(buffer)
}

#' Add volumes to an anatomical hierarchy
#' 
#' Takes a tree of anatomical hierarchies and adds volumes from
#' \link{anatGetAll}
#' 
#' @param hdefs The anatomical hierarchy
#' @param volumes The matrix of volumes
#'   
#' @return The input anatomical hierarchy with volumes added
#' @export
#' 
#' @details Each node in the tree has two new attributes: 
#' \itemize{
#' \item{volumes -}{ The matrix of volumes}
#' \item{meanVolume -}{ The mean of the volumes of that node}
#' }
#' 
#' Currently propagates volumes up the tree by summing. Future versions
#'   will add propagating through weighted means or other functions
#'   
#' @examples \dontrun{
#' abijson <- "allen.json"
#' defs <- "Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv"
#' hdefs <- makeMICeDefsHierachical(defs, abijson)
#' allvols <- anatGetAll(gf$filenames)
#' hanat <- addVolumesToHierarchy(hdefs, allvols)
#' }
addVolumesToHierarchy <- function(hdefs, volumes){
  hanat <- Clone(hdefs)
  volLabels <- as.integer(attributes(volumes)$anatIDs)
  hanat$Do(function(x) {
    if (isLeaf(x)) {
      if (x$label %in% volLabels) {
        whichIndex <- which(volLabels == x$label)
        x$volumes <- volumes[, whichIndex]
        x$meanVolume <- mean(x$volumes)
      }
    }
  })
  
  hanat$Do(function(x) x$meanVolume <- Aggregate(x, "meanVolume", function(x) sum(x,na.rm=TRUE)), traversal="post-order")
  hanat$Do(function(x) x$volumes <- Aggregate(x, "volumes", rowSums), 
           traversal="post-order", filterFun = isNotLeaf)
  return(hanat)
}

addHierarchyToDefs <- function(defs, tree) {
  defs$pathString <- NA
  for (i in 1:nrow(defs)) {
    tmp <- FindNode(tree, defs$ABI[i])$Get("path")[[1]]
    tmp <- tmp[tmp != "children2"]
    tmp <- tmp[3:length(tmp)]
    Leaf <- as.character(defs$Leaf[i])
    if (defs$ABI[i] == Leaf | defs$name[i] == Leaf) {
      endpath <- as.character(defs$name[i])
    }
    else {
      endpath <- c(Leaf, as.character(defs$name[i]))
    }
    if (defs$ABI[i] == Leaf & defs$name[i] == defs$ABI[i]) {
      defs$pathString[i] <- paste(tmp, collapse="/")
    }
    else {
      defs$pathString[i] <- paste(c(tmp, endpath), collapse="/")
    }
  }
  return(defs)
}

fixBranching <- function(ldefs, tree) {
  # branching is a problem if a leaf has too many descendants
  for (i in 1:nrow(ldefs)) {
    ldefs$kidCount[i] <- FindNode(tree, ldefs$Leaf[i])$totalCount
    if (ldefs$kidCount[i] > 3) {
      ldefs$Leaf[i] <- paste0(ldefs$Leaf[i], "-other")
      ldefs$name[i] <- paste0(ldefs$name[i], "-other")
    }
  }
  return(ldefs)

}


lateralizeDefs <- function(defs) {
  ldefs <- c()
  
  # some places in the ABI hierarchy are not the leafs of the structures we have. Fix that.
  defs$Leaf <- as.character(defs$ABI)
  dups <- which(duplicated(defs$ABI))
  dupsAll <- which(defs$ABI %in% defs$ABI[dups])
  defs$Leaf[dupsAll] <- as.character(defs$Structure[dupsAll])
  
  # and in some cases a duplicate will also have the same name in the structure column
  # as in the ABI column, so some munging is necessary
  samename <- which(as.character(defs$ABI) == as.character(defs$Structure))
  samenamedups <- samename[samename %in% dupsAll]
  defs$Leaf[samenamedups] <- paste0(as.character(defs$Structure[samenamedups]), "-other")
            
  
  for (i in 1:nrow(defs)) {
    if (defs$left.label[i] == defs$right.label[i]) {
      ldefs <- rbind(ldefs, data.frame(Structure=defs$Structure[i],
                                       side="both",
                                       ABI=defs$ABI[i],
                                       name=defs$Leaf[i],
                                       label=defs$left.label[i],
                                       Leaf=defs$Leaf[i]))
    }
    else {
      ldefs <- rbind(ldefs, data.frame(Structure=defs$Structure[i],
                                       side="left",
                                       ABI=defs$ABI[i],
                                       name=paste("left", defs$Leaf[i]),
                                       label=defs$left.label[i],
                                       Leaf=defs$Leaf[i]))
      ldefs <- rbind(ldefs, data.frame(Structure=defs$Structure[i],
                                       side="right",
                                       ABI=defs$ABI[i],
                                       name=paste("right", defs$Leaf[i]),
                                       label=defs$right.label[i],
                                       Leaf=defs$Leaf[i]))

    }
  }
  
  ldefs$Leaf <- as.character(ldefs$Leaf)
  ldefs$side <- as.character(ldefs$side)
  ldefs$Structure <- as.character(ldefs$Structure)
  ldefs$name <- as.character(ldefs$name)
  ldefs$ABI <- as.character(ldefs$ABI)
  return(ldefs)
}

copyABIinfo <- function(hdefs, abitree) {
  tmpfunc <- function(node) {
    # if node has an ABI name, use it. Otherwise assume that the name of the node
    # is the ABI name
    if (is.null(node$ABI)) {
      node$ABI <- node$name
    }
    abiNode <- FindNode(abitree, node$ABI)
    # if the find failed, then we are probably in the left or right version of a 
    # new name not in the ABI dictionary. So take the first child and use it's 
    # ABI name instead.
    if (is.null(abiNode)) {
      abiNode <- FindNode(abitree, node$children[[1]]$ABI)
    }
    node$ABI <- abiNode$name
    node$color_hex_triplet <- paste0('#', abiNode$color_hex_triplet)
    if (startsWith(node$name, "left")) {
      apref <- "l"
    }
    else if (startsWith(node$name, "right")) {
      apref <- "r"
    }
    else {
      apref <- ""
    }
    node$acronym <- paste0(apref, abiNode$acronym)
  }
  hdefs$Do(fun = tmpfunc)
}


#' Creates a hierarchical anatomical tree
#' 
#' Makes a hierarchical tree using ontogeny from Allen Brain Institute
#' 
#' @param defs CSV file describing anatomical labels. See details about format
#' @param abijson The JSON file from the Allen Brain Institute
#'   
#' @details Takes the label definitions from a CSV file (same as the rest of the
#'   anat family of functions) and places it into a tree (using
#'   \link{data.tree}). The CSV file must thus have an extra column giving the
#'   corresponding name of each structure in the Allen Institute's nomenclature.
#'   
#'   Currently the left and right structures are leaves at the end of the tree,
#'   and so are the first to be combined into the bilateral volumes. It is
#'   possible that in the future keeping hemispheres separate will be an option.
#'   
#' @return a \link{data.tree} object containing the full anatomical hierarchy
#' @export
#' 
#' @examples \dontrun{
#' abijson <- "allen.json"
#' defs <- "Dorr_2008_Steadman_2013_Ullmann_2013_mapping_of_labels.csv"
#' hdefs <- makeMICeDefsHierachical(defs, abijson) 
#' }
makeMICeDefsHierachical <- function(defs, abijson) {
  # read the definitions from the JSON file provided by the Allen institute
  abi <- fromJSON(file=abijson)
  # create a data.tree
  tree <- FromListSimple(abi, check="no-warn")
  # read the CSV defs file
  defs <- read.csv(defs)

  # split the definitions into one per line (i.e. left and right separate)
  ldefs <- lateralizeDefs(defs)
  tree$Do(function(n){ n$name <- gsub("/", " and ", n$name)})
  # create the first version of the anatomical hierarchy
  mh <- addHierarchyToDefs(ldefs, tree)
  treeMH <- as.Node(mh)
  # fix branching odds and ends (i.e. make sure that only leafs are at end of hierarchy)
  ldefs2 <- fixBranching(ldefs, treeMH)
  # recreate the anatomical hierarchy
  mh <- addHierarchyToDefs(ldefs2, tree )
  treeMH <- as.Node(mh)
  # copy important bits from the ABI tree - name and colour especially.
  copyABIinfo(treeMH, tree)
  # propagate the labels so that each structure has a vector of all associated labels
  treeMH$Do(function(x) x$label <- unlist(Aggregate(x, "label", c)), traversal="post-order")

  return(treeMH)
}

colourVector <- function(vals, low, high, symmetric=F,
                         colourScale=colorRampPalette(c("red", "yellow"))(255), 
                         rColourScale=colorRampPalette(c("blue", "turquoise1"))(255), 
                         transparent=NA) {
  svals <- scaleSlice(vals, low, high, underTransparent=T)
  spvals <- scaleSliceToPalette(svals, low, high, colourScale)
  spvals <- colourScale[spvals]
  if (symmetric) {
    svalsN <- scaleSlice(vals, low*-1, high*-1, underTransparent=T)
    spvalsN <- scaleSliceToPalette(svalsN, low*-1, high*-1, rColourScale)
    spvalsN <- rColourScale[spvalsN]
    spvals <- paste0(spvals, spvalsN)
    spvals <- gsub('NA', '', spvals)
    spvals[spvals == ""] <- NA
  }
  spvals[is.na(spvals)] <- transparent
  return(spvals)
}

#' Turns an anatomical tree into a nodes and edges for visNetwork plotting
#'
#' @param hanatTree The input anatomical tree
#' @param colourVariable String containing the variable to colour nodes with
#' @param colourScale Colour scale to use
#' @param rColourScale Reverse colour scale to use if symmetric colours enabled
#' @param low Lower cut-off for colour scale
#' @param high Upper cut-off for colour scale
#' @param symmetric Boolean for whether colour scale is symmetric
#' @param transparent Colour to use for transparent nodes
#' @param edgeColourFromABI Whether to set edge colours based on ABI atlas
#' @param fontsize The font size for \code{hanatView}
#' @param levelSeparation Spacing between hierarchical layers for \code{hanatView}
#' @param ... Extra arguments to \code{hanatToVisGraph} from \code{hanatView}
#' 
#' @return list with nodes and edges data frames
#' @export
#'
## @examples
hanatToVisGraph <- function(hanatTree, 
                            colourVariable="color_hex_triplet", 
                            colourScale=colorRampPalette(c("red", "yellow"))(255), 
                            rColourScale=colorRampPalette(c("blue", "turquoise1"))(255),
                            low=NULL,
                            high=NULL,
                            symmetric=F,
                            transparent="#FDFDFD",
                            edgeColourFromABI=F) {
  tree <- Clone(hanatTree)
  # check if colourVariable is a number; if so, colourize it
  if (is.double(tree[[colourVariable]])) {
    colourVals <- tree$Get(colourVariable)
    tree$Set(plotcolour=colourVector(colourVals, low, high, 
                                     symmetric = symmetric, 
                                     colourScale=colourScale,
                                     rColourScale=rColourScale,
                                     transparent=transparent))
    tree$Set(stat = round(colourVals, 3))
    plotcolour <- "plotcolour"
  }
  # if not a number, assume that it's a colour itself (likely the ABI colours)
  else {
    plotcolour <- colourVariable
  }
  nodes <- as.data.frame(tree, 
                         id="name", 
                         color=plotcolour, 
                         label="name", 
                         stat="stat")
  
  # set the hover over title
  if (is.double(tree[[colourVariable]])) {
    nodes$title <- paste0("<p><b>", nodes$label, "</b><br>", round(nodes$stat, 3), "</p>")
  }
  else {
    nodes$title <- paste0("<p><b>", nodes$label, "</b></p>")
  }
  
  # create the edges
  edges <- data.frame()
  nn <- Traverse(tree, pruneFun = isNotLeaf)
  for (i in 1:length(nn)) {
    for (j in 1:length(nn[[i]]$children)) {
      edf <- data.frame(from=nn[[i]]$name, to=nn[[i]]$children[[j]]$name)
      if (edgeColourFromABI) {
        edf$color <- nn[[i]]$children[[j]]$color_hex_triplet
      }
      edges <- rbind(edges, edf)
    }
  }
  return(list(nodes=nodes, edges=edges))
}

#' @describeIn hanatToVisGraph Directly view the graph
#' @export
hanatView <- function(..., fontsize=14, levelSeparation=500) {
  gv <- hanatToVisGraph(...)
  visNetwork(gv$nodes, gv$edges) %>% 
    visNodes(shape="box", font=list(size=fontsize)) %>% 
    visEdges(width=10, smooth=F, arrows="to") %>% 
    visHierarchicalLayout(direction="LR", levelSeparation = levelSeparation, sortMethod="directed") %>% 
    visPhysics(enabled=F)
}

#' hanatPlot
#'
#' Create a basic plot of the anatomy tree
#' 
#' @param anatTree the input anatomical tree
#' @export
hanatPlot <- function(anatTree) {
  tree <- Clone(anatTree)
  tree$Do(function(x) x$name <- gsub("'", "", x$name))
  SetGraphStyle(tree, layout="dot", rankdir="LR")

  plot(tree)
}

