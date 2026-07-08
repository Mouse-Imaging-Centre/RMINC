# Package index

## Package overview

- [`RMINC`](https://mouse-imaging-centre.github.io/RMINC/reference/RMINC-package.md)
  [`RMINC-package`](https://mouse-imaging-centre.github.io/RMINC/reference/RMINC-package.md)
  : RMINC: R interface to the MINC universe

## Voxel-wise statistics

Massively-univariate models fit independently at every voxel.

- [`mincLm()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLm.md)
  : Linear model at Every Voxel
- [`mincAnova()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincAnova.md)
  : Voxel-wise ANOVA
- [`mincTtest()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincTtest.md)
  : Minc T-test
- [`mincPairedTtest()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincPairedTtest.md)
  : Minc Paired T Test
- [`mincWilcoxon()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincWilcoxon.md)
  : Minc Wilcoxon
- [`mincSummary()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincMean()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincVar()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincSum()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincSd()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincCorrelation()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  : Minc Voxel Summary Functions
- [`mincLogLikRatio()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLogLikRatio.md)
  : run log likelihood ratio tests for different mincLmer objects
- [`mincLogLikRatioParametricBootstrap()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLogLikRatioParametricBootstrap.md)
  : computes parametric bootstrap on mincLogLikRatio output

## Voxel summaries & effect sizes

- [`mincSummary()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincMean()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincVar()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincSum()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincSd()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  [`mincCorrelation()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSummary.md)
  : Minc Voxel Summary Functions
- [`mincTable()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincTable.md)
  : Read a collection of minc volumes
- [`vertexEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  [`mincEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  [`anatEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  : Effect Sizes

## Voxel mixed-effects models

- [`mincLmer()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLmer.md)
  : mincified version of lmer from lme4
- [`mincLmerEstimateDF()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLmerEstimateDF.md)
  : Estimate the degrees of freedom for parameters in a mincLmer model

## Multiple comparisons (voxel)

- [`mincFDR()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFDR.md)
  : False Discovery Rates

- [`mincFDRMask()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFDRMask.md)
  : mincFDRMask

- [`mincFDRThresholdVector()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFDRThresholdVector.md)
  : a utility function to compute thresholds

- [`mincTFCE()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincTFCE.md)
  : Threshold Free Cluster Enhancement

- [`mincRandomize()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincRandomize.md)
  :

  Run a permutation test on a `mincLm` result

- [`thresholds()`](https://mouse-imaging-centre.github.io/RMINC/reference/thresholds.md)
  : Get Probability Thresholds

## Surface / vertex statistics

Vertex-wise models for cortical surface data.

- [`vertexLm()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexLm.md)
  : Calculates statistics and coefficients for linear model of specified
  vertex files
- [`vertexLmer()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexLmer.md)
  : Vertex Mixed Effects Models
- [`vertexLmerEstimateDF()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexLmerEstimateDF.md)
  : Estimate the degrees of freedom for parameters in a vertexLmer model
- [`vertexAnova()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexAnova.md)
  : Performs ANOVA on each vertex point specified
- [`vertexApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexApply.md)
  : Apply function over vertex Files
- [`vertexAtlasApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexAtlasApply.md)
  : Apply a structure summary function across vertices
- [`vertexEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  [`mincEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  [`anatEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  : Effect Sizes
- [`vertexFDR()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexFDR.md)
  : Vertex False Discovery Rates
- [`vertexFindPeaks()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexFindPeaks.md)
  : Vertex find peaks
- [`vertexLookup()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexLookup.md)
  : Vertex Lookup
- [`vertexSelect()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexSelect.md)
  : Select Vertices From a Surface
- [`vertexMean()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexSummaries.md)
  [`vertexSum()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexSummaries.md)
  [`vertexVar()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexSummaries.md)
  [`vertexSd()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexSummaries.md)
  : Create descriptive statistics across a series of vertex files
- [`vertexTFCE()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexTFCE.md)
  : Threshold Free Cluster Enhancement
- [`vertexTable()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexTable.md)
  : Create a table of vertex values
- [`writeVertex()`](https://mouse-imaging-centre.github.io/RMINC/reference/writeVertex.md)
  : Writes vertex data to a file with an optional header

## ROI / anatomy analysis (flat)

- [`anatLm()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatLm.md)
  : Calculates statistics and coefficients for linear model of specified
  anat structure
- [`anatLmer()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatLmer.md)
  : Anatomy Linear Mixed Effects Model
- [`anatLmerEstimateDF()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatLmerEstimateDF.md)
  : Estimate Degrees of freedom for an anatLmer model
- [`anatAnova()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatAnova.md)
  : Performs ANOVA on each region specified
- [`anatApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatApply.md)
  : Apply function over anat structure
- [`anatCombineStructures()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatCombineStructures.md)
  : Combine left and right volumes
- [`anatCreateVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatCreateVolume.md)
  : Create a Statistic Volume
- [`vertexEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  [`mincEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  [`anatEffectSize()`](https://mouse-imaging-centre.github.io/RMINC/reference/vertexEffectSize.md)
  : Effect Sizes
- [`anatFDR()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatFDR.md)
  : Anatomy False Discovery Rates
- [`anatGetAll()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll.md)
  : Faster AnatGet
- [`anatGetAll_old()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll_old.md)
  : Get values given a set of files and an atlas
- [`as.data.frame(`*`<anatMatrix>`*`)`](https://mouse-imaging-centre.github.io/RMINC/reference/as.data.frame.anatMatrix.md)
  : Create Anatomy data.frame
- [`voxel_atlas_defs`](https://mouse-imaging-centre.github.io/RMINC/reference/voxel_atlas_defs.md)
  : Voxel Atlas Definitions
- [`anatMean()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatSummaries.md)
  [`anatSum()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatSummaries.md)
  [`anatVar()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatSummaries.md)
  [`anatSd()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatSummaries.md)
  : anatSummaries
- [`anatSummarize()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatSummarize.md)
  : Summarize anatGetAll results by hierarchy
- [`anatToVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/anatToVolume.md)
  : Converts a column from an anatLm model into a volume for viewing or
  saving

## Hierarchical anatomy atlases

- [`hanatAnova()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatAnova.md)
  : ANOVA applied to every item in an anatomical tree
- [`hanatFDR()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatFDR.md)
  : Compute the False Discovery Rate for an anatomical hierarchy
- [`hanatLm()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatLm.md)
  : linear model applied to every item in an anatomical tree
- [`hanatLmer()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatLmer.md)
  : linear mixed effects model applied to every item in an anatomical
  tree
- [`hanatLmerEstimateDF()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatLmerEstimateDF.md)
  : Estimates degrees of freedom
- [`hanatPlot()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatPlot.md)
  : hanatPlot
- [`hanatToVisGraph()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatToVisGraph.md)
  [`hanatView()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatToVisGraph.md)
  : Turns an anatomical tree into a nodes and edges for visNetwork
  plotting
- [`hanatToVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/hanatToVolume.md)
  : Converts a column from a hierarchical tree into a volume for viewing
  or saving
- [`addVolumesToHierarchy()`](https://mouse-imaging-centre.github.io/RMINC/reference/addVolumesToHierarchy.md)
  : Add volumes to an anatomical hierarchy
- [`makeMICeDefsHierachical()`](https://mouse-imaging-centre.github.io/RMINC/reference/makeMICeDefsHierachical.md)
  : Creates a hierarchical anatomical tree

## Model selection

- [`compare_models()`](https://mouse-imaging-centre.github.io/RMINC/reference/compare_models.md)
  : Compare a set of massively-univariate models
- [`AICc()`](https://mouse-imaging-centre.github.io/RMINC/reference/AICc.md)
  : Compute the Corrected AIC for an object
- [`flexLm()`](https://mouse-imaging-centre.github.io/RMINC/reference/flexLm.md)
  : Flexible Linear Model with Voxel-Varying Covariates

## Parallel & HPC execution

batchtools-based local and cluster parallelism.

- [`mincApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincApply.md)
  : Apply arbitrary R function at every voxel
- [`mincApplyRCPP()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincApplyRCPP.md)
  : Perform Arbitrary calculations on a collection of mincVolumes
- [`pMincApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/pMincApply.md)
  : Parallel MincApply
- [`mcMincApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/mcMincApply.md)
  : Local multicore mincApply
- [`qMincApply()`](https://mouse-imaging-centre.github.io/RMINC/reference/qMincApply.md)
  [`qMincRegistry()`](https://mouse-imaging-centre.github.io/RMINC/reference/qMincApply.md)
  [`qMincMap()`](https://mouse-imaging-centre.github.io/RMINC/reference/qMincApply.md)
  [`qMincReduce()`](https://mouse-imaging-centre.github.io/RMINC/reference/qMincApply.md)
  : True cluster mincApply
- [`runRMINCTestbed()`](https://mouse-imaging-centre.github.io/RMINC/reference/runRMINCTestbed.md)
  : Run Testbed
- [`verboseRun()`](https://mouse-imaging-centre.github.io/RMINC/reference/verboseRun.md)
  : Run function with/without output silenced

## I/O & low-level access

- [`mincGetVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincGetVolume.md)
  : Read a MINC file
- [`mincWriteVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincWriteVolume.md)
  : Write a MINC volume to file
- [`mincArray()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincArray.md)
  : A utility function to give a MINC object spatial dimensions
- [`mincGetVoxel()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincGetVoxel.md)
  : Retrieve Voxel Values
- [`mincGetWorldVoxel()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincGetWorldVoxel.md)
  : World Vector
- [`mincVectorToVoxelCoordinates()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincVectorToVoxelCoordinates.md)
  : converts a vector index to the voxel indices in MINC
- [`mincConvertTagToMincArrayCoordinates()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincConvertTagToMincArrayCoordinates.md)
  : convert tags from world coordinates to mincArray coordinates
- [`mincConvertVoxelToWorld()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincConvertVoxelToWorld.md)
  : Voxel to World
- [`mincConvertWorldToVoxel()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincConvertWorldToVoxel.md)
  : World to Voxel
- [`mincConvertVoxelMatrix()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincConvertVoxelMatrix.md)
  : Convert Voxel to World Coordinates
- [`mincConvertWorldMatrix()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincConvertWorldMatrix.md)
  : Convert World to Voxel Coordinates
- [`minc.dimensions.sizes()`](https://mouse-imaging-centre.github.io/RMINC/reference/minc.dimensions.sizes.md)
  : Dimension Sizes
- [`minc.separation.sizes()`](https://mouse-imaging-centre.github.io/RMINC/reference/minc.separation.sizes.md)
  : Slice Separations
- [`minc.get.hyperslab()`](https://mouse-imaging-centre.github.io/RMINC/reference/minc.get.hyperslab.md)
  : Get a hyperslab from a MINC2 file
- [`minc.get.history()`](https://mouse-imaging-centre.github.io/RMINC/reference/minc_history.md)
  [`minc.append.history()`](https://mouse-imaging-centre.github.io/RMINC/reference/minc_history.md)
  : Minc History
- [`mincAttributes()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincAttributes.md)
  [`setMincAttributes()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincAttributes.md)
  : Get and Set Minc Specific Attributes
- [`mincGetTagFile()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincGetTagFile.md)
  : read a tag file
- [`mincSelectRandomVoxels()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincSelectRandomVoxels.md)
  : selects a few random indices from a volume
- [`mincFindPeaks()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFindPeaks.md)
  : Finds peak values in a MINC file or volume
- [`mincLabelPeaks()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincLabelPeaks.md)
  : label peaks with the name of the atlas structure they are in
- [`mincFileCheck()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincFileCheck.md)
  : Minc File Check
- [`mincGetMask()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincGetMask.md)
  : Minc Masks
- [`mincGetVector()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincGetVector.md)
  : Voxel Vector
- [`extract_column()`](https://mouse-imaging-centre.github.io/RMINC/reference/extract_column.md)
  : Internal function to read in a table or csv file potentially with
  many columns and extract just the one we need
- [`volume.combineMaskVolumes()`](https://mouse-imaging-centre.github.io/RMINC/reference/volume.combineMaskVolumes.md)
  : Combine Multiple Mask Volumes into a Single Mask Volume
- [`volume.explodeLabelVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/volume.explodeLabelVolume.md)
  : Explode a Label Volume into its Components

## MINC object helpers

- [`as.minc()`](https://mouse-imaging-centre.github.io/RMINC/reference/as.minc.md)
  : Coerce to RMINC object
- [`is.minc()`](https://mouse-imaging-centre.github.io/RMINC/reference/is.minc.md)
  : Check for RMINC objects
- [`likeVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/likeVolume.md)
  [`` `likeVolume<-`() ``](https://mouse-imaging-centre.github.io/RMINC/reference/likeVolume.md)
  : Get or set a likeVolume
- [`maskFile()`](https://mouse-imaging-centre.github.io/RMINC/reference/maskFile.md)
  [`` `maskFile<-`() ``](https://mouse-imaging-centre.github.io/RMINC/reference/maskFile.md)
  : Get or set a mask file
- [`simplify2minc()`](https://mouse-imaging-centre.github.io/RMINC/reference/simplify2minc.md)
  : Collate Minc
- [`simplify_masked()`](https://mouse-imaging-centre.github.io/RMINC/reference/simplify_masked.md)
  : Simplify Masked Results
- [`setRMINCMaskedValue()`](https://mouse-imaging-centre.github.io/RMINC/reference/setRMINCMaskedValue.md)
  : Set Masked Value
- [`RMINC_MASKED_VALUE`](https://mouse-imaging-centre.github.io/RMINC/reference/RMINC_MASKED_VALUE.md)
  : RMINC Masked Value

## Visualization

- [`mincPlotSliceSeries()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincPlotSliceSeries.md)
  : MINC Slice Series
- [`mincTriplanarSlicePlot()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincTriplanarSlicePlot.md)
  : Plot Slice Along Each Axis
- [`mincPlotAnatAndStatsSlice()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincPlotAnatAndStatsSlice.md)
  : Anatomy and Statistics Slice
- [`mincPlotPeak()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincPlotPeak.md)
  : Plotting of peaks
- [`mincImage()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincImage.md)
  : Plot a slice from a MINC volume
- [`mincContour()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincContour.md)
  : Draw contour lines from a MINC volume
- [`mincRayTraceStats()`](https://mouse-imaging-centre.github.io/RMINC/reference/mincRayTraceStats.md)
  : Create an image of a statistical peak.
- [`minc.ray.trace()`](https://mouse-imaging-centre.github.io/RMINC/reference/minc.ray.trace.md)
  : Call ray_trace to get an image of a rendered slice
- [`add_colour_bar()`](https://mouse-imaging-centre.github.io/RMINC/reference/add_colour_bar.md)
  : Add a colour bar for a mesh
- [`add_opacity()`](https://mouse-imaging-centre.github.io/RMINC/reference/add_opacity.md)
  : Add opacity to a mesh
- [`colour_mesh()`](https://mouse-imaging-centre.github.io/RMINC/reference/colour_mesh.md)
  : Colourize a mesh
- [`map_to_colours()`](https://mouse-imaging-centre.github.io/RMINC/reference/map_to_colours.md)
  : Generate a vector of colours from a map
- [`lut_to_palette()`](https://mouse-imaging-centre.github.io/RMINC/reference/lut_to_palette.md)
  : A tool that returns a color function/palette from color lookup files
- [`obj_montage()`](https://mouse-imaging-centre.github.io/RMINC/reference/obj_montage.md)
  : Brain Montage
- [`obj_to_graph()`](https://mouse-imaging-centre.github.io/RMINC/reference/obj_to_graph.md)
  : Object to Graph
- [`plot(`*`<bic_lines>`*`)`](https://mouse-imaging-centre.github.io/RMINC/reference/plot.bic_lines.md)
  : Plot A bic_lines object
- [`plot(`*`<bic_obj>`*`)`](https://mouse-imaging-centre.github.io/RMINC/reference/plot.bic_obj.md)
  : Plot a BIC obj
- [`plot(`*`<obj_mesh>`*`)`](https://mouse-imaging-centre.github.io/RMINC/reference/plot.obj_mesh.md)
  : Plot an BIC obj mesh
- [`launch_shinyRMINC()`](https://mouse-imaging-centre.github.io/RMINC/reference/launch_shinyRMINC.md)
  : launch a shiny based inspector

## Meshes & manifold smoothing

- [`create_mesh()`](https://mouse-imaging-centre.github.io/RMINC/reference/create_mesh.md)
  : Create surface mesh for a bic_obj
- [`read_obj()`](https://mouse-imaging-centre.github.io/RMINC/reference/read_obj.md)
  : Read BIC .obj files
- [`read_line_obj()`](https://mouse-imaging-centre.github.io/RMINC/reference/read_line_obj.md)
  : Read a BIC-obj line file
- [`closestVertex()`](https://mouse-imaging-centre.github.io/RMINC/reference/closestVertex.md)
  : Find Closest Vertex
- [`connected_components()`](https://mouse-imaging-centre.github.io/RMINC/reference/connected_components.md)
  : Find connected components of a object graph
- [`components_to_map()`](https://mouse-imaging-centre.github.io/RMINC/reference/components_to_map.md)
  : Components to map
- [`clean_up_manifold_triangles()`](https://mouse-imaging-centre.github.io/RMINC/reference/clean_up_manifold_triangles.md)
  : Clean up manifold mesh
- [`line_obj_to_voxel()`](https://mouse-imaging-centre.github.io/RMINC/reference/line_obj_to_voxel.md)
  : Convert Lines to Voxel Coordinates
- [`laplace_beltrami_operator()`](https://mouse-imaging-centre.github.io/RMINC/reference/laplace_beltrami_operator.md)
  : Compute Laplace-Beltrami operator
- [`laplace_beltrami_smoothing()`](https://mouse-imaging-centre.github.io/RMINC/reference/laplace_beltrami_smoothing.md)
  : Smoothing on manifold
- [`check_step_sizes()`](https://mouse-imaging-centre.github.io/RMINC/reference/check_step_sizes.md)
  : Check file step sizes

## CIVET

- [`civet.CreateBrainViewFile()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.CreateBrainViewFile.md)
  : Create a brain view file
- [`civet.CreateBrainViewROI()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.CreateBrainViewROI.md)
  : civet.CreateBrainViewROI
- [`civet.checkVersion()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.checkVersion.md)
  : Check for a known Civet Version
- [`civet.computeNativeToStxRescalingFactor()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.computeNativeToStxRescalingFactor.md)
  : Compute Native to Stereotactic Rescaling Factor
- [`civet.computeStxTissueVolumes()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.computeTissueVolumes.md)
  [`civet.computeNativeTissueVolumes()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.computeTissueVolumes.md)
  : Compute GM, WM, and CSF Tissue Volumes
- [`civet.flattenForDplyr()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.flattenForDplyr.md)
  : Flatten CIVET results for dplyr
- [`civet.getAllFilenames()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getAllFilenames.md)
  : civet.getAllFilenames
- [`civet.getFilenameClassify()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameGrayMatterPve()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameWhiteMatterPve()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameCsfPve()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameStxT1()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameCerebrumMask()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameSkullMask()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameGrayMatterSurfaces()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameWhiteMatterSurfaces()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameMidSurfaces()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenamesCorticalThickness()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenamesCorticalArea()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenamesCorticalVolume()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameMeanSurfaceCurvature()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameLinearTransform()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  [`civet.getFilenameNonlinearTransform()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.getFilename.md)
  : Get Selected Civet Filenames
- [`civet.organizeCivetDatFilesAtlas()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.organizeCivetDatFilesAtlas.md)
  : Organizes CIVET .dat files based on an Atlas
- [`civet.readAllCivetFiles()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.readAllCivetFiles.md)
  : Read all CIVET files into R
- [`civet.readCBRAIN()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.readCBRAIN.md)
  : Read a CBRAIN project in R
- [`civet.readCivetDatFiles()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.readCivetDatFiles.md)
  : Read Civet-Generated Dat Files
- [`civet.readQC()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.readQC.md)
  : Read CIVET QC data
- [`civet.vertexFilenames()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.vertexFilenames.md)
  : Assemble vertex files for a CIVET run
- [`civet.vertexTable()`](https://mouse-imaging-centre.github.io/RMINC/reference/civet.vertexTable.md)
  : Create a table of vertex measures

## Test data & helpers

- [`getRMINCTestData()`](https://mouse-imaging-centre.github.io/RMINC/reference/getRMINCTestData.md)
  : Download Example Data
- [`pt2()`](https://mouse-imaging-centre.github.io/RMINC/reference/pt2.md)
  : Two tailed version of pt
