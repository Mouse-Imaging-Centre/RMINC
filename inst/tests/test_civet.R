library("testthat")
context("civet.readAllCivetFiles")

getRMINCTestData()
dataPath <- file.path(tempdir(), "rminctestdata/")

gf <- read.csv(file.path(dataPath, "CIVET_TEST.csv"))
gf <- civet.getAllFilenames(gf,"ID","TEST",file.path(dataPath, "CIVET"),"TRUE","1.1.12")
gf <- civet.readAllCivetFiles(file.path(dataPath, "AAL.csv"),gf)

dataFile <- gf$lobeThickness
AALAtlas <- read.csv(file.path(dataPath, "AAL.csv"))
verticesL <- read.csv(file.path(dataPath, "AAL_atlas_left.txt"),header = FALSE)
verticesR <- read.csv(file.path(dataPath, "AAL_atlas_right.txt"),header = FALSE)
reducedVertices <- c(verticesL[0:40961,1],verticesR[0:40961,1])

atlasIndex <- pmatch(names(dataFile[1,1]),AALAtlas[,3])
reducedVerticesIndices <- which(reducedVertices == AALAtlas[atlasIndex,1],arr.ind=FALSE)
meanThicknessFromVertexFile <- mean(gf$nativeRMS_RSLtlink[1,reducedVerticesIndices])

test_that("Mean Thickness from Vertex File is the same as thickness from Anat File",{
	expect_that(meanThicknessFromVertexFile,equals(gf$lobeThickness[[1,1]],tolerance = 0.01))
})

