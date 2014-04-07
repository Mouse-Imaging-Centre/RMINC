context("civet.readAllCivetFiles")

gf = read.csv("/home/dcassel/Projects/POND/MR160/SubjectInfo/POND-imaging.csv")
gf = civet.getAllFilenames(gf,"POND.ID","POND","/home/dcassel/Projects/POND/MR160/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/home/dcassel/Atlases/AAL/AAL.csv",gf)
dataFile = gf$lobeThickness
AALAtlas = read.csv("/home/dcassel/Atlases/AAL/AAL.csv")
verticesL = read.csv("/home/dcassel/Atlases/AAL/AAL_atlas_left.txt")
verticesR = read.csv("/home/dcassel/Atlases/AAL/AAL_atlas_right.txt")
reducedVertices = c(verticesL[0:40961,1],verticesR[0:40961,1])

atlasIndex = pmatch(names(dataFile[1,1]),AALAtlas[,3])
reducedVerticesIndices = which(reducedVertices == AALAtlas[atlasIndex,1],arr.ind=FALSE)
meanThicknessFromVertexFile = mean(gf$nativeRMS_RSLtlink[1,reducedVerticesIndices])
meanThicknessFromVertexFile
test_that("Mean Thickness from Vertex File is the same as thickness from Anat File",{
	expect_that(meanThicknessFromVertexFile,is_equivalent_to(gf$lobeThickness[1,1]))
})

