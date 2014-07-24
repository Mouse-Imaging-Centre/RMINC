context("writeVertex")

gf = read.csv("/home/dcassel/Projects/POND/MR160/SubjectInfo/POND-imaging.csv")
gf = civet.getAllFilenames(gf,"POND.ID","POND","~/Projects/POND/MR160/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/home/dcassel/resource/Atlases/AAL/AAL.csv",gf)

