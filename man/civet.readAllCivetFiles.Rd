\name{civet.readAllCivetFiles}
\alias{civet.readAllCivetFiles}
\title{Read  all CIVET files into R}
\usage{
  civet.readAllCivetFiles(atlasFile, gf)
}
\arguments{
  \item{gf}{Data Frame containing list of all CIVET file
  names, and where results will be stored}

  \item{atlasFile}{Full path to the atlas on which CIVET
  was run ;}
}
\value{
  gf is returned with CIVET values
}
\description{
  Parses outputs from CIVET pipeline
}
\details{
  Prior to running, civet.getAllFilenames may be called to
  generate the input argument gf . This function will
  extract the following information from the CIVET
  pipeline: Lobe Area (40 mm), Lobe Thickness, Lobe Volume,
  GI, Mid Surface Native Area, Surface Native Volume,
  Native RMS RSL tLink (20mm), Native RMS tLink (20 mm)
}
\examples{
gf = read.csv("~/SubjectTable.csv")
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
}
\seealso{
  civet.getAllFilenames
}

