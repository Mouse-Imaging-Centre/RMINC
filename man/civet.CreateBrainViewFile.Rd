\name{civet.CreateBrainViewFile}
\alias{civet.CreateBrainViewFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
civet.CreateBrainViewFile
}
\description{
Creates a text file that can be loaded into brain-view2
}
\usage{
civet.CreateBrainViewFile(dataFile,AtlasFile,AtlasVertexFile,outputFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataFile}{
Either the name of a file with atlas labeling or an R array with atlas labeling
}
  \item{atlasFile}{
Text file with map between atlas labels and numbers
}
  \item{atlasVertices}{
Text file with map between vertex points and atlas numbers
}
  \item{outputFileName}{
path to file where output will be saved 
}
  \item{civetVersion}{
Version of CIVET used (Default 1.1.12)
}
}
\details{
This function will create a txt file that can be loaded into brain-view2, in order to visualize the results from CIVET. This function either accepts a text file or R variable as input. If using an R variable the rows/columns must be labeled with the Atlas Names. The names are then matched to numbers as given by the AtlasFile, and then the numbers are matched to vertices in the AtlasVertexFile. 
}
\value{
A file is generated (specified by outputFile)
}
\references{
}
\author{
Daniel Cassel (daniel.cassel@sickkids.ca)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
gf = read.csv("~/SubjectTable.csv") 
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
civet.CreateBrainViewFile(gf$lobeThickness[1,],"/Atlases/AAL/AAL.csv","/Atlases/AAL/AAL_atlas_left.txt","leftLobeThickness.txt")
}

