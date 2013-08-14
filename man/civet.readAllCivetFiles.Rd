\name{civet.readAllCivetFiles}
\alias{civet.readAllCivetFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read  all CIVET files into R
}
\description{
Parses outputs from CIVET pipeline 
}
\usage{
civet.readAllCivetFiles(atlasFile, gf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{atlasFile}{
Full path to the atlas on which CIVET was run ;
}
  \item{gf}{
Data Frame containing list of all CIVET file names, and where results will be stored
}
}
\details{
Prior to running, civet.getAllFilenames may be called to generate the input argument gf . This function will extract the following information from the CIVET pipeline: Lobe Area (40 mm), Lobe Thickness, Lobe Volume, GI, Mid Surface Native Area, Surface Native Volume, Native RMS RSL tLink (20mm), Native RMS tLink (20 mm)
}
\value{
gf is returned with CIVET values
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Daniel Cassel (daniel.cassel@sickkids.ca)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 civet.getAllFilenames
}
\examples{
gf = read.csv("~/SubjectTable.csv") 
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
}
