\name{civet.getAllFilenames}
\alias{civet.getAllFilenames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
civet.getAllFilenames
}
\description{
Generates list of filenames output by CIVET
}
\usage{
civet.getAllFilenames(gf, idvar, prefix, basedir, append = TRUE, civetVersion = "1.1.9")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gf}{
Data Frame with subject information
}
  \item{idvar}{
column name in gf with subject IDs
}
  \item{prefix}{
Prefix specified when CIVET was run
}
  \item{basedir}{
directory where all CIVET output was stored
}
  \item{append}{
Whether to append the results to the input gf
}
  \item{civetVersion}{
Version of CIVET 
}
}
\details{
Prior to running, read.csv  may be called to generate the input argument gf. The results will be stored under the column name CIVETFILES either in the input gf (if append = TRUE) or in a new gf. Currently only CIVET versions 1.1.9 and 1.1.12 are supported.

}
\value{
gf is returned with CIVET filenames 
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
civet.readAllCivetFiles
}
\examples{
gf = read.csv("~/SubjectTable.csv") 
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12") 
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
}

