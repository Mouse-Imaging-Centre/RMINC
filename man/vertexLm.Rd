\name{vertexLm}
\alias{vertexLm}
\title{Calculates statistics and coefficients for linear model of specified vertex files}
\usage{
  vertexLm(formula, data, subset = NULL)
}
\arguments{
  \item{formula}{a model formula}

  \item{data}{a data.frame containing variables in formula}

  \item{subset}{rows to be used, by default all are used}
}
\value{
  Returns an object containing the F and t statistcs that
  can be passed directly into vertexFDR. The coefficients
  can be found in the attribute 'coefficients'
}
\description{
  Calculates statistics and coefficients for linear model
  of specified vertex files
}
\examples{
gf = read.csv("~/SubjectTable.csv")
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
gf$vertexFiles = as.factor(gf$CIVETFILES$nativeRMStlink20mmleft)
result = vertexLm(vertexFiles~Primary.Diagnosis,gf)
vertexFDR(result)
coefficients <- attr(results,'coefficients')
}
\seealso{
  mincLm,anatLm,vertexFDR
}

