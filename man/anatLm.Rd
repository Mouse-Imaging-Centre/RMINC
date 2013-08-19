\name{anatLm}
\alias{anatLm}
\title{Calculates statistics and coefficients for linear model of specified anat structure}
\usage{
  anatLm(formula, data, anat, subset = NULL)
}
\arguments{
  \item{formula}{a model formula}

  \item{data}{a data.frame containing variables in formula}

  \item{anat}{an array of atlas labels vs subject data}

  \item{subset}{rows to be used, by default all are used}
}
\value{
  Returns an object containing the F and t statistcs that
  can be passed directly into anatFDR. The coefficients can
  be found in the attribute 'coefficients'
}
\description{
  Calculates statistics and coefficients for linear model
  of specified anat structure
}
\examples{
gf = read.csv("~/SubjectTable.csv")
civet.getAllFilenames(gf,"ID","ABC123","~/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("~/Atlases/AAL/AAL.csv",gf)
result = anatLm(~Primary.Diagnosis,gf,gf$lobeVolume)
anatFDR(result)
coefficients <- attr(results,'coefficients')
}
\seealso{
  mincLm,anatLm,anatFDR
}

