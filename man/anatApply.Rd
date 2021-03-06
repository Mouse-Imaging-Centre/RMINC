% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minc_anatomy.R
\name{anatApply}
\alias{anatApply}
\title{Apply function over anat structure}
\usage{
anatApply(vols, grouping = NULL, method = mean, ...)
}
\arguments{
\item{vols}{anatomy volumes}

\item{grouping}{a factor grouping with which to perform operations}

\item{method}{The function which to apply [default mean]}

\item{...}{Extra arguments to the function passed to method}
}
\value{
out: The output will be a single vector containing as many
         elements as there are regions in the input variable by the number of groupings
}
\description{
This function is used to compute an arbitrary function of every region in an anat structure.
}
\examples{
\dontrun{
getRMINCTestData() 
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)
vm <- anatApply(gf$lobeThickness,gf$Primary.Diagnosis)
}
}
