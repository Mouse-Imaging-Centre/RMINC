\name{rminc.checkForExternalProgram}
\Rdversion{1.1}
\alias{rminc.checkForExternalProgram}

\title{
    Check the Availability of an External Program or Script
}
\description{
    Test to see whether a required program or script is installed and available on the current system. Program or script status is reflected by a returned boolean.
}

\usage{
rminc.checkForExternalProgram(program, test_string, prog_options="")
}

\arguments{
	\item{program}{A string specifying the program or script name.  If the program/script is not expected to be on the users PATH, then the program/script name must contain a fully-qualified path.}
	\item{test_string}{A sub-string contained within the output of the given program/script, when executed with the provided option(s), if all is functioning properly.}
	\item{prog_options}{A string, usually containing a simple program/script option, which would simply test whether the program/script was found and can execute. The most obvious options would include: ``-version'', and ``-help''.}
}

\details{
This function is passed the name of a program or script that is s'posed 
to be on the user's path, along with an option that generates a known
response (the test_string).  If the passed test_string is not found in
the returned output, we send a warning message and then return FALSE.
The user, given a FALSE, can then cobble together a fitting response to
the user.
}
\value{A boolean: TRUE if the test_string was found in the test execution, else FALSE.}

\author{Jim Nikelski \email{nikelski@bic.mni.mcgill.ca}}


\examples{
\dontrun{
library(RMINC)

program <- "xfm2param_nonexisting"
progOptions <- "-version"
test_string <- "mni_autoreg"
status <- rminc.checkForExternalProgram(program, test_string, progOptions)
if ( !status ) { ... }

}}







