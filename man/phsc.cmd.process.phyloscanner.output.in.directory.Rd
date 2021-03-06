% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloscan.fun.private.R
\name{phsc.cmd.process.phyloscanner.output.in.directory}
\alias{phsc.cmd.process.phyloscanner.output.in.directory}
\title{Generate bash commands to process phyloscanner output}
\usage{
phsc.cmd.process.phyloscanner.output.in.directory(
  tmp.dir,
  file.patients,
  pty.args
)
}
\arguments{
\item{tmp.dir}{Directory with phyloscanner output.}

\item{file.patients}{File name of the file that contains the list of unique individuals/units that the bam files correspond to, possibly after re-naming.}

\item{pty.args}{List of phyloscanner input variables.}
}
\value{
character string of bash commands.
}
\description{
This function generates bash commands that combine the various Rscripts in the phyloscanner toolkit
}
\keyword{internal}
