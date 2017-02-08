\name{diffLogoTable}
\alias{diffLogoTable}
\title{Draw DiffLogo-table}
\usage{
  diffLogoTable(PWMs, alphabet = DNA,
    configuration = list(), ...)
}
\arguments{
  \item{PWMs}{a list/vector of position weight matrices
  (PWMs) each of type pwm, data.frame, or matrix}
  \item{configuration}{list of (probably part of) of
  configuration options. See diffLogoTableConfiguration.}
  \item{...}{set of parameters passed to the function 'axis'
  for plotting}
}
\description{
  Draws a table of DiffLogos.
}
\examples{
motif_folder= "extdata/pwm"
motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
motifs = list()
for (name in motif_names) {
  fileName = paste(motif_folder,"/",name,".pwm",sep="")
  file = system.file(fileName, package = "DiffLogo")
  motifs[[name]] = getPwmFromPwmFile(file)
}

diffLogoTable(motifs)
}
\author{
  Martin Nettling
}

