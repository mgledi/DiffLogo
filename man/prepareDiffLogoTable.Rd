% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffSeqLogo.R
\name{prepareDiffLogoTable}
\alias{prepareDiffLogoTable}
\title{Prepare a table of difflogo objects}
\usage{
prepareDiffLogoTable(PWMs, alphabet = DNA, configuration = list())
}
\arguments{
\item{PWMs}{a list/vector of position weight matrices (PWMs) each of type pwm, data.frame, or matrix}

\item{alphabet}{the alphabet of the given PWMs}

\item{configuration}{list of (probably part of) of configuration options. See diffLogoTableConfiguration.}
}
\value{
matrix of difference logos
}
\description{
Prepares a DiffLogoTable and generates an object that contains the hirarchical clustering and a matrix of prepared difference logos.
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
sampleSizes <- c(100, 150, 200, 250)

diffLogoTableObj = prepareDiffLogoTable(motifs);
}
\author{
Martin Nettling
}
