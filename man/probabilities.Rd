\name{probabilities}
\alias{probabilities}
\title{probabilities}
\usage{
  probabilities(p)
}
\arguments{
  \item{p}{probability vector representing the symbol
  distribution}
}
\value{
  the given vector
}
\description{
  the given probabilities
}
\examples{
motif_folder= "extdata/pwm"
motif_name = "HepG2"
fileName = paste(motif_folder,"/",motif_name,".txt",sep="")
file = system.file(fileName, package = "DiffLogo")
motif = as.matrix(read.delim(file,header=FALSE))
seqLogo(pwm = motif, baseDistribution = probabilities)
}
\author{
  Martin Nettling
}

