\name{informationContent}
\alias{informationContent}
\title{information content}
\usage{
  informationContent(p)
}
\arguments{
  \item{p}{probability vector representing the symbol
  distribution}
}
\value{
  an object consisting of height a ylab
}
\description{
  the information content for the given probability vector
}
\examples{
motif_folder= "extdata/pwm"
motif_name = "HepG2"
fileName = paste(motif_folder,"/",motif_name,".txt",sep="")
file = system.file(fileName, package = "DiffLogo")
motif = as.matrix(read.delim(file,header=FALSE))
seqLogo(pwm = motif, stackHeight = informationContent)
}
\author{
  Martin Nettling
}

