\docType{data}
\name{DNA}
\alias{DNA}
\title{DNA alphabet}
\format{List of 3
 $ chars: chr [1:4] "A" "C" "G" "T"
 $ cols : chr [1:4] "green4" "blue" "orange" "red"
 $ size : int 4
 - attr(*, "class")= chr "Alphabet"}
\usage{
  DNA
}
\description{
  the DNA alphabet, i.e. A, C, G, T
}
\examples{
motif_folder= "extdata/pwm"
motif_name = "HepG2"
fileName = paste(motif_folder,"/",motif_name,".txt",sep="")
file = system.file(fileName, package = "DiffLogo")
motif = as.matrix(read.delim(file,header=FALSE))
seqLogo(pwm = motif, alphabet=DNA)
}
\author{
  Martin Nettling
}
\keyword{datasets}

