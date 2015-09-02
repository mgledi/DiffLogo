\docType{data}
\name{ASN}
\alias{ASN}
\title{ASN alphabet}
\format{List of 3
 $ chars: chr [1:20] "A" "C" "D" "E" ...
 $ cols : chr [1:20] "Light green" "Green" "Dark green" "Dark green" ...
 $ size : int 20
 - attr(*, "class")= chr "Alphabet"}
\usage{
  ASN
}
\description{
  the amino acid alphabet (20 symbols), i.e. A, C, D, E, F, G,
  H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
}
\examples{
motif_folder= "extdata/alignments"
motif_name = "calamodulin_1"
fileName = paste(motif_folder,"/",motif_name,".txt",sep="")
file = system.file(fileName, package = "DiffLogo")
motif = getPwmFromAlignment(readLines(file), ASN, 1)
seqLogo(pwm = motif, alphabet=ASN)
}
\author{
  Martin Nettling
}
\keyword{datasets}

