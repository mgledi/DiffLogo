\name{createDiffLogoObject}
\alias{createDiffLogoObject}
\title{DiffLogo object}
\usage{
  createDiffLogoObject(pwm1, pwm2,
    stackHeight = shannonDivergence,
    baseDistribution = normalizedDifferenceOfProbabilities,
    alphabet = DNA)
}
\arguments{
  \item{pwm1}{representation of the first position weight
  matrix (PWM) of type pwm, data.frame, or matrix} \item{pwm2}
  {representation of the second position weight matrix (PWM)
  of type pwm, data.frame, or matrix} \item{stackHeight}
  {function for the height of a stack at position i}
  \item{baseDistribution}{function for the heights of the
  individual bases} \item{alphabet}{of type Alphabet}
}
\value{
  DiffLogo object
}
\description{
  Creates a DiffLogo object
}
\examples{
motif_folder= "extdata/pwm"
motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
motifs = list()
for (name in motif_names) {
  fileName = paste(motif_folder,"/",name,".txt",sep="")
  file = system.file(fileName, package = "DiffLogo")
  motifs[[name]] = as.matrix(read.delim(file,header=FALSE))
}

pwm1 = motifs[[motif_names[[1]]]]
pwm2 = motifs[[motif_names[[2]]]]

diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
diffLogo(diffLogoObj)
}
\author{
  Martin Nettling
}

