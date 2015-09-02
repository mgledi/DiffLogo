\name{diffLogoFromPwm}
\alias{diffLogoFromPwm}
\title{Draw DiffLogo from PWM}
\usage{
  diffLogoFromPwm(pwm1, pwm2, ymin = 0,
    ymax = 0, stackHeight = shannonDivergence,
    baseDistribution = normalizedDifferenceOfProbabilities,
    sparse = FALSE, alphabet = DNA)
}
\arguments{
  \item{pwm1}{representation of the first position weight
  matrix (PWM) of type pwm, data.frame, or matrix} \item{pwm2}
  {representation of the second position weight matrix (PWM)
  of type pwm, data.frame, or matrix} \item{ymin}{minimum
  value on the y-axis} \item{ymax}{maximum value on the y-
  axis} \item{stackHeight}{function for the height of a
  stack at position i} \item{baseDistribution}{function for
  the heights of the individual bases} \item{sparse}{if TRUE
  margins are reduced and tickmarks are removed from the logo}
  \item{alphabet}{of type Alphabet}
}
\description{
  Draws the difference of two sequence logos.
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

diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)
}
\author{
  Martin Nettling
}

