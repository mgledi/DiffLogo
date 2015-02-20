\name{getPwmFromAlignment}
\alias{getPwmFromAlignment}
\title{Create PWM from alignment}
\usage{
  getPwmFromAlignment(alignment, alphabet, pseudoCount)
}
\arguments{
  \item{alignment}{a vector or list of sequences each with
  equal length}

  \item{alphabet}{of type Alphabet}

  \item{pseudoCount}{the number of pseudo-observations for
  each character in the alphabet}
}
\value{
  PWM as matrix
}
\description{
  Creates a matrix-representation of a PWM from a set of
  sequences
}
\author{
  Hendrik Treutler
}

