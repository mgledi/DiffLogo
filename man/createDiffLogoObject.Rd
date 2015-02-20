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
  matrix (PWM) of type pwm, data.frame, or matrix}

  \item{pwm2}{representation of the second position weight
  matrix (PWM) of type pwm, data.frame, or matrix}

  \item{stackHeight}{function for the height of a stack at
  position i}

  \item{baseDistribution}{function for the heights of the
  individual bases}

  \item{alphabet}{of type Alphabet}
}
\value{
  DiffLogo object
}
\description{
  Creates a DiffLogo object
}
\author{
  Martin Nettling
}

