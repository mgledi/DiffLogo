\name{diffLogoFromPwm}
\alias{diffLogoFromPwm}
\title{Draw DiffLogo from PWM}
\usage{
  diffLogoFromPwm(pwm1, pwm2, ymin = 0, ymax = 0,
    stackHeight = shannonDivergence,
    baseDistribution = normalizedDifferenceOfProbabilities,
    sparse = FALSE, alphabet = DNA)
}
\arguments{
  \item{pwm1}{representation of the first position weight
  matrix (PWM) of type pwm, data.frame, or matrix}

  \item{pwm2}{representation of the second position weight
  matrix (PWM) of type pwm, data.frame, or matrix}

  \item{ymin}{minimum value on the y-axis}

  \item{ymax}{maximum value on the y-axis}

  \item{stackHeight}{function for the height of a stack at
  position i}

  \item{baseDistribution}{function for the heights of the
  individual bases}

  \item{sparse}{if TRUE margins are reduced and tickmarks
  are removed from the logo}

  \item{alphabet}{of type Alphabet}
}
\description{
  Draws the difference of two sequence logos.
}
\author{
  Martin Nettling
}

