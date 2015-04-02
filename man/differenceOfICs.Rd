\name{differenceOfICs}
\alias{differenceOfICs}
\title{normalized information content differences}
\usage{
  differenceOfICs(p1, p2)
}
\arguments{
  \item{p1}{probability vector representing the first
  symbol distribution}

  \item{p2}{probability vector representing the second
  symbol distribution}
}
\value{
  a vector with one result for each symbol
}
\description{
  information content differences normalized by the sum of
  absolute information content differences for the given
  pair of probability vectors
}
\author{
  Martin Nettling
}

