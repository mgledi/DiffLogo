\name{diffLogoTable}
\alias{diffLogoTable}
\title{Draw DiffLogo-table}
\usage{
  diffLogoTable(PWMs, stackHeight = shannonDivergence,
    baseDistribution = normalizedDifferenceOfProbabilities,
    uniformYaxis = TRUE, sparse = TRUE,
    showSequenceLogosTop = TRUE, enableClustering = TRUE,
    treeHeight = 0.5, margin = 0.02, ratio = 1,
    alphabet = DNA, ...)
}
\arguments{
  \item{PWMs}{a list/vector of position weight matrices
  (PWMs) each of type pwm, data.frame, or matrix}

  \item{stackHeight}{function for the height of a stack at
  position i}

  \item{baseDistribution}{function for the heights of the
  individual bases}

  \item{uniformYaxis}{if TRUE each DiffLogo is plotted with
  the same scaling of the y-axis}

  \item{sparse}{if TRUE margins are reduced and tickmarks
  are removed from the logo}

  \item{showSequenceLogosTop}{if TRUE the classical
  sequence logos are drawn above each column of the table}

  \item{treeHeight}{the height of the plotted cluster tree
  above the columns of the table; set equal to zero to omit
  the cluster tree}

  \item{enableClustering}{if TRUE the motifs are reordered,
  so that similar motifs have a small vertical and
  horizontal distance in the table}

  \item{margin}{the space reseverved for labels}

  \item{ratio}{the ratio of the plot; this is needed to
  determine the margin sizes correctly}

  \item{alphabet}{of type Alphabet}

  \item{...}{TODO}
}
\description{
  Draws a table of DiffLogos.
}
\author{
  Martin Nettling
}

