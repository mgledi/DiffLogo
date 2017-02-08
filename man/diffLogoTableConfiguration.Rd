\name{diffLogoTableConfiguration}
\alias{diffLogoTableConfiguration}
\title{Configuration object for diffLogoTable}
\usage{
  diffLogoTableConfiguration(alphabet,
    stackHeight = shannonDivergence,
    baseDistribution = normalizedDifferenceOfProbabilities,
    uniformYaxis = TRUE, sparse = TRUE,
    showSequenceLogosTop = TRUE, enableClustering = TRUE,
    treeHeight = 0.5, margin = 0.02, ratio = 1,
    align_pwms = F, multiple_align_pwms = T,
    unaligned_penalty = divergencePenaltyForUnaligned,
    try_reverse_complement = T, length_normalization = F)
}
\arguments{
  \item{stackHeight}{function for the height of a stack at
  position i} \item{baseDistribution}{function for the heights
  of the individual bases} \item{uniformYaxis}{if TRUE each
  DiffLogo is plotted with the same scaling of the y-axis}
  \item{sparse}{if TRUE margins are reduced and tickmarks are
  removed from the logo} \item{showSequenceLogosTop}{if TRUE
  the classical sequence logos are drawn above each column
  of the table} \item{treeHeight}{the height of the plotted
  cluster tree above the columns of the table; set equal
  to zero to omit the cluster tree} \item{enableClustering}
  {if TRUE the motifs are reordered, so that similar motifs
  have a small vertical and horizontal distance in the
  table} \item{margin}{the space reseverved for labels}
  \item{ratio}{the ratio of the plot; this is needed to
  determine the margin sizes correctly} \item{alphabet}
  {of type Alphabet} \item{align_pwms}{if True, will
  align and extend pwms in each sell of diffLogoTable
  independently.} \item{unaligned_penalty}{is a function
  for localPwmAlignment.} \item{try_reverse_complement}
  {if True, alignment will try reverse complement pwms}
  \item{length_normalization}{if True, divergence between pwms
  is divided by length of pwms.} \item{...}{set of parameters
  passed to the function 'axis' for plotting}
}
\description{
  Default configuration list for diffLogoTable
}
\examples{
diffLogoTableConfiguration(DNA)
}
\author{
  Lando Andrey
}

