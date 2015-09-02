\name{diffLogoTable}
\alias{diffLogoTable}
\title{Draw DiffLogo-table}
\usage{
  diffLogoTable(PWMs, names2, stackHeight = shannonDivergence,
    baseDistribution = normalizedDifferenceOfProbabilities,
    uniformYaxis = TRUE, sparse = TRUE,
    showSequenceLogosTop = TRUE, enableClustering = TRUE,
    treeHeight = 0.5, margin = 0.02, ratio = 1,
    alphabet = DNA, ...)
}
\arguments{
  \item{PWMs}{a list/vector of position weight matrices (PWMs)
  each of type pwm, data.frame, or matrix} \item{stackHeight}
  {function for the height of a stack at position i}
  \item{baseDistribution}{function for the heights of
  the individual bases} \item{uniformYaxis}{if TRUE each
  DiffLogo is plotted with the same scaling of the y-axis}
  \item{sparse}{if TRUE margins are reduced and tickmarks are
  removed from the logo} \item{showSequenceLogosTop}{if TRUE
  the classical sequence logos are drawn above each column
  of the table} \item{treeHeight}{the height of the plotted
  cluster tree above the columns of the table; set equal
  to zero to omit the cluster tree} \item{enableClustering}
  {if TRUE the motifs are reordered, so that similar motifs
  have a small vertical and horizontal distance in the table}
  \item{margin}{the space reseverved for labels} \item{ratio}
  {the ratio of the plot; this is needed to determine the
  margin sizes correctly} \item{alphabet}{of type Alphabet}
  \item{...}{set of parameters passed to the function 'axis'
  for plotting}
}
\description{
  Draws a table of DiffLogos.
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

diffLogoTable(motifs)
}
\author{
  Martin Nettling
}

