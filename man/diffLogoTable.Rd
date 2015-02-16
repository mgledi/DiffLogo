\name{diffLogoTable}
\alias{diffLogoTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw DiffLogo-table
}
\description{
Draws a table of DiffLogos.
}
\usage{
diffLogoTable(PWMs, stackHeight = shannonDivergence, baseDistribution = normalizedDifferenceOfProbabilities, uniformYaxis = T, sparse = TRUE, showSequenceLogosTop = TRUE, treeHeight = 0.5, margin = 0.03, ratio = 16/10, alphabet = DNA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PWMs}{
a list/vector of position weight matrices (PWMs) each of type pwm, data.frame, or matrix
}
  \item{stackHeight}{
function for the height of a stack at position i
}
  \item{baseDistribution}{
function for the heights of the individual bases
}
  \item{uniformYaxis}{
if TRUE each DiffLogo is plotted with the same scaling of the y-axis
}
  \item{sparse}{
if TRUE margins are reduced and tickmarks are removed from the logo
}
  \item{showSequenceLogosTop}{
if TRUE the classical sequence logos are drawn above each column of the table
}
  \item{treeHeight}{
the height of the plotted cluster tree above the columns of the table; set equal to zero to omit the cluster tree
}
  \item{margin}{
the space reseverved for labels
}
  \item{ratio}{
the ratio of the plot; this is needed to determine the margin sizes correctly
}
  \item{alphabet}{
of type Alphabet
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Martin Gleditzsch
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (PWMs, stackHeight = shannonDivergence, baseDistribution = normalizedDifferenceOfProbabilities, 
    uniformYaxis = T, sparse = TRUE, showSequenceLogosTop = TRUE, 
    treeHeight = 0.5, margin = 0.03, ratio = 16/10, alphabet = DNA, 
    ...) 
{
    plot.new()
    dim = length(PWMs)
    st = 0
    if (showSequenceLogosTop) {
        st = 0.5
    }
    marDiffLogo = marSeqLogo = c(1, 1.5, 0.1, 0.1)
    if (sparse) {
        marDiffLogo = c(0.3, 1.2, 0.1, 0.1)
        marSeqLogo = c(0.3, 1.2, 0, 0.1)
    }
    similarities = matrix(0, dim, dim)
    diffLogos = list()
    palette = colorRampPalette(c(rgb(0.9, 1, 0.9), rgb(1, 0.9, 
        0.9)))(100)
    names = names(PWMs)
    ymin = 0
    ymax = 0
    for (i in 1:dim) {
        motif_i = names[i]
        for (k in 1:dim) {
            motif_k = names[k]
            similarities[i, k] = NA
            if (i != k) {
                diffLogoObj = createDiffLogoObject(PWMs[[motif_i]], 
                  PWMs[[motif_k]], stackHeight = stackHeight, 
                  baseDistribution = baseDistribution, alphabet = alphabet)
                if (uniformYaxis) {
                  ymin = min(diffLogoObj$ylim.negMax, ymin)
                  ymax = max(diffLogoObj$ylim.posMax, ymax)
                }
                similarities[i, k] = diffLogoObj$distance
            }
        }
    }
    colors = matrix(palette[cut(similarities, 100)], dim, dim)
    hc = hclust(dist(similarities))
    leaveOrder = hc$order
    dimV = c(dim, dim, dim + st + treeHeight, dim + st + treeHeight)
    for (i in 1:dim) {
        motif_i = names[leaveOrder[i]]
        for (k in 1:dim) {
            if (i != k) {
                motif_k = names[leaveOrder[k]]
                subplotcoords = c(i - 1, i, dim - k, dim - k + 
                  1)
                print(paste("plotting ", motif_i, " and ", motif_k))
                par(fig = (subplotcoords/dimV) * c(1 - margin, 
                  1 - margin, 1 - margin * ratio, 1 - margin * 
                    ratio) + c(margin, margin, 0, 0), new = TRUE, 
                  mar = c(0, 0, 0, 0))
                plot(NA, ylim = c(0, 1), xlim = c(0, 1), xaxt = "n", 
                  yaxt = "n", xaxs = "i", yaxs = "i", bty = "n")
                rect(0, 0, 1, 1, col = colors[leaveOrder[i], 
                  leaveOrder[k]], border = NA)
                par(fig = (subplotcoords/dimV) * c(1 - margin, 
                  1 - margin, 1 - margin * ratio, 1 - margin * 
                    ratio) + c(margin, margin, 0, 0), new = TRUE, 
                  mar = marDiffLogo)
                diffLogoObj = createDiffLogoObject(PWMs[[motif_i]], 
                  PWMs[[motif_k]], stackHeight = stackHeight, 
                  baseDistribution = baseDistribution, alphabet = alphabet)
                diffLogo(diffLogoObj, sparse = sparse, ymin = ymin, 
                  ymax = ymax)
            }
        }
        if (showSequenceLogosTop) {
            subplotcoords = c(i - 1, i, dim, dim + st)
            par(fig = (subplotcoords/dimV) * c(1 - margin, 1 - 
                margin, 1 - margin * ratio, 1 - margin * ratio) + 
                c(margin, margin, margin * ratio, margin * ratio), 
                new = TRUE, mar = marSeqLogo)
            seqLogo(PWMs[[motif_i]], sparse = sparse)
        }
    }
    if (treeHeight > 0) {
        par(fig = c(0 + 0.5, dim - 0.5, dim + st, dim + st + 
            treeHeight)/dimV * c(1 - margin, 1 - margin, 1 - 
            margin * ratio, 1 - margin * ratio) + c(margin, margin, 
            margin * ratio, margin * ratio), new = TRUE, mar = c(0, 
            0, 0.1, 0))
        plot(hc, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
            bty = "n", labels = rep("", dim), main = "")
    }
    par(fig = (c(0, dim, 0, dim)/dimV) * c(1 - margin, 1 - margin, 
        1 - margin * ratio, 1 - margin * ratio) + c(margin, margin, 
        0, 0), new = TRUE, mar = c(0, 0, 0, 0))
    plot(NA, ylim = c(0, dim), xlim = c(0, dim), xaxs = "i", 
        xaxt = "n", yaxt = "n", yaxs = "i", bty = "n", mar = c(0, 
            0, 0, 0))
    axis(2, pos = 0, at = (1:dim) - 0.5, labels = rev(names[leaveOrder]), 
        tick = F, mgp = c(3, 0, 0), ...)
    axis(3, pos = dim, at = (1:dim) - 0.5, labels = names[leaveOrder], 
        tick = F, mgp = c(3, 0, 0), ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
