\name{diffLogo}
\alias{diffLogo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw DiffLogo
}
\description{
Draws the difference of two sequence logos. 
}
\usage{
diffLogo(diffLogoObj, ymin = 0, ymax = 0, sparse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diffLogoObj}{
a DiffLogoObject created by the function createDiffLogoObject
}
  \item{ymin}{
minimum value on the y-axis
}
  \item{ymax}{
maximum value on the y-axis
}
  \item{sparse}{
if TRUE margins are reduced and tickmarks are removed from the logo
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
function (diffLogoObj, ymin = 0, ymax = 0, sparse = FALSE) 
{
    if (class(diffLogoObj) != "DiffLogo") {
        msg = paste("Expected DiffLogo, but got ", class(diffLogoObj), 
            ". Use #createDiffLogoObject to get an DiffLogo from two PWMs.", 
            sep = "")
        stop(msg)
    }
    yAbsMax = 2
    if (ymin == 0) {
        ymin = min(diffLogoObj$ylim.posMax * 1, yAbsMax)
    }
    if (ymax == 0) {
        ymax = max(diffLogoObj$ylim.negMax * 1, -yAbsMax)
    }
    ylab = diffLogoObj$ylab
    if (sparse) {
        plot(NA, xlim = c(0.5, diffLogoObj$npos + 0.5), ylim = c(ymin, 
            ymax), xaxt = "n", ylab = "", mgp = c(0, 0.35, 0), 
            tck = -0.02, cex.axis = 0.8, frame.plot = F, xlab = "")
    }
    else {
        plot(NA, xlim = c(0.5, diffLogoObj$npos + 0.5), ylim = c(ymin, 
            ymax), xaxt = "n", ylab = ylab, xlab = "Position", 
            frame.plot = F, )
    }
    yLabs = c("", "", "")
    yAt = c(-yAbsMax, 0, yAbsMax)
    if (sparse) {
        axis(1, labels = c(rep("", diffLogoObj$npos)), at = (1:diffLogoObj$npos), 
            tck = -0.02)
        axis(2, labels = yLabs, at = yAt, mgp = c(0, 0.35, 0), 
            tck = -0.02, cex.axis = 0.8)
    }
    else {
        axis(1, labels = c(1:diffLogoObj$npos), at = (1:diffLogoObj$npos))
        axis(2, labels = c("", ""), at = c(-yAbsMax, yAbsMax))
    }
    polygon(diffLogoObj$letters, col = diffLogoObj$letters$col, 
        border = F)
    lines(c(0, diffLogoObj$npos), c(0, 0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
