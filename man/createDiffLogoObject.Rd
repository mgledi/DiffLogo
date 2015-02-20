\name{createDiffLogoObject}
\alias{createDiffLogoObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DiffLogo object
}
\description{
Creates a DiffLogo object
}
\usage{
createDiffLogoObject(pwm1, pwm2, stackHeight = shannonDivergence, baseDistribution = normalizedDifferenceOfProbabilities, alphabet = DNA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pwm1}{
representation of the first position weight matrix (PWM) of type pwm, data.frame, or matrix
}
  \item{pwm2}{
representation of the second position weight matrix (PWM) of type pwm, data.frame, or matrix
}
  \item{stackHeight}{
function for the height of a stack at position i
}
  \item{baseDistribution}{
function for the heights of the individual bases
}
  \item{alphabet}{
of type Alphabet
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
DiffLogo object
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Martin Nettling
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
function (pwm1, pwm2, stackHeight = shannonDivergence, baseDistribution = normalizedDifferenceOfProbabilities, 
    alphabet = DNA) 
{
    pwm1 = preconditionTransformPWM(pwm1, alphabet)
    pwm2 = preconditionTransformPWM(pwm2, alphabet)
    preconditionPWM(pwm1)
    preconditionPWM(pwm2)
    preconditionPWMSameSize(pwm1, pwm2)
    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm1)
    eps = 5e-04
    ylim.negMax = 0
    ylim.posMax = 0
    npos = ncol(pwm1)
    wt = 1
    x.pos = 0.5
    heights = c()
    ymins = c()
    ymaxs = c()
    for (j in 1:npos) {
        heightObj = stackHeight(pwm1[, j], pwm2[, j])
        heights[j] = heightObj$height
        ylab = heightObj$ylab
        hts = heights[j] * baseDistribution(pwm1[, j], pwm2[, 
            j])
        letterOrder = order(abs(hts))
        yneg.pos = 0
        ypos.pos = 0
        for (i in 1:length(hts)) {
            ht = hts[letterOrder[i]]
            if (ht >= 0) {
                y.pos = ypos.pos
                ypos.pos = ypos.pos + ht + eps
            }
            else if (ht < 0) {
                y.pos = yneg.pos
                yneg.pos = yneg.pos + ht - eps
            }
            char = alphabet$chars[(letterOrder[i] - 1)\%\%alphabet$size + 
                1]
            col = alphabet$cols[(letterOrder[i] - 1)\%\%alphabet$size + 
                1]
            letters = addLetter(letters, letterPolygons[[char]], 
                x.pos, y.pos, ht, wt, col = col)
        }
        x.pos = x.pos + wt
        ylim.negMax = min(ylim.negMax, yneg.pos)
        ylim.posMax = max(ylim.posMax, ypos.pos)
        ymins[j] = yneg.pos
        ymaxs[j] = ypos.pos
    }
    diffObj = list()
    diffObj$ylab = ylab
    diffObj$npos = npos
    diffObj$heights = heights
    diffObj$ymins = ymins
    diffObj$ymaxs = ymaxs
    diffObj$ylim.negMax = ylim.negMax
    diffObj$ylim.posMax = ylim.posMax
    diffObj$letters = letters
    diffObj$pwm1 = pwm1
    diffObj$pwm2 = pwm2
    diffObj$distance = sum(heights)
    diffObj$alphabet = alphabet
    class(diffObj) = "DiffLogo"
    return(diffObj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
