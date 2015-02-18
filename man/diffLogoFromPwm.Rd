\name{diffLogoFromPwm}
\alias{diffLogoFromPwm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw DiffLogo from PWM
}
\description{
Draws the difference of two sequence logos. 
}
\usage{
diffLogoFromPwm(pwm1, pwm2, ymin = 0, ymax = 0, stackHeight = shannonDivergence, baseDistribution = normalizedDifferenceOfProbabilities, sparse = FALSE, alphabet = DNA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pwm1}{
representation of the first position weight matrix (PWM) of type pwm, data.frame, or matrix
}
  \item{pwm2}{
representation of the second position weight matrix (PWM) of type pwm, data.frame, or matrix
}
  \item{ymin}{
minimum value on the y-axis
}
  \item{ymax}{
maximum value on the y-axis
}
  \item{stackHeight}{
function for the height of a stack at position i
}
  \item{baseDistribution}{
function for the heights of the individual bases
}
  \item{sparse}{
if TRUE margins are reduced and tickmarks are removed from the logo
}
  \item{alphabet}{
of type Alphabet
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
function (pwm1, pwm2, ymin = 0, ymax = 0, stackHeight = shannonDivergence, 
    baseDistribution = normalizedDifferenceOfProbabilities, sparse = FALSE, 
    alphabet = DNA) 
{
    diffLogoObj = createDiffLogoObject(pwm1, pwm2, stackHeight = stackHeight, 
        baseDistribution = baseDistribution, alphabet = alphabet)
    diffLogo(diffLogoObj, ymin = ymin, ymax = ymax, sparse = sparse)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
