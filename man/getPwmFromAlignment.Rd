\name{getPwmFromAlignment}
\alias{getPwmFromAlignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create PWM from alignment
}
\description{
Creates a matrix-representation of a PWM from a set of sequences
}
\usage{
getPwmFromAlignment(alignment, alphabet, pseudoCount)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{
a vector or list of sequences each with equal length
}
  \item{alphabet}{
of type Alphabet
}
  \item{pseudoCount}{
the number of pseudo-observations for each character in the alphabet
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
PWM as matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hendrik Treutler
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
function (alignment, alphabet, pseudoCount) 
{
    alphabetSize = alphabet$size
    alignmentLength = nchar(alignment[[1]])
    numberOfSequences = length(alignment)
    pwm = matrix(nrow = alphabetSize, ncol = alignmentLength)
    colnames(pwm) = 1:alignmentLength
    rownames(pwm) = alphabet$chars
    for (charIdx in 1:alphabetSize) {
        for (posIdx in 1:alignmentLength) {
            pwm[charIdx, posIdx] = (length(which(substr(alignment, 
                posIdx, posIdx) == alphabet$chars[[charIdx]])) + 
                pseudoCount)/(numberOfSequences + pseudoCount * 
                alphabetSize)
        }
    }
    return(pwm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
