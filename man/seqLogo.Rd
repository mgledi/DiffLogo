\name{seqLogo}
\alias{seqLogo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw sequence logo
}
\description{
Draws the classic sequence logo. 
}
\usage{
seqLogo(pwm, sparse = FALSE, drawLines = 0.5, stackHeight = informationContent, baseDistribution = probabilities, alphabet = DNA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pwm}{
representation of a position weight matrix (PWM) of type pwm, data.frame, or matrix
}
  \item{sparse}{
if TRUE margins are reduced and tickmarks are removed from the logo
}
  \item{drawLines}{
distance between background lines
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
function (pwm, sparse=FALSE, drawLines=0.5, stackHeight=informationContent, baseDistribution=probabilities, alphabet=DNA) { 
    pwm = preconditionTransformPWM(pwm,alphabet);
    preconditionPWM(pwm);

    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm)

    ylim.negMax = 0;
    ylim.posMax = 0;
   
    wt = 1.0
    x.pos = 0.5 # initial position on x axis is 0.5; Letter is one right from this point
    eps = 0; # spacer between two bases in one stack
    heights = c(); ymins=c(); ymaxs=c()
    for (j in 1:npos) {
        column = pwm[, j]
        sh = stackHeight(column);
        hts = baseDistribution(column) * sh$height ;
        letterOrder = order(abs(hts)) # reorder letters
        yneg.pos = 0 
        ypos.pos = 0
        for (i in 1:alphabet$size) {
            ht = hts[letterOrder[i]]
            y.pos = ypos.pos;
            ht = ht - eps;
            ypos.pos = ypos.pos + ht + eps
            char = alphabet$chars[letterOrder[i]]
            col = alphabet$cols[letterOrder[i]];
            letters = addLetter(letters, letterPolygons[[char]], x.pos, y.pos, ht, wt*0.99, col=col)
        }
        x.pos = x.pos + wt
    }
    if(sparse) {
        plot(NA, xlim=c(0.5,x.pos), ylim=c(0,log2(alphabet$size)),xaxt="n", ylab="",
        mgp=c(0, .35, 0),tck=-0.02, cex.axis=0.8, frame.plot=F,xlab="")
    } else {
        plot(NA, xlim=c(0.5,x.pos), ylim=c(0,log2(alphabet$size)), xaxt="n", ylab=sh$ylab, frame.plot=F,xlab="Position")
    }

    for(y in seq(0,log2(alphabet$size),drawLines)) {
        abline(y,0,col="gray");
    }
    
    if(sparse) {
        axis(1,labels=c("",rep("",npos),""), at=c(0,1:npos,npos+1),tck=-0.02)
    } else {
        axis(1,labels=c("",1:npos,""),at=c(0,1:npos,npos+1))
    }
    polygon(letters, col=letters$col, border=NA)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
