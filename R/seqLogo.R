##' Draws the classic sequence logo. 
##'
##' @title Draw sequence logo
##' @param pwm representation of a position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param sparse if TRUE margins are reduced and tickmarks are removed from the logo
##' @param drawLines a vector of y-values where to draw horizontal gray lines
##' @param stackHeight function for the height of a stack at position i
##' @param baseDistribution function for the heights of the individual bases
##' @param alphabet of type Alphabet
##' @export
##' @author Martin Nettling
seqLogo = function (pwm, sparse=FALSE, drawLines=c(0.5,1.0,1.5,2.0), stackHeight=informationContent, baseDistribution=probabilities, alphabet=DNA) { 
    pwm = preconditionTransformPWM(pwm,alphabet);
    preconditionPWM(pwm);

    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm)

    ylim.negMax = 0;
    ylim.posMax = 0;
   
    wt = 1.0
    x.pos = 0.5 # initial position on x axis is 0.5; Letter is one right from this point
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
            ypos.pos = ypos.pos + ht + 0.0005
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

    for(y in drawLines) {
        abline(y,0,col="gray");
    }
    
    if(sparse) {
        axis(1,labels=c("",rep("",npos),""), at=c(0,1:npos,npos+1),tck=-0.02)
    } else {
        axis(1,labels=c("",1:npos,""),at=c(0,1:npos,npos+1))
    }
    polygon(letters, col=letters$col, border=NA)
}


# appends the letter which to the object letters
addLetter = function (letters, letterPolygon, x.pos, y.pos, ht, wt, col="black") 
{
    x = x.pos + wt * letterPolygon$x
    y = y.pos + ht * letterPolygon$y
    polygons = sum(is.na(x))+1  # a letter can consist of more then one polygon
    letters$x = c(letters$x, NA, x)
    letters$y = c(letters$y, NA, y)
    letters$col = c(letters$col, rep(col,polygons))
    letters
}
