##' Draws the classic sequence logo. 
##'
##' @title Draw sequence logo
##' @param pwm representation of a position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param sparse if TRUE margins are reduced and tickmarks are removed from the logo
##' @param drawLines distance between background lines
##' @param stackHeight function for the height of a stack at position i
##' @param baseDistribution function for the heights of the individual bases
##' @param alphabet of type Alphabet
##' @param main the main title for the plot
##' @return none (draws sequence logo)
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_name = "HepG2"
##' fileName = paste(motif_folder,"/",motif_name,".pwm",sep="")
##' file = system.file(fileName, package = "DiffLogo")
##' motif = getPwmFromPwmFile(file)
##' seqLogo(pwm = motif)
seqLogo = function (pwm, sparse=FALSE, drawLines=0.5, stackHeight=informationContent, baseDistribution=probabilities, alphabet=DNA, main=NULL) { 
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
        mgp=c(0, .35, 0),tck=-0.02, cex.axis=0.8, frame.plot=FALSE,xlab="", main=main)
    } else {
        plot(NA, xlim=c(0.5,x.pos), ylim=c(0,log2(alphabet$size)), xaxt="n", ylab=sh$ylab, frame.plot=FALSE,xlab="Position", main=main)
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


# appends the letter which to the object letters
addLetter = function (letters, letterPolygon, x.pos, y.pos, ht, wt, col="black") {
    x = x.pos + wt * letterPolygon$x
    y = y.pos + ht * letterPolygon$y
    polygons = sum(is.na(x))+1  # a letter can consist of more then one polygon
    letters$x = c(letters$x, NA, x)
    letters$y = c(letters$y, NA, y)
    letters$col = c(letters$col, rep(col, polygons))
    letters
}
