source("./R/preconditions.R"); # contains functions that check preconditions
source("./R/alphabet.R"); # contains functions for drawing
source("./R/seqLogo.R"); # contains functions for drawing
source("./R/stackHeights.R"); # contains functions to calculate the stackheihts in a diffLogo
source("./R/baseDistrs.R"); # contains functions to calculate the proportions for each base in a stack of a diffLogo

##' Creates a DiffLogo object
##'
##' @title DiffLogo object
##' @param pwm1 representation of the first position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param pwm2 representation of the second position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param stackHeight function for the height of a stack at position i
##' @param baseDistribution function for the heights of the individual bases
##' @param alphabet of type Alphabet
##' @return DiffLogo object
##' @export
##' @author Martin Nettling
createDiffLogoObject = function (pwm1, pwm2, stackHeight=shannonDivergence, baseDistribution=normalizedDifferenceOfProbabilities,alphabet=DNA) {
    pwm1 = preconditionTransformPWM(pwm1,alphabet);
    pwm2 = preconditionTransformPWM(pwm2,alphabet);
    preconditionPWM(pwm1);
    preconditionPWM(pwm2);
    preconditionPWMSameSize(pwm1,pwm2);

    # init needed variables
    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm1)
    eps = 0.0000; # spacer between two bases in one stack
    ylim.negMax = 0;
    ylim.posMax = 0;

    npos = ncol(pwm1)
    wt = 1.0 # the width of one letter
    x.pos = 0.5 # initial position on x axis is 0.5; Letter is one right from this point
    heights = c(); ymins=c(); ymaxs=c()

    # determine intermediate values
    for (j in 1:npos) {
	    heightObj = stackHeight(pwm1[,j],pwm2[,j]);
	    heights[j] = heightObj$height;
	    ylab = heightObj$ylab;
	    hts = heights[j] * baseDistribution(pwm1[,j],pwm2[,j]);
        letterOrder = order(abs(hts)) # reorder letters

	    yneg.pos = 0 
	    ypos.pos = 0
        # adds all letters as polygons to the list of letters
	    for (i in 1:length(hts)) {
	        ht = hts[letterOrder[i]]
	        if (ht >= 0){ 
	            y.pos = ypos.pos;
                ypos.pos = ypos.pos + ht + eps
	        } else if(ht < 0 ) {
	            y.pos = yneg.pos;
		        yneg.pos = yneg.pos + ht - eps
	        }
            char = alphabet$chars[ (letterOrder[i]-1)%%alphabet$size+1 ]
            col = alphabet$cols[ (letterOrder[i]-1)%%alphabet$size+1 ];
            letters = addLetter(letters, letterPolygons[[char]], x.pos, y.pos, ht, wt, col=col)
	    }
	
        x.pos = x.pos + wt
	    ylim.negMax = min(ylim.negMax, yneg.pos)
	    ylim.posMax = max(ylim.posMax, ypos.pos)
	    # remember values for plotting
	    ymins[j] = yneg.pos
	    ymaxs[j] = ypos.pos
    }

    diffObj = list();
    diffObj$ylab=ylab; # should be set by the stackHeight function
    diffObj$npos = npos;
    diffObj$heights = heights;
    diffObj$ymins = ymins
    diffObj$ymaxs = ymaxs
    diffObj$ylim.negMax = ylim.negMax
    diffObj$ylim.posMax = ylim.posMax
    diffObj$letters = letters;
    diffObj$pwm1 = pwm1
    diffObj$pwm2 = pwm2
    diffObj$distance = sum(heights) # TODO as function
    diffObj$alphabet = alphabet
    
    class(diffObj) = "DiffLogo"
    return(diffObj);
}

##' Draws the difference of two sequence logos. 
##'
##' @title Draw DiffLogo
##' @param diffLogoObj a DiffLogoObject created by the function createDiffLogoObject
##' @param ymin minimum value on the y-axis
##' @param ymax maximum value on the y-axis
##' @param sparse if TRUE margins are reduced and tickmarks are removed from the logo
##' @export
##' @author Martin Nettling
diffLogo = function (diffLogoObj, ymin=0, ymax=0, sparse=FALSE) {
    if(class(diffLogoObj) != "DiffLogo") {
        msg = paste("Expected DiffLogo, but got ", class(diffLogoObj), ". Use #createDiffLogoObject to get an DiffLogo from two PWMs.",sep="")
        stop(msg)
    }

    yAbsMax=2 # this variable defines the possible maximum and the minimum of the y-axis

    if(ymin == 0) {
        ymin = min(diffLogoObj$ylim.posMax*1.0,yAbsMax);
    }	    
    if(ymax == 0) {
        ymax = max(diffLogoObj$ylim.negMax*1.0,-yAbsMax);
    }

    # set ylab
     ylab = diffLogoObj$ylab

    if(sparse) {
        # the sparse plot has small ticks, small y-labels, no x-labels, no xlab, no ylab
        plot(NA, xlim=c(0.5,diffLogoObj$npos + 0.5), ylim=c(ymin,ymax), xaxt="n", ylab="", mgp=c(0, .35, 0), tck=-0.02, cex.axis=0.8, frame.plot=F,xlab="")
    } else {
        plot(NA, xlim=c(0.5,diffLogoObj$npos + 0.5), ylim=c(ymin,ymax), xaxt="n", ylab=ylab, xlab="Position", frame.plot=F, )
    }

    yLabs = c("","","");
    yAt = c(-yAbsMax,0,yAbsMax);


    if(sparse) {
        axis(1,labels=c(rep("",diffLogoObj$npos)), at=(1:diffLogoObj$npos),tck=-0.02)
        axis(2,labels=yLabs,at=yAt,mgp=c(0, .35, 0),tck=-0.02, cex.axis=0.8)
    } else {
        axis(1,labels=c(1:diffLogoObj$npos),at=(1:diffLogoObj$npos))
        axis(2,labels=c("",""),at=c(-yAbsMax,yAbsMax))
    }
    
    polygon(diffLogoObj$letters, col=diffLogoObj$letters$col, border=F)
    lines(c(0,diffLogoObj$npos), c(0,0) ) # the line at y = 0
}

##' Draws the difference of two sequence logos. 
##'
##' @title Draw DiffLogo from PWM
##' @param pwm1 representation of the first position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param pwm2 representation of the second position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param ymin minimum value on the y-axis
##' @param ymax maximum value on the y-axis
##' @param stackHeight function for the height of a stack at position i
##' @param baseDistribution function for the heights of the individual bases
##' @param sparse if TRUE margins are reduced and tickmarks are removed from the logo
##' @param alphabet of type Alphabet
##' @export
##' @author Martin Nettling
diffLogoFromPwm = function (pwm1, pwm2, ymin=0, ymax=0,stackHeight=shannonDivergence, baseDistribution=normalizedDifferenceOfProbabilities, sparse=FALSE, alphabet=DNA) {
    diffLogoObj = createDiffLogoObject(pwm1,pwm2,stackHeight=stackHeight, baseDistribution=baseDistribution, alphabet=alphabet);
    diffLogo(diffLogoObj,ymin=ymin, ymax=ymax, sparse=sparse)
}

##' Draws a table of DiffLogos.
##'
##' @title Draw DiffLogo-table
##' @param PWMs a list/vector of position weight matrices (PWMs) each of type pwm, data.frame, or matrix
##' @param stackHeight function for the height of a stack at position i
##' @param baseDistribution function for the heights of the individual bases
##' @param uniformYaxis if TRUE each DiffLogo is plotted with the same scaling of the y-axis
##' @param sparse if TRUE margins are reduced and tickmarks are removed from the logo
##' @param showSequenceLogosTop if TRUE the classical sequence logos are drawn above each column of the table
##' @param treeHeight the height of the plotted cluster tree above the columns of the table; set equal to zero to omit the cluster tree
##' @param margin TODO
##' @param ratio TODO
##' @param alphabet of type Alphabet
##' @export
##' @author Martin Nettling
diffLogoTable = function (PWMs, stackHeight=shannonDivergence, baseDistribution=normalizedDifferenceOfProbabilities, uniformYaxis=T, sparse=TRUE, showSequenceLogosTop=TRUE, treeHeight=0.5, margin=0.03, ratio=16/10, alphabet=DNA,...) {
    plot.new();
    dim = length(PWMs);

    st = 0;
    if ( showSequenceLogosTop ) {
        st = 0.5;
    }

    marDiffLogo = marSeqLogo = c(1,1.5,0.1,0.1);
    if(sparse) {
        marDiffLogo = c(0.3,1.2,0.1,0.1);
        marSeqLogo = c(0.3,1.2,0.0,0.1);
    }

    similarities = matrix(0,dim,dim);
    diffLogos = list();
    palette = colorRampPalette(c(rgb(0.9,1,0.9),rgb(1,0.9,0.9)))(100)
    names = names(PWMs);
    ymin = 0;
    ymax = 0;
    for ( i in 1:dim) {
        motif_i = names[i];
        for ( k in 1:dim) {
            motif_k = names[k];
	        similarities[i,k] = NA
            if( i != k ) {
		        diffLogoObj = createDiffLogoObject(PWMs[[ motif_i ]],PWMs[[ motif_k ]],stackHeight=stackHeight, baseDistribution=baseDistribution, alphabet=alphabet);
                if(uniformYaxis) {
		            ymin = min(diffLogoObj$ylim.negMax,ymin)
		            ymax = max(diffLogoObj$ylim.posMax,ymax)
		        }
		        similarities[i,k] = diffLogoObj$distance;
	        }
        }
    }
    colors = matrix(palette[cut(similarities,100)],dim,dim)
    hc = hclust(dist(similarities));
    leaveOrder = hc$order;

    # draw DiffLogos
    dimV = c(dim, dim,dim+st+treeHeight,dim+st+treeHeight);
    for ( i in 1:dim) {
        motif_i = names[leaveOrder[i]];
        for ( k in 1:dim) {
            if( i != k ) {
                motif_k = names[leaveOrder[k]];
                subplotcoords = c(i-1,i,dim-k,dim-k+1)
                
                print(paste("plotting ",motif_i," and ",motif_k));
                par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,0,0), new=TRUE, mar=c(0,0,0,0))
		
		        plot(NA,ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
        		rect(0,0,1,1,col=colors[leaveOrder[i],leaveOrder[k]],border=NA);
		
                par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,0,0), new=TRUE, mar=marDiffLogo)
                diffLogoObj = createDiffLogoObject(PWMs[[ motif_i ]],PWMs[[ motif_k ]],stackHeight=stackHeight, baseDistribution=baseDistribution, alphabet=alphabet); 
                diffLogo(diffLogoObj,sparse=sparse,ymin=ymin,ymax=ymax)
            }
        }
        if(showSequenceLogosTop) {
            subplotcoords = c(i-1, i, dim, dim + st)
            #par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=c(0,0,0,0)) 
            #plot(NA,ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
            #rect(0,0,1,1,col="gray",border=NA);            
            par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=marSeqLogo)       
            seqLogo(PWMs[[ motif_i ]],sparse=sparse)
        }
    }


    if(treeHeight > 0) {
        par(fig=c(0 + 0.5, dim - 0.5, dim + st, dim + st + treeHeight) / dimV * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=c(.0,.0,.1,.0))
        plot(hc, xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n", labels=rep("",dim), main="" )            
        #plot(NA,ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
        #rect(0,0,1,1,col="gray",border=NA);    
    }

    # add names
    par(fig=(c(0,dim,0,dim) / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio)+ c(margin,margin,0,0), new=TRUE, mar=c(0,0,0,0))

    plot(NA,ylim=c(0,dim),xlim=c(0,dim),xaxs="i",xaxt="n",yaxt="n",yaxs="i", bty="n", mar=c(0,0,0,0)) 
    axis(2, pos=0, at= (1:dim) - 0.5, labels = rev(names[leaveOrder]), tick = F, mgp = c(3, 0, 0), ...)
    axis(3, pos=dim, at= (1:dim) - 0.5, labels = names[leaveOrder], tick = F, mgp = c(3, 0, 0), ...)
}
