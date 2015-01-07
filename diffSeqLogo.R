# draws the letter A at x=x.pos and y=y.pow with height=ht and width=wt
letterA = function (x.pos, y.pos, ht, wt) {
	x = c(0,  4,  6, 10,  8,  5, 2, 0, NA,2.2,2.6,7.4,7.8,2.2)
	y = c(0, 10, 10,  0,  0,7.5, 0, 0, NA,  3,  4,  4,  3,  3)
	x = 0.1 * x
	y = 0.1 * y
	x = x.pos + wt * x
	y = y.pos + ht * y
	fill = c("green4", "green4")
	list(x = x, y = y, fill = fill)
}

# draws the letter T at x=x.pos and y=y.pow with height=ht and width=wt
letterT = function (x.pos, y.pos, ht, wt) {
	x = c(0, 10, 10, 6, 6, 4, 4, 0)
	y = c(10, 10, 9, 9, 0, 0, 9, 9)
	x = 0.1 * x
	y = 0.1 * y
	x = x.pos + wt * x
	y = y.pos + ht * y
	fill = "red"
	list(x = x, y = y, fill = fill)
}

# draws the letter C at x=x.pos and y=y.pow with height=ht and width=wt
letterC = function (x.pos, y.pos, ht, wt) {
	angle1 = seq(0.3 + pi/2, pi, length = 100)
	angle2 = seq(pi, 1.5 * pi, length = 100)
	x.l1 = 0.5 + 0.5 * sin(angle1)
	y.l1 = 0.5 + 0.5 * cos(angle1)
	x.l2 = 0.5 + 0.5 * sin(angle2)
	y.l2 = 0.5 + 0.5 * cos(angle2)
	x.l = c(x.l1, x.l2)
	y.l = c(y.l1, y.l2)
	x = c(x.l, rev(x.l))
	y = c(y.l, 1 - rev(y.l))
	x.i1 = 0.5 + 0.35 * sin(angle1)
	y.i1 = 0.5 + 0.35 * cos(angle1)
	x.i1 = x.i1[y.i1 <= max(y.l1)]
	y.i1 = y.i1[y.i1 <= max(y.l1)]
	y.i1[1] = max(y.l1)
	x.i2 = 0.5 + 0.35 * sin(angle2)
	y.i2 = 0.5 + 0.35 * cos(angle2)
	x.i = c(x.i1, x.i2)
	y.i = c(y.i1, y.i2)
	x1 = c(x.i, rev(x.i))
	y1 = c(y.i, 1 - rev(y.i))
	x = c(x, rev(x1))
	y = c(y, rev(y1))
	x = x.pos + wt * x
	y = y.pos + ht * y
	fill = "blue"
	list(x = x, y = y, fill = fill)
}

# draws the letter G at x=x.pos and y=y.pow with height=ht and width=wt
letterG = function (x.pos, y.pos, ht, wt) {
	C = letterC(0, 0, 1, 1)
	x = C$x
	y = C$y
	r1 = max(x)
	h1 = 0.4
	x = c(x, NA, r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
	y = c(y, NA, h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
	x = x.pos + wt * x
	y = y.pos + ht * y
	fill = c("orange", "orange")
	list(x = x, y = y, fill = fill)
}

# calculates the information content of the given PWM
pwm2ic = function(pwm) {
	npos = ncol(pwm)
	ic = numeric(length = npos)
	for (i in 1:npos) {
		ic[i] = 2 + sum(sapply(pwm[, i], function(x) {
			if (x > 0) {
				x * log2(x)
			} else {
				0
			}
		}))
	}
	ic
}

# calculates the information content for each base in each column
myPwm2ic = function(pwm) {
	ic = matrix(0,nrow(pwm),ncol(pwm))
	for(k in 1:ncol(pwm)) {
		for(a in 1:nrow(pwm)) {
			ic[a,k] = pwm[a,k]*log2(pwm[a,k])
		}
	}
	return(ic);
}

# appends the letter which to the object letters
addLetter = function (letters, which, x.pos, y.pos, ht, wt) 
{
	if (which == "A") {
		letter = letterA(x.pos, y.pos, ht, wt)
	} else if (which == "C") {
		letter = letterC(x.pos, y.pos, ht, wt)
	} else if (which == "G") {
		letter = letterG(x.pos, y.pos, ht, wt)
	} else if (which == "T") {
		letter = letterT(x.pos, y.pos, ht, wt)
	} else {
		stop("which must be one of A,C,G,T")
	}
	letters$x = c(letters$x, NA, letter$x)
	letters$y = c(letters$y, NA, letter$y)
	letters$col = c(letters$col, letter$fill)
	letters
}

###
# Draws the classic sequence logo. 
#
# pwm:  the pwm for that the sequence logo should be drawn
# sparse: if TRUE margins are reduced and tickmarks are removed from the logo
# drawLines: a vector of y-values where to draw gray lines
#
seqLogo = function (pwm, sparse=FALSE, drawLines=c(0.5,1.0,1.5,2.0)) { 
    if (class(pwm) == "pwm") {
        pwm = pwm@pwm
    } else if (class(pwm) == "data.frame") {
        pwm = as.matrix(pwm)
    } else if (class(pwm) != "matrix") {
        print("pwm must be of class matrix or data.frame. Trying to convert")
	pwm = matrix(pwm,4,length(pwm)/4)
    }

    chars = c("A", "C", "G", "T")
    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm)

    ylim.negMax = 0;
    ylim.posMax = 0;
    facs = pwm2ic(pwm);
    
    ylab = "Information content [bits]"
   
    wt = 1.0
    x.pos = 0.5 # initial position on x axis is 0.5; Letter is one right from this point
    heights = c(); ymins=c(); ymaxs=c()
    for (j in 1:npos) {
	column = pwm[, j]
        hts = column * facs[j] 
	letterOrder = order(abs(hts)) # reorder letters
	yneg.pos = 0 
	ypos.pos = 0
	for (i in 1:4) {
	    letter = chars[letterOrder[i]]
	    ht = hts[letterOrder[i]]
   	    y.pos = ypos.pos;
	    ypos.pos = ypos.pos + ht + 0.0005
	    letters = addLetter(letters, letter, x.pos, y.pos, ht, wt)
	}
        x.pos = x.pos + wt
    }
  
    if(sparse) {
	plot(NA, xlim=c(0.5,x.pos), ylim=c(0,2),xaxt="n", ylab="",
	    mgp=c(0, .35, 0),tck=-0.02, cex.axis=0.8, frame.plot=F,xlab="")
    } else {
	plot(NA, xlim=c(0.5,x.pos), ylim=c(0,2), xaxt="n", ylab=ylab, frame.plot=F,xlab="Position")
    }

    for(y in drawLines) {
	abline(y,0,col="gray");
    }
    
    if(sparse) {
	axis(1,labels=c("",rep("",npos),""), at=c(0,1:npos,npos+1),tck=-0.02)
    } else {
         axis(1,labels=c("",1:npos,""),at=c(0,1:npos,npos+1))
    }
    polygon(letters, col=letters$col, border=F)
}

###
# Draws the difference of two sequence logos. 
#
# pwm1: the minuend pwm 
# pwm2: the subtrahend pwm
# type: a value between 1 and 6. Indicates the type of difference
# showSums: if TRUE the sum of a column is shown over the nucleotides. E.g. the overall difference of IC
# sparse: if TRUE margins are reduced and tickmarks are removed from the logo
#
diffSeqLogo = function (pwm1, pwm2, ymin=0, ymax=0, type=1, showSums=FALSE, sparse=FALSE) {
    if(type < 1 || type > 6) {
        stop("Unknown type. Type must be 1 <= type <= 6")
    }    

    if (class(pwm1) == "pwm") {
        pwm1 = pwm@pwm
    } else if (class(pwm1) == "data.frame") {
        pwm1 = as.matrix(pwm1)
    } else if (class(pwm1) != "matrix") {
	print("pwm1 must be of class matrix or data.frame. Trying to convert")
	pwm1 = matrix(pwm1,4,length(pwm1)/4)
    }

    if (class(pwm2) == "pwm") {
        pwm2 = pwm@pwm
    } else if (class(pwm2) == "data.frame") {
        pwm2 = as.matrix(pwm2)
    } else if (class(pwm2) != "matrix") {
        print("pwm2 must be of class matrix or data.frame. Trying to convert")
	pwm2 = matrix(pwm2,4,length(pwm2)/4)
    }


    if (any(abs(1 - apply(pwm1, 2, sum)) > 0.01)) { # check if the sum of each columns is near 1
        stop("Columns of PWM1 must add up to 1.0")
    }

    if (any(abs(1 - apply(pwm2, 2, sum)) > 0.01)) { # check if the sum of each columns is near 1
        stop("Columns of PWM2 must add up to 1.0")
    }

    if(ncol(pwm1) != ncol(pwm2)) {  # check if the two PWMs have the same length
        stop("The two given PWMs must have the same dimension");
    }

    # init needed variables
    diffPWM = pwm1-pwm2
    chars = c("A", "C", "G", "T")
    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(diffPWM)

    ylim.negMax = 0;
    ylim.posMax = 0;
    facs = numeric(npos) + 1; # heights are initialised with 1
    facs1 = pwm2ic(pwm1) # calc entropies for pwm1
    facs2 = pwm2ic(pwm2) # calc entropies for pwm2
	
    # set ylab
    if( type==1 ) {
	ylab = "Difference of probabilities"
    } else if( type==2 ) {
	ylab = "Difference of information content [bits]"
    } else if( type==3 ) {
	ylab = "Pairwise Difference of nucleotide information content [bits]"
    } else if( type==4 ) {
	ylab = "Relative difference of probabilities [%]"
    } else if( type==5 ) {
	ylab = "Relative difference of information content [%]"
    } else if( type==6 ) {
	ylab = "Relative pairwise difference of nucleotide information content [%]"
    } 

    # for type 1 and 5 recalculate heights
    if( type == 1 ) {
	facs = abs(facs1- facs2); # heights are set to the difference of ICs of the two given PWMs
	sums = colSums(abs(diffPWM))
	for( i in 1:npos) {
	    diffPWM[,i] = diffPWM[,i]/sums[i]	# normalize diffPWM
	}
    } else if (type==5) {
	facs = abs((facs1 - facs2)) / facs1 * 100; # heights are set to the ratio of ICs of the two given PWMs
	sums = colSums(abs(diffPWM))
	for( i in 1:npos) {
	    diffPWM[,i] = diffPWM[,i]/sums[i]	# normalize diffPWM
	}
    }
     
    wt = 1.0
    x.pos = 0.5 # initial position on x axis is 0.5; Letter is one right from this point
    heights = c(); ymins=c(); ymaxs=c()
    for (j in 1:npos) {
	column = diffPWM[, j]

	if(type==1 || type==2 || type==5) { # calculate heights of letters relative to overall hight
	    hts = -1.0* column * facs[j] 
	} else if( type==3 ) {
	    hts = -1.0* (pwm1[,j]*facs1[j] - pwm2[,j]*facs2[j])
	} else if( type==4 ) {
	    hts = -1.0* (pwm1[,j] - pwm2[,j]) / pwm1[,j] * 100
	} else if( type==6 ) {
	    hts = -1.0* (pwm1[,j]*facs1[j] - pwm2[,j]*facs2[j]) / pwm1[,j]*facs1[j] * 100
	}

	letterOrder = order(abs(hts)) # reorder letters
	yneg.pos = 0 
	ypos.pos = 0
	for (i in 1:4) {
	    letter = chars[letterOrder[i]]
	    ht = hts[letterOrder[i]]
	    if (ht >= 0){ 
	 	y.pos = ypos.pos;
		ypos.pos = ypos.pos + ht + 0.0005
	    } else if(ht < 0 ) {
		y.pos = yneg.pos;
		yneg.pos = yneg.pos + ht - 0.0005
	    }
	    letters = addLetter(letters, letter, x.pos, y.pos, ht, wt)
	}
	
        x.pos = x.pos + wt
	ylim.negMax = min(ylim.negMax, yneg.pos)
	ylim.posMax = max(ylim.posMax, ypos.pos)
	# remember values for plotting

	heights[j] = sum(abs(hts))
	ymins[j] = yneg.pos
	ymaxs[j] = ypos.pos
    }

    yAbsMax=2
    if( type==4 || type==5 || type==6) {
    	yAbsMax=400;
    }

    if(sparse) {
	plot(NA, xlim=c(0.5,x.pos), ylim=c(min(ylim.posMax*1.15,yAbsMax),max(ylim.negMax*1.15,-yAbsMax))*1.1, mgp=c(0, .35, 0),tck=-0.02, cex.axis=0.8,xaxt="n",yaxt="n", ylab="", frame.plot=F,xlab="")
    } else {
	plot(NA, xlim=c(0.5,x.pos), ylim=c(min(ylim.posMax*1.1,yAbsMax),max(ylim.negMax*1.1,-yAbsMax)), xaxt="n", ylab=ylab, frame.plot=F,xlab="Position")
    }

    # add sums over columns if wanted
    if( showSums ) {
	if(type==1 || type==2 || type==3) {
	   text(x=c(1:npos),y=ymins,paste(round(heights,2),"",sep=""),cex=0.8,font=3,pos=3)
	} else if(type==4 || type==5 || type==6) {
	   text(x=c(1:npos)+.20,y=ymins*1.01,paste(round(heights,0),"%",sep=""),cex=0.8,font=3,pos=3,srt=30)
	}
    }

    if( type==5 ) {
      yAt = (-5:5)*15
      yLabs = paste(yAt,"%",sep="")
    } else {
	yLabs = c("","0","");
	yAt = c(-yAbsMax,0,yAbsMax);
    }
	
    if(sparse) {
	axis(1,labels=c("",rep("",npos),""), at=c(0,1:npos,npos+1),tck=-0.02)
        axis(2,labels=yLabs,at=yAt,mgp=c(0, .35, 0),tck=-0.02, cex.axis=0.8)
    } else {
	axis(1,labels=c("",1:npos,""),at=c(0,1:npos,npos+1))
        axis(2,labels=c("",""),at=c(-yAbsMax,yAbsMax))
    }
    
    polygon(letters, col=letters$col, border=F)
    lines(c(0,x.pos), c(0,0) )
}



