# calculates the information content of the given PWM
pwm2ic_col = function(pwm) {
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
pwm2ic_base = function(pwm) {
	ic = matrix(0,nrow(pwm),ncol(pwm))
	for(k in 1:ncol(pwm)) {
		for(a in 1:nrow(pwm)) {
			ic[a,k] = pwm[a,k]*log2(pwm[a,k])
		}
	}
	return(ic);
}


###
# Draws the classic sequence logo. 
#
# pwm:  the pwm for that the sequence logo should be drawn
# sparse: if TRUE margins are reduced and tickmarks are removed from the logo
# drawLines: a vector of y-values where to draw gray lines
#
seqLogo = function (pwm, sparse=FALSE, drawLines=c(0.5,1.0,1.5,2.0),alphabet=DNA) { 
    pwm = preconditionTransformPWM(pwm);
    preconditionPWM(pwm);

    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm)

    ylim.negMax = 0;
    ylim.posMax = 0;
    facs = pwm2ic_col(pwm);
    
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
	    ht = hts[letterOrder[i]]
   	    y.pos = ypos.pos;
	    ypos.pos = ypos.pos + ht + 0.0005
	    char = DNA$chars[letterOrder[i]]
	    letters = addLetter(letters, DNA$letters[[char]], x.pos, y.pos, ht, wt)
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
    polygon(letters, col=letters$col, border=letters$col)
}


# appends the letter which to the object letters
addLetter = function (letters, letter, x.pos, y.pos, ht, wt) 
{
	x = x.pos + wt * letter$x
	y = y.pos + ht * letter$y
	polygons = sum(is.na(x))+1
	letters$x = c(letters$x, NA, x)
	letters$y = c(letters$y, NA, y)
	letters$col = c(letters$col, rep(letter$col,polygons))
	letters
}
