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
seqLogo = function (pwm, sparse=FALSE, drawLines=c(0.5,1.0,1.5,2.0)) { 
    pwm = preconditionTransformPWM(pwm);
    preconditionPWM(pwm);

    chars = c("A", "C", "G", "T")
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