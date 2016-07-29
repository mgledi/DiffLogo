##' Creates a DiffLogo object
##'
##' @title DiffLogo object
##' @param pwm1 representation of the first position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param pwm2 representation of the second position weight matrix (PWM) of type pwm, data.frame, or matrix
##' @param stackHeight function for the height of a stack at position i
##' @param baseDistribution function for the heights of the individual bases
##' @param alphabet of type Alphabet
##' @param align_pwms if True, will aling and extend pwms.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @return DiffLogo object
##' @export
##' @exportClass DiffLogo
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".txt",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = as.matrix(read.delim(file,header=FALSE))
##' }
##'
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##'
##' diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
##' diffLogo(diffLogoObj)
createDiffLogoObject = function (pwm1, pwm2, stackHeight=shannonDivergence,
                                 baseDistribution=
                                     normalizedDifferenceOfProbabilities,
                                 alphabet=DNA, align_pwms=F,
                                 unaligned_penalty=divergencePenaltyForUnaligned,
                                 try_reverse_complement=T) {
    pwm1 = preconditionTransformPWM(pwm1,alphabet);
    pwm2 = preconditionTransformPWM(pwm2,alphabet);
    preconditionPWM(pwm1);
    preconditionPWM(pwm2);
    if (align_pwms) {
        aligned_extended_pwm = alignExtendPwms(pwm1, pwm2,
                                               divergence=stackHeight,
                                               unaligned_penalty=
                                                   divergencePenaltyForUnaligned,
                                               try_reverse_complement=T)
        pwm1 = aligned_extended_pwm[[1]]
        pwm2 = aligned_extended_pwm[[2]]
    }
    preconditionPWMSameSize(pwm1,pwm2);

    # init needed variables
    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm1)
    eps = 0; # spacer between two bases in one stack
    ylim.negMax = 0;
    ylim.posMax = 0;

    npos = ncol(pwm1)
    wt = 1.0 # the width of one letter
    x.pos = 0.5 # initial position on x axis is 0.5; Letter is one right from this point
    heights = c(); ymins=c(); ymaxs=c()

    # determine intermediate values
    for (j in 1:npos) {
        heightObj = stackHeight(pwm1[,j],pwm2[,j]);
        preconditionStackHeight(heightObj); # check for correctness
        heights[j] = heightObj$height;
        ylab = heightObj$ylab;

        distr = baseDistribution(pwm1[,j],pwm2[,j]);
        preconditionBaseDistribution(distr); # check for correctness

        hts = distr*heights[j];
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
            letters = addLetter(letters, letterPolygons[[char]], x.pos, y.pos, ht, wt*0.99, col=col)
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
    diffObj$distance = sum(abs(heights)) # TODO as function
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
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".txt",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = as.matrix(read.delim(file,header=FALSE))
##' }
##'
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##'
##' diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
##' diffLogo(diffLogoObj)
diffLogo = function (diffLogoObj, ymin=0, ymax=0, sparse=FALSE) {
    if(class(diffLogoObj) != "DiffLogo") {
        msg = paste("Expected DiffLogo, but got ", class(diffLogoObj), ". Use #createDiffLogoObject to get an DiffLogo from two PWMs.",sep="")
        stop(msg)
    }

    if(ymin == 0) {
        ymin = diffLogoObj$ylim.posMax*1.0
    }
    if(ymax == 0) {
        ymax = diffLogoObj$ylim.negMax*1.0
    }

    # set ylab
    ylab = diffLogoObj$ylab

    if(sparse) {
        # the sparse plot has small ticks, small y-labels, no x-labels, no xlab, no ylab
        plot(NA, xlim=c(0.5,diffLogoObj$npos + 0.5), ylim=c(ymin,ymax), xaxt="n", ylab="", mgp=c(0, .35, 0), tck=-0.02, cex.axis=0.8, frame.plot=FALSE,xlab="")
    } else {
        plot(NA, xlim=c(0.5,diffLogoObj$npos + 0.5), ylim=c(ymin,ymax), xaxt="n", ylab=ylab, xlab="Position", frame.plot=FALSE, )
    }

    if(sparse) {
        axis(1,labels=c(rep("",diffLogoObj$npos)), at=(1:diffLogoObj$npos),tck=-0.02)
        axis(1,labels=c("",""), at=c(0,(diffLogoObj$npos+1)),tck=-0.00)
    } else {
        axis(1,labels=c(1:diffLogoObj$npos),at=(1:diffLogoObj$npos))
        axis(1,labels=c("",""), at=c(0,(diffLogoObj$npos+1)),tck=-0.00)
    }

    polygon(diffLogoObj$letters, col=diffLogoObj$letters$col, border=FALSE)
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
##' @param align_pwms if true, DiffLogo will align pwms before plotting
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".txt",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = as.matrix(read.delim(file,header=FALSE))
##' }
##'
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##'
##' diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2)
diffLogoFromPwm = function (
                      pwm1, pwm2, ymin=0, ymax=0,
                      stackHeight=shannonDivergence,
                      baseDistribution=normalizedDifferenceOfProbabilities,
                      sparse=FALSE, alphabet=DNA, align_pwms=F,
                      unaligned_penalty=divergencePenaltyForUnaligned,
                      try_reverse_complement=T) {
    diffLogoObj = createDiffLogoObject(
                      pwm1,pwm2,stackHeight=stackHeight,
                      baseDistribution=baseDistribution,
                      alphabet=alphabet,
                      align_pwms=align_pwms,
                      unaligned_penalty=unaligned_penalty,
                      try_reverse_complement=try_reverse_complement);
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
##' @param enableClustering if TRUE the motifs are reordered, so that similar motifs have a small vertical and horizontal distance in the table
##' @param margin the space reseverved for labels
##' @param ratio the ratio of the plot; this is needed to determine the margin sizes correctly
##' @param alphabet of type Alphabet
##' @param align_pwms if True, will aling and extend pwms.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param ... set of parameters passed to the function 'axis' for plotting
##' @export
##' @importFrom cba order.optimal
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".txt",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = as.matrix(read.delim(file,header=FALSE))
##' }
##'
##' diffLogoTable(motifs)
diffLogoTable = function (
            PWMs,
            names2,
            stackHeight=shannonDivergence,
            baseDistribution=normalizedDifferenceOfProbabilities,
            uniformYaxis=TRUE,
            sparse=TRUE,
            showSequenceLogosTop=TRUE,
            enableClustering=TRUE,
            treeHeight=0.5,
            margin=0.02,
            ratio=1,
            alphabet=DNA,
            align_pwms=F,
            unaligned_penalty=divergencePenaltyForUnaligned,
            try_reverse_complement=T,
            ...
) {
    plot.new();
    dim = length(PWMs);

    st = 0;
    if ( showSequenceLogosTop ) {
        st = 0.5;
    }
    if( !enableClustering ) {
        treeHeight=0;
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
            similarities[i,k] = 0
            if( i != k ) {
                diffLogoObj = createDiffLogoObject(PWMs[[ motif_i ]], PWMs[[ motif_k ]],
                                                   stackHeight=stackHeight,
                                                   baseDistribution=baseDistribution,
                                                   alphabet=alphabet, align_pwms=align_pwms,
                                                   unaligned_penalty=unaligned_penalty,
                                                   try_reverse_complement=try_reverse_complement);
                if(uniformYaxis) {
                    ymin = min(diffLogoObj$ylim.negMax,ymin)
                    ymax = max(diffLogoObj$ylim.posMax,ymax)
                }
                similarities[i,k] = diffLogoObj$distance;
            }
        }
    }
    colors = matrix(palette[cut(similarities,100)],dim,dim)
    leafOrder=1:dim;
    if(enableClustering) {
        distance = dist(similarities);
        hc = hclust(distance, "average");
        opt = order.optimal(distance,hc$merge)
        hc$merge=opt$merge
        hc$order=opt$order
        leafOrder = hc$order;
    }

    # draw DiffLogos
    dimV = c(dim, dim,dim+st+treeHeight,dim+st+treeHeight);
    for ( i in 1:dim) {
        motif_i = names[leafOrder[i]];
        for ( k in 1:dim) {
            if( i != k ) {
                motif_k = names[leafOrder[k]];
                subplotcoords = c(i-1,i,dim-k,dim-k+1)

                par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,0,0), new=TRUE, mar=c(0,0,0,0))

                plot(NA,ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
                rect(0,0,1,1,col=colors[leafOrder[i],leafOrder[k]],border=NA);

                par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,0,0), new=TRUE, mar=marDiffLogo)
                diffLogoObj = createDiffLogoObject(PWMs[[ motif_i ]],PWMs[[ motif_k ]],
                                                   stackHeight=stackHeight,
                                                   baseDistribution=baseDistribution,
                                                   alphabet=alphabet,
                                                   align_pwms=align_pwms,
                                                   unaligned_penalty=unaligned_penalty,
                                                   try_reverse_complement=try_reverse_complement);
                diffLogo(diffLogoObj,sparse=sparse,ymin=ymin,ymax=ymax)
            }
        }
        if(showSequenceLogosTop) {
            subplotcoords = c(i-1, i, dim, dim + st)
            #par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=c(0,0,0,0))
            #plot(NA,ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
            #rect(0,0,1,1,col="gray",border=NA);
            par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=marSeqLogo)
            seqLogo(PWMs[[ motif_i ]],sparse=sparse, alphabet=alphabet)
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
    axis(2, pos=0, at= (1:dim) - 0.5, labels = rev(names2[leafOrder]), tick = FALSE, mgp = c(3, 0, 0), ...)
    axis(3, pos=dim, at= (1:dim) - 0.5, labels = names2[leafOrder], tick = FALSE, mgp = c(3, 0, 0), ...)
}

##' Counts PWM divergence as sum of divergencies of their columns.
##'
##' @title PWM divergence
##' @param pwm_left is a PWM representation in type of matrix
##' @param pwm_right is a PWM representation in type of matrix. The result is symmetric on pwm_left and pwm_right
##' @param divergence is a Divergence function on columns.
##' @return float - sum of divergences
##' @export
pwmDivergence = function(pwm_left, pwm_right, divergence=shannonDivergence) {
    stopifnot(ncol(pwm_left) == ncol(pwm_right))
    return(sum(sapply(1:ncol(pwm_left), function(i){
        divergence(pwm_left[,i], pwm_right[,i])$height
    })))
}

uniformPwm = function(pwm_length, alphabet_length) {
    return(matrix(rep(rep(1.0/alphabet_length, each=alphabet_length),
                      each=length, pwm_length), nrow=alphabet_length))
}

divergencePenaltyForUnaligned = function(pwm, unaligned_at_left,
                                         aligned_length,
                                         unaligned_at_right, divergence) {
    shift_divergence = 0
    if (unaligned_at_left > 0) {
       shift_divergence = shift_divergence+
                          pwmDivergence(matrix(
                              pwm[,1:unaligned_at_left], nrow=nrow(pwm)),
                              uniformPwm(unaligned_at_left,
                                         nrow(pwm)),
                              divergence)
    }
    if (unaligned_at_right > 0) {
       aligned_length = ncol(pwm) - unaligned_at_right
       shift_divergence = shift_divergence+
                          pwmDivergence(
                              matrix(pwm[,(aligned_length+1):ncol(pwm)],
                                     nrow=nrow(pwm)),
                              uniformPwm(ncol(pwm)-aligned_length,
                                  nrow(pwm)), divergence)
    }
    return(shift_divergence)
}

findBestShiftForPwms = function(static_pwm, shifted_pwm, divergence,
                                unaligned_penalty) {
    best_divergence = Inf
    best_shift = 0
    for (shift in 0:(ncol(static_pwm)-1)){
        intersection_length = min(ncol(static_pwm) - shift, ncol(shifted_pwm))
        shift_divergence = pwmDivergence(
                               matrix(
                                   static_pwm[,
                                       (shift+1):(shift+intersection_length)],
                                                nrow=nrow(static_pwm)),
                                    matrix(
                                        shifted_pwm[,1:intersection_length],
                                        nrow=nrow(shifted_pwm)),
                                    divergence=divergence)
        shift_divergence = shift_divergence+unaligned_penalty(
                                      static_pwm,
                                      unaligned_at_left = shift,
                                      intersection_length,
                                      unaligned_at_right =
                                         max(0, ncol(static_pwm)
                                                   -(shift+intersection_length)),
                                      divergence)
        shift_divergence = shift_divergence+unaligned_penalty(
                                      shifted_pwm,
                                      unaligned_at_left = 0,
                                      intersection_length,
                                      unaligned_at_right =
                                         max(0, ncol(shifted_pwm)
                                                   -intersection_length),
                                      divergence)
        shift_divergence = shift_divergence #/ max(ncol(static_pwm),
                                            #      shift + ncol(shifted_pwm))
        if (shift_divergence < best_divergence) {
           best_divergence = shift_divergence
           best_shift = shift
        }
    }
    result = list()
    result$divergence = best_divergence
    result$shift = best_shift
    return(result)
}

revCompPwm = function (pwm) {
    result = pwm[nrow(pwm):1, ncol(pwm):1]
    rownames(result) = rownames(pwm)
    return(result)
}

##' Finds best local alignment for PWMs.
##'
##' @title Align pwms
##' @param PWM is a matrix of type matrix
##' @param
##' @param divergence is a measure of difference between two pwm columns. Smaller is more similar. If you want to use non-uniform background distribution, provide your own function.
##' @param distance for unaligned columns at edges of matrixes. See divergencePenaltyForUnaligned as an example for providing your own function
##' @param If false the alignment will not be performed on reverse complements. If true, the input pwms should have column order of ACTG.
##' @export
##' @author Lando Andrey
localPwmAlignment = function(pwm_left, pwm_right, divergence=shannonDivergence,
                             unaligned_penalty=divergencePenaltyForUnaligned,
                             try_reverse_complement=T) {
    no_change = list("shift"=0, "direction"="forward")
    result = list()
    best_divergence = Inf

    alignment = findBestShiftForPwms(pwm_left, pwm_right, divergence=divergence,
                                     unaligned_penalty=unaligned_penalty)
    if (alignment$divergence < best_divergence) {
       result[[1]] = no_change
       result[[2]] = list("shift" = alignment$shift, "direction"="forward")
       result$divergence = alignment$divergence
       best_divergence = alignment$divergence
    }

    alignment = findBestShiftForPwms(pwm_right, pwm_left, divergence=divergence,
                             unaligned_penalty=unaligned_penalty)
    if (alignment$divergence < best_divergence) {
       result[[1]] = list("shift" = alignment$shift, "direction"="forward")
       result[[2]] = no_change
       result$divergence = alignment$divergence
       best_divergence = alignment$divergence
    }
    if (try_reverse_complement) {
        alignment = findBestShiftForPwms(pwm_left, revCompPwm(pwm_right),
                                         divergence=divergence,
                                         unaligned_penalty=unaligned_penalty)
        if (alignment$divergence < best_divergence) {
           result[[1]] = no_change
           result[[2]] = list("shift" = alignment$shift, "direction"="reverse")
           result$divergence = alignment$divergence
           best_divergence = alignment$divergence
        }

        alignment = findBestShiftForPwms(revCompPwm(pwm_right), pwm_left,
                                         divergence=divergence,
                                         unaligned_penalty=unaligned_penalty)
        if (alignment$divergence < best_divergence) {
           result[[1]] = list("shift" = alignment$shift, "direction"="forward")
           result[[2]] = list("shift" = 0, "direction"="reverse")
           result$divergence = alignment$divergence
           best_divergence = alignment$divergence
        }
    }
    return(result)
}

addBefore = function (matrix, to_add_length) {
    return(cbind(uniformPwm(to_add_length, nrow(matrix)), matrix))
}

addAfter = function (matrix, to_add_length) {
    return(cbind(matrix, uniformPwm(to_add_length, nrow(matrix))))
}

##' Finds best local alignment for PWMs and extends them from both sides to equal length
##'
##' @title Align pwms and extend them
##' @param PWM is a matrix of type matrix
##' @param
##' @param divergence is a measure of difference between two pwm columns. Smaller is more similar. If you want to use non-uniform background distribution, provide your own function.
##' @param distance for unaligned columns at edges of matrixes. See divergencePenaltyForUnaligned as an example for providing your own function
##' @param If false the alignment will not be performed on reverse complements. If true, the input pwms should have column order of ACTG.
##' @export
##' @author Lando Andrey
alignExtendPwms = function(left_pwm, right_pwm, divergence=shannonDivergence,
                           unaligned_penalty=divergencePenaltyForUnaligned,
                           try_reverse_complement=T) {
    alignment = localPwmAlignment(left_pwm, right_pwm, divergence=divergence,
                           unaligned_penalty=unaligned_penalty,
                           try_reverse_complement=try_reverse_complement)
    if (alignment[[1]]$direction == 'reverse') {
        left_pwm = revCompPwm(left_pwm)
    }
    left_pwm = addBefore(left_pwm, alignment[[1]]$shift)
    if (alignment[[2]]$direction == 'reverse') {
        right_pwm = revCompPwm(right_pwm)
    }
    right_pwm = addBefore(right_pwm, alignment[[2]]$shift)
    if (ncol(right_pwm) > ncol(left_pwm)) {
       left_pwm = addAfter(left_pwm, ncol(right_pwm)-ncol(left_pwm))
    }
    if (ncol(right_pwm) < ncol(left_pwm)) {
       right_pwm = addAfter(right_pwm, ncol(left_pwm)-ncol(right_pwm))
    }
    return(list(left_pwm, right_pwm))
}
