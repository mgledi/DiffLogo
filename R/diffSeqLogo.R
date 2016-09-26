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
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @param If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
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
                                 try_reverse_complement=T,
                                 base_distribution=NULL,
                                 length_normalization = F) {
    pwm1 = preconditionTransformPWM(pwm1,alphabet);
    pwm2 = preconditionTransformPWM(pwm2,alphabet);
    preconditionPWM(pwm1);
    preconditionPWM(pwm2);
    if (align_pwms) {
        alignment = localPwmAlignment(
            pwm1, pwm2,
            divergence=stackHeight,
            unaligned_penalty=
                unaligned_penalty,
            try_reverse_complement=T,
            base_distribution=base_distribution,
            length_normalization = length_normalization)
        aligned_extended_pwms = extendPwmsFromAlignmentVector(list(pwm1, pwm2),
                                                              alignment$vector,
                                                              base_distribution)
        pwm1 = aligned_extended_pwms[[1]]
        pwm2 = aligned_extended_pwms[[2]]
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
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @param If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
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
                      try_reverse_complement=T, base_distribution=NULL,
                      length_normalization = F) {
    diffLogoObj = createDiffLogoObject(
                      pwm1,pwm2,stackHeight=stackHeight,
                      baseDistribution=baseDistribution,
                      alphabet=alphabet,
                      align_pwms=align_pwms,
                      unaligned_penalty=unaligned_penalty,
                      try_reverse_complement=try_reverse_complement,
                      base_distribution=NULL, length_normalization = length_normalization);
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
            multiple_align_pwms=F,
            align_pwms=F,
            unaligned_penalty=divergencePenaltyForUnaligned,
            try_reverse_complement=T,
            length_normalization = F,
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
    if (is.null(names)) {
        names = 1:length(PWMs)
    }

    ymin = 0;
    ymax = 0;
    if (multiple_align_pwms) {
        multiple_pwms_alignment = multipleLocalPwmsAlignment(PWMs)
        PWMs = extendPwmsFromAlignmentVector(
                   PWMs,
                   multiple_pwms_alignment$alignment$vector)
        stopifnot(dim==length(PWMs))
    }
    for ( i in 1:dim) {
        motif_i = names[i];
        for ( k in 1:dim) {
            motif_k = names[k];
            similarities[i,k] = 0
            if( i != k ) {
                diffLogoObj = createDiffLogoObject(PWMs[[ motif_i ]],
                                                   PWMs[[ motif_k ]],
                                                   stackHeight=stackHeight,
                                                   baseDistribution=
                                                       baseDistribution,
                                                   alphabet=alphabet,
                                                   align_pwms=align_pwms,
                                                   unaligned_penalty=
                                                       unaligned_penalty,
                                                   try_reverse_complement=
                                                       try_reverse_complement,
                                                   length_normalization =
                                                       length_normalization);

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
        if (multiple_align_pwms) {
            hc = list()
            distance_matrix = multiple_pwms_alignment$distance_matrix
            for (i in 1:ncol(distance_matrix)) {
                distance_matrix[[i,i]] = 0
                for (j in 1:i) {
                    distance_matrix[[i, j]] = distance_matrix[[j, i]]
                }
            }
            opt = order.optimal(as.dist(distance_matrix),
                                multiple_pwms_alignment$merge)
            hc$merge = opt$merge
            leafOrder = opt = hc$order = opt$order

            hc$merge = multiple_pwms_alignment$merge
            hc$height = multiple_pwms_alignment$height
            hc$labels = names(PWMs)
            class(hc) = 'hclust'
        } else {
            distance = dist(similarities);
            hc = hclust(distance, "average");
            opt = order.optimal(distance,hc$merge)
            hc$merge=opt$merge
            hc$order=opt$order
            leafOrder = hc$order;
        }
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
                diffLogoObj = createDiffLogoObject(PWMs[[ motif_i ]],
                                                   PWMs[[ motif_k ]],
                                                   stackHeight=stackHeight,
                                                   baseDistribution =
                                                       baseDistribution,
                                                   alphabet=alphabet,
                                                   align_pwms=align_pwms,
                                                   unaligned_penalty=
                                                       unaligned_penalty,
                                                   try_reverse_complement=
                                                       try_reverse_complement,
                                                   length_normalization=
                                                       length_normalization);
                diffLogo(diffLogoObj, sparse=sparse, ymin=ymin, ymax=ymax)
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

baseDistributionPwm = function(pwm_length, alphabet_length, base_distribution=NULL) {
    if (is.null(base_distribution)) {
       base_distribution = rep(1.0/alphabet_length, each=alphabet_length)
    }
    return(matrix(rep(base_distribution,
                      each=length, pwm_length), nrow=alphabet_length))
}

divergencePenaltyForUnaligned = function(pwm, unaligned_at_left,
                                         aligned_length,
                                         unaligned_at_right, divergence,
                                         base_distribution) {
    shift_divergence = 0
    if (unaligned_at_left > 0) {
       shift_divergence = shift_divergence+
                          pwmDivergence(matrix(
                              pwm[,1:unaligned_at_left], nrow=nrow(pwm)),
                              baseDistributionPwm(unaligned_at_left,
                                                  nrow(pwm),
                                                  base_distribution=
                                                      base_distribution),
                              divergence)
    }
    if (unaligned_at_right > 0) {
       aligned_length = ncol(pwm) - unaligned_at_right
       shift_divergence = shift_divergence+
                          pwmDivergence(
                              matrix(pwm[,(aligned_length+1):ncol(pwm)],
                                     nrow=nrow(pwm)),
                              baseDistributionPwm(ncol(pwm)-aligned_length,
                                                  nrow(pwm),
                                                  base_distribution=
                                                      base_distribution),
                              divergence)
    }
    return(shift_divergence)
}

findBestShiftForPwms = function(static_pwm, shifted_pwm, divergence,
                                unaligned_penalty, base_distribution,
                                length_normalization) {
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
                                      divergence,
                                      base_distribution)
        shift_divergence = shift_divergence+unaligned_penalty(
                                      shifted_pwm,
                                      unaligned_at_left = 0,
                                      intersection_length,
                                      unaligned_at_right =
                                         max(0, ncol(shifted_pwm)
                                                   -intersection_length),
                                      divergence,
                                      base_distribution)
        if (length_normalization) {
           shift_divergence = shift_divergence / max(ncol(static_pwm),
                                                  shift + ncol(shifted_pwm))
        }
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


##' Finds best local alignment for PWMs.
##'
##' @title Align pwms
##' @param PWM is a matrix of type matrix
##' @param
##' @param divergence is a measure of difference between two pwm columns. Smaller is more similar. If you want to use non-uniform background distribution, provide your own function.
##' @param distance for unaligned columns at edges of matrixes. See divergencePenaltyForUnaligned as an example for providing your own function
##' @param If false the alignment will not be performed on reverse complements. If true, the input pwms should have column order of ACTG.
##' @param If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
##' @export
##' @author Lando Andrey
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
##' @param If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
##' @export
##' @author Lando Andrey
localPwmAlignment = function(pwm_left, pwm_right, divergence=shannonDivergence,
                             unaligned_penalty=divergencePenaltyForUnaligned,
                             try_reverse_complement=T, base_distribution=NULL,
                             length_normalization = F) {
    no_change = list("shift"=0, "direction"="forward")
    result_alignment_vector = list()
    best_divergence = Inf

    alignment = findBestShiftForPwms(pwm_left, pwm_right, divergence=divergence,
                                     unaligned_penalty=unaligned_penalty,
                                     base_distribution,
                                     length_normalization=length_normalization)
    if (alignment$divergence < best_divergence) {
       result_alignment_vector[[1]] = no_change
       result_alignment_vector[[2]] = list("shift" = alignment$shift,
                                           "direction"="forward")
       result_divergence = alignment$divergence
       best_divergence = alignment$divergence
    }

    alignment = findBestShiftForPwms(pwm_right, pwm_left, divergence=divergence,
                                     unaligned_penalty=unaligned_penalty,
                                     base_distribution=base_distribution,
                                     length_normalization=length_normalization)
    if (alignment$divergence < best_divergence) {
        result_alignment_vector[[1]] = list("shift" = alignment$shift,
                                            "direction"="forward")
        result_alignment_vector[[2]] = no_change
        result_divergence = alignment$divergence
        best_divergence = alignment$divergence
    }
    if (try_reverse_complement) {
        alignment = findBestShiftForPwms(pwm_left, revCompPwm(pwm_right),
                                         divergence=divergence,
                                         unaligned_penalty=unaligned_penalty,
                                         base_distribution=base_distribution,
                                         length_normalization=length_normalization)
        if (alignment$divergence < best_divergence) {
            result_alignment_vector[[1]] = no_change
            result_alignment_vector[[2]] = list("shift" = alignment$shift,
                                                "direction"="reverse")
            result_divergence = alignment$divergence
            best_divergence = alignment$divergence
        }

        alignment = findBestShiftForPwms(revCompPwm(pwm_right), pwm_left,
                                         divergence=divergence,
                                         unaligned_penalty=unaligned_penalty,
                                         base_distribution=base_distribution,
                                         length_normalization=length_normalization)
        if (alignment$divergence < best_divergence) {
            result_alignment_vector[[1]] = list("shift" = alignment$shift,
                                                "direction"="forward")
            result_alignment_vector[[2]] = list("shift" = 0,
                                                "direction"="reverse")
            result_divergence = alignment$divergence
            best_divergence = alignment$divergence
        }
    }
    return(list("vector"=result_alignment_vector,
                "divergence"=result_divergence))
}

addBefore = function (matrix, to_add_length, base_distribution=NULL) {
    return(cbind(baseDistributionPwm(to_add_length,
                                     nrow(matrix),
                                     base_distribution), matrix))
}

addAfter = function (matrix, to_add_length, base_distribution=NULL) {
    return(cbind(matrix,
        baseDistributionPwm(to_add_length, nrow(matrix),
                            base_distribution)))
}

##' Extends pwms by adding base_distribution to both sides, so that they keep aligned, but have equal length.
##'
##' @title Extend pwms with respect to alignment
##' @param pwms is a list of matrixes
##' @param alignment_vector is a list of shifts ($shift) and orientations ($direction)
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @export
##' @author Lando Andrey
extendPwmsFromAlignmentVector = function(pwms, alignment_vector,
                                         base_distribution=NULL) {
    stopifnot(length(pwms)==length(alignment_vector))
    max_length = 0
    for (i in 1:length(pwms)) {
        if (alignment_vector[[i]]$direction == 'reverse') {
            pwms[[i]] = revCompPwm(pwms[[i]])
        }
        pwms[[i]] = addBefore(pwms[[i]], alignment_vector[[i]]$shift,
                              base_distribution)
        max_length = max(max_length, ncol(pwms[[i]]))
    }
    for (i in 1:length(pwms)) {
        pwms[[i]] = addAfter(pwms[[i]], max_length-ncol(pwms[[i]]),
                             base_distribution)
    }
    return(pwms)
}

##' Computes average pwm divergence from alignment vector for two datasets. This equals to average divergence between all pairs where one pwm comes from left set, and other comes from right
##'
##' @title Average divergence between two sets.
##' @param left_pwms_list is a list of matrixes
##' @param left_pwms_alignment is a list of shifts ($shift) and orientations ($direction)
##' @param right_pwms_list is a list of matrixes
##' @param right_pwms_alignment is a list of shifts ($shift) and orientations ($direction)
##' @param divergence divergence measure.
##' @export
##' @return float
##' @author Lando Andrey
twoSetsAveragePwmDivergenceFromAlignmentVector = function(left_pwms_list,
                                                          left_pwms_alignment,
                                                          right_pwms_list,
                                                          right_pwms_alignment,
                                                          divergence=
                                                              shannonDivergence){
    stopifnot(length(left_pwms_list)==length(left_pwms_alignment))
    stopifnot(length(right_pwms_list)==length(right_pwms_alignment))
    pwms_divergences = c()
    for (left_pwm_index in 1:length(left_pwms_list)) {
        for (right_pwm_index in 1:length(right_pwms_list)) {
            stopifnot(!is.null(left_pwms_list[left_pwm_index]))
            stopifnot(!is.null(right_pwms_list[right_pwm_index]))

            stopifnot(!is.null(left_pwms_alignment[left_pwm_index]))
            stopifnot(!is.null(right_pwms_alignment[right_pwm_index]))
            pwms = extendPwmsFromAlignmentVector(c(left_pwms_list[left_pwm_index],
                                                   right_pwms_list[right_pwm_index]),
                                                 c(left_pwms_alignment[left_pwm_index],
                                                   right_pwms_alignment[right_pwm_index]))
            pwms_divergences = c(pwms_divergences,
                                 pwmDivergence(pwms[[1]], pwms[[2]],
                                               divergence=divergence))
        }
    }
    return(mean(pwms_divergences))
}

findBestShiftForPwmSets = function(static_pwms_list, static_pwms_alignment,
                                   shifted_pwms_list, shifted_pwms_alignment,
                                   divergence=shannonDivergence,
                                   unaligned_penalty=
                                       divergencePenaltyForUnaligned) {
    stopifnot(length(static_pwms_list)==length(static_pwms_alignment$vector))
    stopifnot(length(shifted_pwms_list)==length(shifted_pwms_alignment$vector))
    best_divergence = Inf
    best_shift = 0
    max_possible_shift = 0
    result = list()
    for (i in 1:length(static_pwms_list)) {
        max_possible_shift = max(max_possible_shift,
                                 ncol(static_pwms_list[[i]])
                                     +static_pwms_alignment$vector[[i]]$shift)
    }
    for (shift in 0:max_possible_shift) {
        divergence_for_shift = twoSetsAveragePwmDivergenceFromAlignmentVector(
            static_pwms_list, static_pwms_alignment$vector,
            shifted_pwms_list, shifted_pwms_alignment$vector,
            divergence=divergence)
        if (divergence_for_shift < best_divergence) {
            best_divergence = divergence_for_shift
            result = list("vector" = c(static_pwms_alignment$vector,
                                       shifted_pwms_alignment$vector),
                          "divergence"=best_divergence)
        }
        for (i in 1:length(shifted_pwms_alignment$vector)) {
            shifted_pwms_alignment$vector[[i]]$shift =
                shifted_pwms_alignment$vector[[i]]$shift+1
        }
    }
    return(result)
}

##' If 'forward' return 'reverse'
##' If 'reverse' return 'forward'
switchDirection = function(direction) {
    if (direction=='forward') {
       return('reverse')
    }
    if (direction=='reverse') {
       return('forward')
    }
    stopifnot(F);
}

##' Returns alignment vector as if all pwm were reverted.
##'
##' @title Reverse for alignment vector
##' @param alignment_vector list of list which $shift and $orientation
##' @param pwms list of matrixes.
##' @export
##' @return list - reversed alignment vector
##' @author Lando Andrey
reverseAlignmentVector = function(alignment_vector, pwms) {
    stopifnot(length(alignment_vector)==length(pwms))
    alignment_length = 0
    for (i in 1:length(alignment_vector)) {
        alignment_vector[[i]]$direction = switchDirection(
            alignment_vector[[i]]$direction)
        stopifnot(alignment_vector[[i]]$shift >= 0)
        alignment_length = max(alignment_length,
                               ncol(pwms[[i]]) + alignment_vector[[i]]$shift)
    }
    min_shift = 0
    for (i in 1:length(alignment_vector)) {
        alignment_vector[[i]]$shift =
            alignment_length - ncol(pwms[[i]]) - alignment_vector[[i]]$shift
        min_shift = min(min_shift, alignment_vector[[i]]$shift)
    }
    for (i in 1:length(alignment_vector)) {
        alignment_vector[[i]]$shift = -min_shift + alignment_vector[[i]]$shift
    }
    return(alignment_vector)
}


##' Align two sets of pwms
##'
##' @title Multiple PWMs alignment
##' @param left_pwms_set list of pwms(matrixes)
##' @param left_aligment alignment of left_pwms_set. 
##' @param right_pwms_set list of pwms;
##' @param right_alignment alignment of right_pwms_set.
##' @param try_reverse_complement if true(default), also try reverse complement.
##' @return list - alignment of concatination of left_pwms_set and right_pwms_set
##' @export
##' @exportClass DiffLogo
##' @author Lando Andrey
alignPwmSets = function(left_pwms_set, left_alignment,
                        right_pwms_set, right_alignment) {
    stopifnot(length(left_pwms_set)==length(left_alignment$vector))
    stopifnot(length(right_pwms_set)==length(right_alignment$vector))
    result = list()
    best_divergence = Inf
    alignment = findBestShiftForPwmSets(left_pwms_set, left_alignment,
                                        right_pwms_set, right_alignment)
    if (alignment$divergence < best_divergence) {
       result = alignment
       best_divergence = alignment$divergence
    }

    alignment = findBestShiftForPwmSets(right_pwms_set, right_alignment,
                                        left_pwms_set, left_alignment)
    if (alignment$divergence < best_divergence) {
       best_divergence = alignment$divergence
       result$vector = c(alignment$vector[(length(right_pwms_set)+1):(length(alignment$vector))],
                            alignment$vector[1:length(right_pwms_set)])
       result$divergence = best_divergence
    }
    right_alignment$vector = reverseAlignmentVector(right_alignment$vector,
                                                       right_pwms_set)
    alignment = findBestShiftForPwmSets(left_pwms_set, left_alignment,
                                        right_pwms_set, right_alignment)
    if (alignment$divergence < best_divergence) {
       result = alignment
       best_divergence = alignment$divergence
    }

    alignment = findBestShiftForPwmSets(right_pwms_set, right_alignment, 
                                        left_pwms_set, left_alignment)
    if (alignment$divergence < best_divergence) {
       best_divergence = alignment$divergence
       result$vector = c(alignment$vector[(length(right_pwms_set)+1):(length(alignment$vector))],
                            alignment$vector[1:length(right_pwms_set)])
       result$divergence = best_divergence
    }
    return(result)
}

##' Creates a distance matrix for pwms
##'
##' @title Multiple PWMs alignment
##' @param pwms list of pwms
##' @param diagonal_value value to put on diagonal.
##' @param bottom_default_value value to put on bottom triangle. Set to NULL to get symmetric distance matrix.
##' @param divergence divergence measure.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @param length_normalization is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @return list
##' @export
##' @exportClass DiffLogo
##' @author Lando Andrey
pwmsDistanceMatrix = function(pwms,
                              diagonal_value = 0,
                              bottom_default_value = NULL,
                              divergence=shannonDivergence,
                              unaligned_penalty=divergencePenaltyForUnaligned,
                              try_reverse_complement=T, base_distribution=NULL,
                              length_normalization = F){
    distance_matrix = matrix(0, ncol=length(pwms), nrow=length(pwms))
    for (i in 1:(length(pwms)-1)) {
        distance_matrix[[i, i]] = diagonal_value
        for (j in (i+1):length(pwms)) {
            distance_matrix[[i, j]] = localPwmAlignment(
                pwms[[i]], pwms[[j]],
                divergence=divergence,
                unaligned_penalty=unaligned_penalty,
                try_reverse_complement=try_reverse_complement,
                base_distribution=base_distribution,
                length_normalization = length_normalization)$divergence
            if (is.null(bottom_default_value)) {
                distance_matrix[[j, i]] = distance_matrix[[i, j]]
            } else {
                distance_matrix[[j, i]] = bottom_default_value
            }
        }
    }
    distance_matrix[[length(pwms), length(pwms)]] = diagonal_value
    return(distance_matrix)
}

minimumMatrixElementIndexes = function(matrix) {
    minimum_element = which.min(matrix)
    col_of_minimum = round(minimum_element/ncol(matrix)-0.5)+1
    row_of_minimum = minimum_element - ((col_of_minimum-1)*(
        ncol(matrix)))
    first_index = min(col_of_minimum, row_of_minimum)
    second_index = max(col_of_minimum, row_of_minimum)
    return(list(first_index, second_index))
}

addLastTreeNodeToDistanceMatrix = function(distance_matrix, tree_nodes) {
    distance_matrix = cbind(distance_matrix, rep(Inf,
                                                 ncol(distance_matrix)))
    distance_matrix = rbind(distance_matrix, rep(Inf,
                                                 ncol(distance_matrix)))
    last = ncol(distance_matrix)
    if (last > 0) {
        for (i in 1:(last-1)) {
            distance_matrix[[i, last]] = alignPwmSets(
                tree_nodes[[last]]$pwms, tree_nodes[[last]]$pwms_alignment,
                tree_nodes[[i]]$pwms,    tree_nodes[[i]]$pwms_alignment)$divergence
        }
    }
    return(distance_matrix)
}

joinTwoNodesInAlignmentTree = function(distance_matrix, tree_nodes, join_counter) {
    min_pwm_row_col = minimumMatrixElementIndexes(distance_matrix)
    first_pwm_index = min_pwm_row_col[[1]]
    second_pwm_index = min_pwm_row_col[[2]]
    pwms_alignment = alignPwmSets(tree_nodes[[first_pwm_index]]$pwms,
                                  tree_nodes[[first_pwm_index]]$pwms_alignment,
                                  tree_nodes[[second_pwm_index]]$pwms,
                                  tree_nodes[[second_pwm_index]]$pwms_alignment)
    stopifnot(length(c(tree_nodes[[first_pwm_index]]$pwms,
                       tree_nodes[[second_pwm_index]]$pwms)) == length(pwms_alignment$vector))
    tree_nodes = c(tree_nodes,
                   list(
                       list("pwms"=c(tree_nodes[[first_pwm_index]]$pwms,
                                     tree_nodes[[second_pwm_index]]$pwms),
                            "left_child"=tree_nodes[[first_pwm_index]],
                            "right_child"=tree_nodes[[second_pwm_index]],
                            "num_leafs"=tree_nodes[[first_pwm_index]]$num_leaf+
                                tree_nodes[[second_pwm_index]]$num_leaf,
                            "node_number"=join_counter,
                            "pwms_alignment" = pwms_alignment
                          )))
    join_counter = join_counter+1
    tree_nodes[[first_pwm_index]] = NULL
    tree_nodes[[second_pwm_index-1]] = NULL
    distance_matrix = matrix(
        distance_matrix[-c(first_pwm_index, second_pwm_index),
                        -c(first_pwm_index, second_pwm_index)],
        nrow = nrow(distance_matrix)-2, ncol = ncol(distance_matrix)-2)
    if (ncol(distance_matrix) != 0) {
        distance_matrix = addLastTreeNodeToDistanceMatrix(distance_matrix, tree_nodes)
    }
    return(list("distance_matrix"=distance_matrix,
                "tree_nodes"=tree_nodes))
}

alignmentTree = function(pwms, distance_matrix) {
    tree_nodes = list()
    for (i in 1:length(pwms)) {
        tree_nodes = c(tree_nodes,
                       list(
                           list("pwms"=list(pwms[[i]]),
                                "left_child"=NULL,
                                "right_child"=NULL,
                                "num_leafs"=1,
                                "pwm_index"=i,
                                "pwms_alignment"=list(
                                    "vector"=list(list("direction"="forward",
                                                          "shift"=0))))))
    }
    join_counter = 1
    for (step in 1:(length(pwms)-1)) {
        joined = joinTwoNodesInAlignmentTree(distance_matrix, tree_nodes, join_counter)
        distance_matrix = joined$distance_matrix
        tree_nodes = joined$tree_nodes
        join_counter = join_counter+1
    }
    return(tree_nodes[[1]])
}

alignmentTreeLeftRightTriversal = function(node) {
    merge = list()
    order = list()
    height = list()
    leftRightTreeTriversalHelper = function(node){
        if (is.null(node$left_child)) {
            order <<- c(order, node$pwm_index)
            return(-(node$pwm_index))
        } else {
            left_index = leftRightTreeTriversalHelper(node$left_child)
            right_index = leftRightTreeTriversalHelper(node$right_child)
            merge <<- c(merge, list(c(left_index, right_index)))
            height <<- c(height, list(node$pwms_alignment$divergence))
            #height <<- c(height, list(node$node_number))
        }
        return(length(merge))
    }
    leftRightTreeTriversalHelper(node)
    return(list("merge"=merge, "order"=unlist(order), "height"=height))
}


##' Creates a multiple alignment of pwms
##'
##' @title Multiple PWMs alignment
##' @param pwms list of pwms
##' @param divergence Divergence measure.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param stackHeight function for the height of a stack at position i
##' @param base_distribution is a vector of length nrow(pwm) that is added to unaligned columns of pwms for comparing. If NULL, uniform distribution is used
##' @return list
##' @export
##' @exportClass DiffLogo
##' @author Lando Andrey
multipleLocalPwmsAlignment = function(
  pwms,
  divergence=shannonDivergence,
  unaligned_penalty=divergencePenaltyForUnaligned,
  try_reverse_complement=T,
  base_distribution=NULL,
  length_normalization = F)
{
    distance_matrix = pwmsDistanceMatrix(
        pwms,
        bottom_default_value = Inf,
        diagonal_value = Inf,
        divergence=divergence,
        unaligned_penalty=unaligned_penalty,
        try_reverse_complement=try_reverse_complement,
        base_distribution=base_distribution,
        length_normalization = length_normalization)
    alignment_tree_nodes = alignmentTree(pwms, distance_matrix)
    traversal_result = alignmentTreeLeftRightTriversal(alignment_tree_nodes)
    alignment_tree_nodes$pwms_alignment$vector =
        alignment_tree_nodes$pwms_alignment$vector[
                                                sort.list(
                                                    unlist(
                                                        traversal_result$order))]
    return(list("alignment"=alignment_tree_nodes$pwms_alignment,
                "order"=unlist(traversal_result$order),
                "merge"=t(matrix(unlist(traversal_result$merge), 2, byrow=F)),
                "height"=unlist(traversal_result$height),
                "raw_tree"=list(alignment_tree_nodes),
                "distance_matrix"=distance_matrix
                ))
}
