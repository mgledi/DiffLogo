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
##' @param length_normalization If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
##' @param unaligned_from_left the number of unaligned positions on the left
##' @param unaligned_from_right the number of unaligned positions on the right
##' @return DiffLogo object
##' @export
##' @exportClass DiffLogo
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##'
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##'
##' diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
##' diffLogo(diffLogoObj)
createDiffLogoObject = function (pwm1, pwm2, 
                                 stackHeight=shannonDivergence,
                                 baseDistribution= normalizedDifferenceOfProbabilities,
                                 alphabet=DNA, align_pwms=FALSE,
                                 unaligned_penalty=divergencePenaltyForUnaligned,
                                 try_reverse_complement=TRUE,
                                 base_distribution=NULL,
                                 length_normalization = FALSE,
                                 unaligned_from_left = 0,
                                 unaligned_from_right = 0) {
    pwm1 = preconditionTransformPWM(pwm1,alphabet);
    pwm2 = preconditionTransformPWM(pwm2,alphabet);
    preconditionPWM(pwm1);
    preconditionPWM(pwm2);
    if (align_pwms) {
        alignment = localPwmAlignment(
            pwm1, pwm2,
            divergence = stackHeight,
            unaligned_penalty = unaligned_penalty,
            try_reverse_complement = try_reverse_complement && alphabet$supportReverseComplement,
            base_distribution = base_distribution,
            length_normalization = length_normalization);
        aligned_extended_pwms = extendPwmsFromAlignmentVector(list(pwm1, pwm2), alignment$vector, base_distribution);        
        unaligned_from_left = max(alignment$vector[[1]]$shift, alignment$vector[[2]]$shift);
        alignment_length = max(alignment$vector[[1]]$shift+ncol(pwm1), alignment$vector[[2]]$shift+ncol(pwm2));
        unaligned_from_right = max( alignment_length-ncol(pwm1)-alignment$vector[[1]]$shift, alignment_length-ncol(pwm2)-alignment$vector[[2]]$shift);

        pwm1 = aligned_extended_pwms[[1]]
        pwm2 = aligned_extended_pwms[[2]]
    }
    preconditionPWMSameSize(pwm1,pwm2);


    # init needed variables
    letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos = ncol(pwm1)
    ylim.negMax = 0;
    ylim.posMax = 0;

    npos = ncol(pwm1)
    wt = 1.0 # the width of one letter
    x.pos = 0.5 + unaligned_from_left # initial position on x axis is 0.5; Letter is one right from this point
    heights = rep(0,npos); ymins=rep(0,npos); ymaxs=rep(0,npos); 

    for (j in (unaligned_from_left+1):(npos - unaligned_from_right)) {
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
            ht_total = sum(hts);
            if (ht >= 0){
                y.pos = ypos.pos;
                ypos.pos = ypos.pos + ht
            } else if(ht < 0 ) {
                y.pos = yneg.pos;
                yneg.pos = yneg.pos + ht
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
    diffObj$unaligned_from_left = unaligned_from_left
    diffObj$unaligned_from_right = unaligned_from_right
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
##' @param diffLogoConfiguration list of configuration parameters (see function diffLogoTableConfiguration(...))
##' @return none (draws difference logo)
##' @export
##' @importFrom graphics plot axis polygon rect text lines
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##'
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##'
##' diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
##' diffLogo(diffLogoObj)
diffLogo = function (diffLogoObj, ymin=0, ymax=0, sparse=FALSE, diffLogoConfiguration = list()) {
    if(!is(diffLogoObj, "DiffLogo")) {
        msg = paste("Expected DiffLogo, but got ", class(diffLogoObj), ". Use #createDiffLogoObject to get an DiffLogo from two PWMs.",sep="")
        stop(msg)
    }

    if(ymin == 0) {
        ymin = diffLogoObj$ylim.posMax
    }
    if(ymax == 0) {
        ymax = diffLogoObj$ylim.negMax
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
    if (!is.null(diffLogoObj$unaligned_from_left) && diffLogoObj$unaligned_from_left>0) {
        rect(0.5, -ymin, diffLogoObj$unaligned_from_left+0.5, -ymax, col="gray", border="gray")
    }

    if (!is.null(diffLogoObj$unaligned_from_right) && diffLogoObj$unaligned_from_right>0) {
        rect(diffLogoObj$npos-diffLogoObj$unaligned_from_right+0.5, -ymin, diffLogoObj$npos+0.5, -ymax, col="gray", border="gray")
    }

    if(!is.null(diffLogoObj$pvals)) {
        leftOffset = 0;
        if(!is.null(diffLogoObj$unaligned_from_left)) {
            leftOffset = diffLogoObj$unaligned_from_left;
        }
        if(!is.null(diffLogoObj$unaligned_from_right)) {
            rightOffset = diffLogoObj$unaligned_from_right;
        }
        npos = ncol(diffLogoObj$pwm1);
        for (j in (leftOffset+1):(npos - rightOffset)) {
            if( diffLogoObj$pvals[j] < 0.05  ) {
                text(j, ymin, "*");
            }
        }
    }
    
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
##' @param length_normalization If true, will minimize the average divergence between PWMs. Otherwise will minimize the sum of divergences between positions. In both cases unalignes positions are compared to base_distribution and are counted when computing the alignment length.
##' @return none (draws difference logo)
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
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
                      sparse=FALSE, alphabet=DNA, align_pwms=FALSE,
                      unaligned_penalty=divergencePenaltyForUnaligned,
                      try_reverse_complement=TRUE, base_distribution=NULL,
                      length_normalization = FALSE) {
    diffLogoObj = createDiffLogoObject(
                      pwm1, pwm2, stackHeight=stackHeight,
                      baseDistribution=baseDistribution,
                      alphabet=alphabet,
                      align_pwms=align_pwms,
                      unaligned_penalty=unaligned_penalty,
                      try_reverse_complement=try_reverse_complement,
                      base_distribution=NULL, length_normalization = length_normalization);
    diffLogo(diffLogoObj, ymin=ymin, ymax=ymax, sparse=sparse)
}

##' Prepares a DiffLogoTable and generates an object that contains the hirarchical clustering and a matrix of prepared difference logos. 
##'
##' @title Prepare a table of difflogo objects
##' @param PWMs a list/vector of position weight matrices (PWMs) each of type pwm, data.frame, or matrix
##' @param alphabet the alphabet of the given PWMs
##' @param configuration list of (probably part of) of configuration options. See diffLogoTableConfiguration.
##' @return matrix of difference logos
##' @importFrom utils modifyList
##' @importFrom stats dist
##' @importFrom stats hclust
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##' sampleSizes <- c(100, 150, 200, 250)
##' 
##' diffLogoTableObj = prepareDiffLogoTable(motifs);
prepareDiffLogoTable = function (
            PWMs,
            alphabet=DNA,
            configuration=list()
) {
    configuration        = modifyList(diffLogoTableConfiguration(alphabet), configuration)
    enableClustering     = configuration$enableClustering
    stackHeight            = configuration$stackHeight
    baseDistribution       = configuration$baseDistribution
    multiple_align_pwms    = configuration$multiple_align_pwms
    align_pwms             = configuration$align_pwms
    unaligned_penalty      = configuration$unaligned_penalty;
    try_reverse_complement = (configuration$try_reverse_complement && alphabet$supportReverseComplement)
    length_normalization   = configuration$length_normalization
    numberOfPermutations   = configuration$numberOfPermutations
    dim = length(PWMs);
    similarities = matrix(0,dim,dim);
    names = extractNames(PWMs);
    diffLogoTableObj = list();
    for ( i in 1:dim) {
        PWMs[[names[i]]] = preconditionTransformPWM(PWMs[[names[i]]],alphabet);
    }

    if (multiple_align_pwms) {
        multiple_pwms_alignment = multipleLocalPwmsAlignment(PWMs, try_reverse_complement = try_reverse_complement);
        originalPwmLengths = lapply(PWMs, ncol);
        PWMs = extendPwmsFromAlignmentVector( PWMs, multiple_pwms_alignment$alignment$vector);
        stopifnot(dim==length(PWMs))
        align_pwms = FALSE
    }
    
    diffLogoObjMatrix = list();
    for ( i in 1:dim) {
        motif_i = names[i];
        diffLogoObjMatrix[[ motif_i ]] = list()

        # alignment size in columns is always equal
        alignment_length = ncol(PWMs[[motif_i]]);

        for ( k in 1:dim) {
            if( i != k ) {
                motif_k = names[k];
            
                unaligned_from_left = 0
                unaligned_from_right = 0
                if (multiple_align_pwms) {
                    # the number of unaligned left positions is the maximum of the two shifts of the motifs
                    unaligned_from_left = max( multiple_pwms_alignment$alignment$vector[[i]]$shift,
                                               multiple_pwms_alignment$alignment$vector[[k]]$shift)
                    
                    # the number of unaligned right positions is the length of the alignment minus the orignal pwm length
                    rightShiftMotif_i = alignment_length - originalPwmLengths[[ motif_i ]] - multiple_pwms_alignment$alignment$vector[[i]]$shift;
                    rightShiftMotif_k = alignment_length - originalPwmLengths[[ motif_k ]] - multiple_pwms_alignment$alignment$vector[[k]]$shift;
                    unaligned_from_right = max(rightShiftMotif_i, rightShiftMotif_k)
                }
                diffLogoObjMatrix[[motif_i]][[motif_k]] = 
                        createDiffLogoObject(  PWMs[[ motif_i ]],
                                               PWMs[[ motif_k ]],
                                               stackHeight = stackHeight,
                                               baseDistribution = baseDistribution,
                                               alphabet = alphabet,
                                               align_pwms = align_pwms,
                                               unaligned_penalty = unaligned_penalty,
                                               try_reverse_complement = try_reverse_complement,
                                               length_normalization = length_normalization,
                                               unaligned_from_left = unaligned_from_left,
                                               unaligned_from_right = unaligned_from_right);
                similarities[i, k] = diffLogoObjMatrix[[motif_i]][[motif_k]]$distance;
                similarities[k, i] = diffLogoObjMatrix[[motif_i]][[motif_k]]$distance;
            } 
        }
    }

    leafOrder=1:dim;
    if(enableClustering) {
        distance = dist(similarities);
        hc = hclust(distance, "average");
        opt = order.optimal(distance,hc$merge);
        hc$merge=opt$merge;
        hc$order=opt$order;

        diffLogoTableObj$hc = hc;
        leafOrder = diffLogoTableObj$hc$order;
    }
    
    diffLogoTableObj$diffLogoObjMatrix = diffLogoObjMatrix;
    diffLogoTableObj$PWMs = PWMs;
    diffLogoTableObj$leafOrder = leafOrder;
    diffLogoTableObj$configuration = configuration;
    diffLogoTableObj$alphabet = alphabet;
    return (diffLogoTableObj);
}

##' Draws a table of DiffLogos.
##'
##' @title Draws a table of DiffLogos
##' @param diffLogoTableObj the diffLogoTable-Object created by function prepareDiffLogoTable(...)
##' @param ... optional parameters for functon axis
##' @return none (draws difference logo)
##' @importFrom grDevices colorRampPalette rgb
##' @importFrom graphics plot.new par plot rect axis
##' @export
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##'
##' diffLogoTableObj = prepareDiffLogoTable(motifs)
##' drawDiffLogoTable(diffLogoTableObj)
drawDiffLogoTable = function (
            diffLogoTableObj,
            ...
) {
    alphabet             = diffLogoTableObj$alphabet
    configuration        = diffLogoTableObj$configuration
    enableClustering     = configuration$enableClustering
    uniformYaxis         = configuration$uniformYaxis
    sparse               = configuration$sparse
    showSequenceLogosTop = configuration$showSequenceLogosTop
    treeHeight           = configuration$treeHeight
    margin               = configuration$margin
    ratio                = configuration$ratio

    PWMs = diffLogoTableObj$PWMs;
    diffLogoObjMatrix = diffLogoTableObj[['diffLogoObjMatrix']];

    marDiffLogo = getMarginsSeqLogo(sparse);
    marSeqLogo = getMarginsSeqLogo(sparse);
    ymin = 0; ymax = 0;
    if( !enableClustering ) {
        treeHeight=0;
    }
    st = 0;
    if(showSequenceLogosTop) {
        st = 0.5;
    }
    orderedMotifs = names(diffLogoObjMatrix)[diffLogoTableObj$leafOrder];
    
    dim = length(orderedMotifs);
    similarities = matrix(0,dim,dim);

    # Filling similarity matrix and computing y axes limits.
    # The order of motifs in diffLogoObjMatrix is the same as in orderedMotifs
    for ( i in 1:dim) {
        for ( k in 1:dim) {
            motif_i = orderedMotifs[i];
            motif_k = orderedMotifs[k];
            if(!is.null(diffLogoObjMatrix[[motif_i]][[motif_k]])) { 
                similarities[i,k] = diffLogoObjMatrix[[motif_i]][[motif_k]]$distance;
            } else {
                similarities[i,k]  = 0;
            }
            if(uniformYaxis) {
                ymin = min(diffLogoObjMatrix[[motif_i]][[motif_k]]$ylim.negMax, ymin)
                ymax = max(diffLogoObjMatrix[[motif_i]][[motif_k]]$ylim.posMax, ymax)
            }
        }
    }
    palette = colorRampPalette(c(rgb(0.9,1,0.9),rgb(1,0.9,0.9)))(100)
    colors = matrix(palette[cut(similarities,100)],dim,dim)

    # drawing the matrix of DiffLogos
    plot.new();
    dimV = c(dim, dim, dim + st + treeHeight, dim + st + treeHeight);
    for ( i in 1:dim) {
        for ( k in 1:dim) {
            motif_i = orderedMotifs[i];
            motif_k = orderedMotifs[k];
            if(!is.null(diffLogoObjMatrix[[motif_i]][[motif_k]])) { 
                subplotcoords = c(i-1,i,dim-k,dim-k+1)

                par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,0,0), new=TRUE, mar=c(0,0,0,0))

                plot(NA,ylim=c(0,1),xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n")
                rect(0,0,1,1,col=colors[i,k],border=NA);

                par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,0,0), new=TRUE, mar=marDiffLogo)

                diffLogo(diffLogoObjMatrix[[motif_i]][[motif_k]], sparse=sparse, ymin=ymin, ymax=ymax);
            }
        }
    }
                
    # draw sequence logos
    if(showSequenceLogosTop) {
        for(i in 1:dim) {
            subplotcoords = c(i-1, i, dim, dim + st)
            par(fig=(subplotcoords / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=marSeqLogo)
            seqLogo(PWMs[[  orderedMotifs[i] ]],sparse=sparse, alphabet=alphabet)
        }
    }

    # draw cluster tree
    if(treeHeight > 0 && enableClustering) {
        par(fig=c(0 + 0.5, dim - 0.5, dim + st, dim + st + treeHeight) / dimV * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio) + c(margin,margin,margin*ratio,margin*ratio), new=TRUE, mar=c(.0,.0,.1,.0))
        plot(diffLogoTableObj$hc, xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n", labels=rep("",dim), main="" )
    }

    # add names
    par(fig=(c(0,dim,0,dim) / dimV) * c(1-margin,1-margin,1-margin*ratio,1-margin*ratio)+ c(margin,margin,0,0), new=TRUE, mar=c(0,0,0,0))

    plot(NA,ylim=c(0,dim),xlim=c(0,dim),xaxs="i",xaxt="n",yaxt="n",yaxs="i", bty="n", mar=c(0,0,0,0))
    axis(2, pos=0, at= (1:dim) - 0.5, labels = rev(orderedMotifs), tick = FALSE, mgp = c(3, 0, 0), ...)
    axis(3, pos=dim, at= (1:dim) - 0.5, labels = orderedMotifs, tick = FALSE, mgp = c(3, 0, 0), ...)
}

##' Draws a table of DiffLogos.
##'
##' @title Draw DiffLogo-table
##' @param PWMs a list/vector of position weight matrices (PWMs) each of type pwm, data.frame, or matrix
##' @param sampleSizes the number of sequences behind each PWM
##' @param alphabet the alphabet of the given PWMs
##' @param configuration list of (probably part of) of configuration options. See diffLogoTableConfiguration.
##' @param ... set of parameters passed to the function 'axis' for plotting
##' @return none (draws table of difference logos)
##' @export
##' @importFrom cba order.optimal
##' @importFrom utils modifyList
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##'
##' diffLogoTable(motifs)
diffLogoTable = function (
            PWMs,
            sampleSizes=NULL,
            alphabet=DNA,
            configuration=list(),
            ...
) {
    if(sum(names(configuration) %in% names(diffLogoTableConfiguration(alphabet))) != length(configuration)) {
        stop(paste("Unknown arguments passed to diffLogoTable:",  paste(names(configuration[ !(names(configuration) %in% names(diffLogoTableConfiguration(alphabet)))]), sep=",")))
    }
    if(is.null(names(PWMs))) {
      names(PWMs) = sapply((1:length(PWMs)), toString)
    }
    
    configuration        = modifyList(diffLogoTableConfiguration(alphabet), configuration)
    diffLogoTableObj = prepareDiffLogoTable(PWMs, alphabet, configuration, ...);
    
    if(!is.null(sampleSizes)) {
      diffLogoTableObj$diffLogoObjMatrix = enrichDiffLogoTableWithPvalues(diffLogoTableObj$diffLogoObjMatrix, sampleSizes, configuration$stackHeight, configuration$numberOfPermutations)
    }
    drawDiffLogoTable(diffLogoTableObj, ... );
}
