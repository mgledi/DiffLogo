##' Default configuration list for diffLogoTable
##'
##' @title Configuration object for diffLogoTable
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
##' @param align_pwms if True, will align and extend pwms in each cell of diffLogoTable independently.
##' @param multiple_align_pwms if True, will align and extend pwms in the diffLogoTable jointly.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param length_normalization if True, divergence between pwms is divided by length of pwms.
##' @param numberOfPermutations number of permutations for the permutation test for the calculation of p-values
##' @param ... set of parameters passed to the function 'axis' for plotting
##' @return list of parameters
##' @export
##' @author Lando Andrey
##' @examples
##' diffLogoTableConfiguration(DNA)
diffLogoTableConfiguration = function(
		alphabet,
		stackHeight=shannonDivergence,
		baseDistribution=normalizedDifferenceOfProbabilities,
		uniformYaxis=TRUE,
		sparse=TRUE,
		showSequenceLogosTop=TRUE,
		enableClustering=TRUE,
		treeHeight=0.5,
		margin=0.02,
		ratio=1,
		align_pwms=FALSE,
		multiple_align_pwms=TRUE,
		unaligned_penalty=divergencePenaltyForUnaligned,
		try_reverse_complement=TRUE,
		length_normalization=FALSE,
		numberOfPermutations=100) {
    #if (!alphabet$supportReverseComplement) {
    #   try_reverse_complement = F;
    #}
    return(
      list(
    		uniformYaxis=uniformYaxis,
    		sparse=sparse,
    		showSequenceLogosTop=showSequenceLogosTop,
    		enableClustering=enableClustering,
    		treeHeight=treeHeight,
    		margin=margin,
    		ratio=ratio,
    		stackHeight=stackHeight,
    		baseDistribution=baseDistribution,
    		multiple_align_pwms=multiple_align_pwms,
    		align_pwms=align_pwms,
    		unaligned_penalty=unaligned_penalty,
    		try_reverse_complement=try_reverse_complement && alphabet$supportReverseComplement,
    		length_normalization=length_normalization,
    		numberOfPermutations=numberOfPermutations)
    )
}

getMarginsDiffLogo = function(sparse) {
	if(sparse) {
		return(c(0.3,1.2,0.1,0.1));
	} else {
		return(c(1,1.5,0.1,0.1));
	}
}

getMarginsSeqLogo = function(sparse) {
	if(sparse) {
		return (c(0.3,1.2,0.0,0.1));
	} else {
		return (c(1,1.5,0.1,0.1));
	}
}

extractNames = function(PWMs) {
	names = names(PWMs);
    if (is.null(names)) {
        names = 1:length(PWMs)
    }
    return (names);
}

revCompPwm = function (pwm) {
    if(nrow(pwm) != 4)
      stop("Can only reverse complement DNA and RNA alphabet.")
    result = pwm[nrow(pwm):1, ncol(pwm):1]
    rownames(result) = rownames(pwm)
    return(result)
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
    stopifnot(ncol(pwm_left) == ncol(pwm_right));
    return(sum(sapply(1:ncol(pwm_left), function(i) {
        divergence(pwm_left[,i], pwm_right[,i])$height
    })));
}

##' Generates a PWM consisting of only the uniform distribution or the given base_distribution (if defined).
##' 
##' @title Generates a PWM
##' @param pwm_length the number of positions
##' @param alphabet_length the alphabet size
##' @param base_distribution optional base distribution for each PWM position
##' @return a PWM
baseDistributionPwm = function(pwm_length, alphabet_length, base_distribution=NULL) {
    if (is.null(base_distribution)) {
       base_distribution = rep(1.0/alphabet_length, each=alphabet_length)
    }
    return(matrix(rep(base_distribution, each=length, pwm_length), nrow=alphabet_length))
}

##' Enriches a matrix of difflogo objects with p-values which quantifies the probability that two PWM-positions are from the same distribution
##'
##' @title Enriches a matrix of difflogo objects with p-values
##' @param diffLogoObjMatrix matrix of difflogo objects
##' @param sampleSizes number of sequences behind the pwms behind the given difflogo objects
##' @param stackHeight function for the calculation of a divergence measure for two probability vectors
##' @param numberOfPermutations the number of permutations to perform for the calculation of stackHeights
##' @return matrix of difflogo objects enriched with p-values
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
##' names(sampleSizes) <- motif_names
##' 
##' diffLogoTableObj = prepareDiffLogoTable(motifs);
##' diffLogoTableObj$diffLogoObjMatrix = enrichDiffLogoTableWithPvalues(diffLogoTableObj$diffLogoObjMatrix, sampleSizes)
enrichDiffLogoTableWithPvalues <- function(diffLogoObjMatrix, sampleSizes, stackHeight=shannonDivergence, numberOfPermutations = 100 ) {
  if(is.null(sampleSizes))
    stop("No sample sizes given!")
  
  motifs = names(diffLogoObjMatrix);
  dim = length(motifs);
  for ( i in 1:dim) {
    for ( k in 1:dim) {
      motif_i = motifs[i];
      motif_k = motifs[k];
      if(!is.null(diffLogoObjMatrix[[motif_i]][[motif_k]])) {
        diffLogoObjMatrix[[motif_i]][[motif_k]] = enrichDiffLogoObjectWithPvalues(
                                                          diffLogoObjMatrix[[motif_i]][[motif_k]],
                                                          sampleSizes[[motif_i]],
                                                          sampleSizes[[motif_k]]
                                                    );
      }
    }
  }
  return(diffLogoObjMatrix);
}

##' Enriches a difflogo object with p-values which quantifies the probability that two PWM-positions are from the same distribution
##'
##' @title Enriches a difflogo object with p-values
##' @param diffLogoObj matrix of difflogo objects
##' @param n1 the number of sequences behind the first pwm behind the given difflogo object
##' @param n2 the number of sequences behind the second pwm behind the given difflogo object
##' @param stackHeight function for the calculation of a divergence measure for two probability vectors
##' @param numberOfPermutations the number of permutations to perform for the calculation of stackHeights
##' @return enriched difflogo object
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
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##' n1 <- 100
##' n2 <- 100
##' diffLogoObj = createDiffLogoObject(pwm1 = pwm1, pwm2 = pwm2)
##' diffLogoObj = enrichDiffLogoObjectWithPvalues(diffLogoObj, n1, n2)
enrichDiffLogoObjectWithPvalues <- function(diffLogoObj, n1, n2, stackHeight=shannonDivergence, numberOfPermutations = 100) {
    pwm1 = diffLogoObj$pwm1
    pwm2 = diffLogoObj$pwm2
    npos = ncol(diffLogoObj$pwm1);
    pvals = rep(1,npos);
    for (j in (diffLogoObj$unaligned_from_left+1):(npos - diffLogoObj$unaligned_from_right)) {
      pvals[j] = calculatePvalue(p1 = pwm1[,j], p2 = pwm2[,j], n1 = n1, n2 = n2, stackHeight = stackHeight, numberOfPermutations = numberOfPermutations);
    }
    diffLogoObj$pvals = pvals;
    return(diffLogoObj);
} 


##' Calculates the p-value for the null-hypothesis that two given probability vectors p1, p2 calculated from n1/n2 observations arise from the same distribution
##'
##' @title p-value that two PWM-positions are from the same distribution
##' @param p1 first probability vector with one probability for each symbol of the alphabet
##' @param p2 second probability vector with one probability for each symbol of the alphabet
##' @param n1 number of observations for the calculation of p1
##' @param n2 number of observations for the calculation of p2
##' @param stackHeight function for the calculation of a divergence measure for two probability vectors
##' @param numberOfPermutations the number of permutations to perform for the calculation of stackHeights
##' @param plotGammaDistributionFit if TRUE the fit of a gamma distribution to the sampled stackHeights is plotted
##' @return a numeric p-value
##' @export
##' @author Hendrik Treutler
##' @examples
##' p1 <- c(0.2, 0.3, 0.1, 0.4)
##' p2 <- c(0.2, 0.1, 0.3, 0.4)
##' n1 <- 100
##' n2 <- 200
##' numberOfPermutations = 100
##' plotGammaDistributionFit = TRUE
##' 
##' pValue <- calculatePvalue(p1 = p1, p2 = p2, n1 = n1, n2 = n2, stackHeight = shannonDivergence, numberOfPermutations = numberOfPermutations, plotGammaDistributionFit = plotGammaDistributionFit)
calculatePvalue <- function(p1, p2, n1, n2, stackHeight=shannonDivergence, numberOfPermutations = 100, plotGammaDistributionFit = FALSE){
  ############################################################################################
  ## preconditions
  preconditionVectorSameSize(p1, p2)
  preconditionBaseDistribution(p1)
  preconditionBaseDistribution(p2)
  preconditionProbabilityVector(p1)
  preconditionProbabilityVector(p2)

  if(all(p1==p2)) {
    return(1);
  }

  if(any(is.na(n1), is.na(n2), is.null(n1), is.null(n2), !is.numeric(n1), !is.numeric(n2), n1<=0, n2<=0))
    stop("Given counts are corrupt!")
  
  ############################################################################################
  ## divergence to test
  height = stackHeight(p1 = p1, p2 = p2)
  preconditionStackHeight(height)
  observedDivergence = height$height
  
  ############################################################################################
  ## parameters
  alphabet <- 1:length(p1)
  p = (p1 + p2) / 2
  n = n1 + n2
  a = as.integer(round(p * n))
  # the rounding can cause that sum(a) != n
  n = sum(a);

  ############################################################################################
  ## permutations
  multiplier <- 1:length(alphabet)
  seed <- (sum((1 / (p1+1e-3)) * multiplier) * sum((1 / (p2+1e-3)) * multiplier)) %% 2147483647
  set.seed(seed)
  
  symbols <- unlist(sapply(X = alphabet, FUN = function(x){rep(x = x, times = a[[x]])}))
  classes <- c(rep(x = TRUE, times = n1), rep(x = FALSE, times = n - n1))

  divergences <- vector(mode = "numeric", length = numberOfPermutations)
  for(idx in 1:numberOfPermutations){
    newOrder = classes[sample(n)];
    symbols2 = symbols[newOrder]

    a1 = unlist(lapply(X = alphabet, FUN = function(x){sum(symbols2 == x)}))
    a2 = a - a1
    divergences[[idx]] = stackHeight(p1 = a1 / n1, p2 = a2 / n2)$height
  }
  maximumDivergence = max(max(divergences) * 2, observedDivergence * 2)
  
  ############################################################################################
  ## fit gamma distribution
  med.gam <- mean(divergences)
  var.gam <- var(divergences)
  rate <- med.gam/var.gam
  alpha <- ((med.gam)^2)/var.gam
  scale <- 1/rate
  
  gammaDistX <- seq(from=0, to=maximumDivergence, length.out = 10000)

  gammaDistY <- dgamma(x = gammaDistX, rate = rate, shape = alpha)
  if(gammaDistY[1] == Inf) {
    gammaDistY[1] = 0;
  }
  gammaDistY = gammaDistY / sum(gammaDistY)  
  gammaDistYcum = cumsum(gammaDistY)
  
  ############################################################################################
  ## plot
  if(plotGammaDistributionFit){
    cumSumValuesX <- sort(divergences)
    cumSumValuesY <- (1:numberOfPermutations) / numberOfPermutations
    
    plot(NA, xlim = c(0,maximumDivergence), ylim = c(0, 1), ylab = "p", xlab = "x", main = paste("p = (", paste(p, collapse = ","), "), n1 = ", n1, ", n2 = ", n2, sep = ""))
    lines(x = cumSumValuesX, y = cumSumValuesY, col = "blue")
    lines(x = gammaDistX, y = gammaDistYcum, col = "green")
    legend(x = max(divergences)/2, y = 0.3, legend = c("Empiric", paste("Gamma, k=", round(rate, digits = 2), ", theta=", round(alpha, digits = 2), sep = "")), lty = c(1, 1), col = c("blue", "green"))
  }
  
  ############################################################################################
  ## p-value
  pValue <- 1 - gammaDistYcum[[min(which(gammaDistX > observedDivergence))]]
  
  return(pValue)
}