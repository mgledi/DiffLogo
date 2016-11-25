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
##' @param align_pwms if True, will align and extend pwms in each sell of diffLogoTable independently.
##' @param unaligned_penalty is a function for localPwmAlignment.
##' @param try_reverse_complement if True, alignment will try reverse complement pwms
##' @param length_normalization if True, divergence between pwms is divided by length of pwms.
##' @param ... set of parameters passed to the function 'axis' for plotting
##' @export
##' @author Lando Andrey
##' @examples
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
		align_pwms=F,
		multiple_align_pwms=T,
		unaligned_penalty=divergencePenaltyForUnaligned,
		try_reverse_complement=T,
		length_normalization=F) {
    return(list(
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
		length_normalization=length_normalization))
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