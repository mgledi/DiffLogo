
##' probability differences normalized by the sum of absolute probability differences for the given pair of probability vectors
##'
##' @title normalized probability differences
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return a vector with one result for each symbol
##' @export
##' @author Martin Nettling
normalizedDifferenceOfProbabilities = function(p1,p2) {
    # if p is identical to q, return uniform distribution
    if(all(p1==p2)) {
        return( rep(1/length(p1),length(p1)) )
    }
    tmp = p1-p2;
    return( tmp / sum(abs(tmp)));
}

##' information content differences normalized by the sum of absolute information content differences for the given pair of probability vectors
##'
##' @title normalized information content differences
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return a vector with one result for each symbol
##' @export
##' @author Martin Nettling
differenceOfICs = function(p1,p2) {
    # if p is identical to q, return uniform distribution
    if(all(p1==p2)) {
        return( rep(1/length(p1),length(p1)) )
    }
    H1 = informationContent(p1)$height;
    H2 = informationContent(p1)$height;
    hts1 = H1*p1;
    hts2 = H2*p2
    diff = hts1 - hts2
    return( diff / sum(abs(diff)));
}

##' the given probabilities
##'
##' @title probabilities
##' @param p probability vector representing the symbol distribution
##' @return the given vector
##' @export
##' @author Martin Nettling
probabilities = function(p) {
  return(p);
}
