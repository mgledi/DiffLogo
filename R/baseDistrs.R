
##' TODO
##'
##' @title TODO
##' @param p1 TODO
##' @param p2 TODO
##' @return TODO
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

##' TODO
##'
##' @title TODO
##' @param p1 TODO
##' @param p2 TODO
##' @return TODO
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

##' TODO
##'
##' @title TODO
##' @param p TODO
##' @return TODO
##' @export
##' @author Martin Nettling
probabilities = function(p) {
  return(p);
}
