##############
# returns a vector of size 4
normalizedDifferenceOfProbabilities = function(p1,p2) {
    # if p is identical to q, return uniform distribution
    if(all(p1==p2)) {
        return( rep(1/length(p1),length(p1)) )
    }
    tmp = p1-p2;
    return( tmp / sum(abs(tmp)));
}

##############
# returns a vector of size 4
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

probabilities = function(p) {
  return(p);
}