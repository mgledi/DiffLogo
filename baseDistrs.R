##############
# returns a vector of size 4
normalizedDifferenceOfProbabilities = function(p1,p2) {
    tmp = p1-p2;
    return( tmp / sum(abs(tmp)));
}


##############
# returns a vector of size 4
differenceOfICs = function(p1,p2) {
    H1 = pwm2ic_col(as.matrix(p1));
    H2 = pwm2ic_col(as.matrix(p2))
    hts1 = H1*p1;
    hts2 = H2*p2
    diff = hts1 - hts2
    return( diff / sum(abs(diff)));
}