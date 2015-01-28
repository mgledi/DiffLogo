##############
# sums up the absolute differences of corresponding probability weighted ICs
# returns an object consisting of height a ylab
sumOfAbsICDifferences = function(p1,p2) {
    H1 = pwm2ic_col(as.matrix(p1))
    H2 = pwm2ic_col(as.matrix(p2))
    obj=list();
    obj$height=sum( abs(H1*p1 - H2*p2));
    obj$ylab="Information Content [bits]";
    return(obj);
}

##############
# sums up the absolute differences of corresponding probability weighted ICs
# returns an object consisting of height a ylab
shannonDivergence = function(p,q) {
    m = (p + q) / 2
    obj=list();
    obj$height=0.5*sum( p * (log2(p) - log2(m))) + 0.5*sum( q * (log2(q) - log2(m)));
    obj$ylab="JS divergence";
    return(obj);
}