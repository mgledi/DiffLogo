##############
# sums up the absolute differences of corresponding probability weighted ICs
# returns an object consisting of height a ylab
sumOfAbsICDifferences = function(p1,p2) {
    H1 = pwm2ic_col(as.matrix(p1))
    H2 = pwm2ic_col(as.matrix(p2))
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p1==p2)) {
        obj$height=0;   
    } else {
        obj$height=sum( abs(H1*p1 - H2*p2));
    }
    obj$ylab="Information Content [bits]";
    return(obj);
}

##############
# sums up the absolute differences of corresponding probability weighted ICs
# returns an object consisting of height and a ylab
shannonDivergence = function(p,q) {
    base=log2(length(p))
    m = (p + q) / 2
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p==q)) {
        obj$height=0;   
    } else {
        obj$height=0.5*sum( p * (log2(p)/base - log2(m)/base)) + 0.5*sum( q * (log2(q)/base - log2(m)/base));
    }
    obj$ylab="JS divergence";
    return(obj);
}
