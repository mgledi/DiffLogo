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
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p==q)) {
        obj$height=0;   
    } else {
        px= p[p>0];
        qx= q[q>0];
        mx = (px + qx) / 2
        obj$height=0.5*sum( px * (log2(px) - log2(mx))) + 0.5*sum( qx * (log2(qx) - log2(mx)));
    }
    obj$ylab="JS divergence";
    return(obj);
}


sumProbabilities = function (p) {
    obj=list();
    obj$height=sum(p);
    obj$ylab="Probability";
    return(obj);
}


informationContent = function (p) {
    obj=list();
    x = p[p>0];
    ic = log2(length(p)) + sum(x * log2(x));
    obj$height=ic;
    obj$ylab="Information Content [bits]";
    return(obj);
}
