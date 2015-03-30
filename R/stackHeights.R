
##' TODO
##'
##' @title TODO
##' @param p1 TODO
##' @param p2 TODO
##' @return an object consisting of height a ylab
##' @export
##' @author Martin Nettling
sumOfAbsICDifferences = function(p1,p2) {
    H1 = informationContent(p1)$height
    H2 = informationContent(p2)$height
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

##' TODO
##'
##' @title TODO
##' @param p TODO
##' @param q TODO
##' @return an object consisting of height a ylab
##' @export
##' @author Martin Nettling
shannonDivergence = function(p,q) {
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p==q)) {
        obj$height=0;   
    } else {
        m = (p + q) / 2
        obj$height=0.5*sum( p * (log2(p) - log2(m)),na.rm=T) + 0.5*sum( q * (log2(q) - log2(m)),na.rm=T);
    }
    obj$ylab="JS divergence";
    return(obj);
}

##' TODO
##'
##' @title TODO
##' @param p TODO
##' @return an object consisting of height a ylab
##' @export
##' @author Martin Nettling
sumProbabilities = function (p) {
    obj=list();
    obj$height=sum(p);
    obj$ylab="Probability";
    return(obj);
}

##' TODO
##'
##' @title TODO
##' @param p TODO
##' @return an object consisting of height a ylab
##' @export
##' @author Martin Nettling
informationContent = function (p) {
    obj=list();
    ic = log2(length(p)) + sum(x * log2(x),na.rm=T);
    obj$height=ic;
    obj$ylab="Information Content [bits]";
    return(obj);
}
