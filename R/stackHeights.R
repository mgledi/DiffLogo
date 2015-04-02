
##' the sum of absolute information content differences for the given pair of probability vectors
##'
##' @title sum of absolute information content differences
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return an object consisting of height and ylab
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

##' the shannon divergence for the given pair of probability vectors
##'
##' @title shannon divergence
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return an object consisting of height and ylab
##' @export
##' @author Martin Nettling
shannonDivergence = function(p1,p2) {
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p1==p2)) {
        obj$height=0;   
    } else {
        m = (p1 + p2) / 2
        obj$height=0.5*sum( p1 * (log2(p1) - log2(m)),na.rm=T) + 0.5*sum( p2 * (log2(p2) - log2(m)),na.rm=T);
    }
    obj$ylab="JS divergence";
    return(obj);
}

##' the sum of probabilities for the given probability vector, i.e. 1.0
##'
##' @title sum of probabilities, i.e. 1.0
##' @param p probability vector representing the symbol distribution
##' @return an object consisting of height and ylab
##' @export
##' @author Martin Nettling
sumProbabilities = function (p) {
    obj=list();
    obj$height=sum(p);
    obj$ylab="Probability";
    return(obj);
}

##' the information content for the given probability vector
##'
##' @title information content
##' @param p probability vector representing the symbol distribution
##' @return an object consisting of height a ylab
##' @export
##' @author Martin Nettling
informationContent = function (p) {
    obj=list();
    ic = log2(length(p)) + sum(p * log2(p),na.rm=T);
    obj$height=ic;
    obj$ylab="Information Content [bits]";
    return(obj);
}
