

##' the sum of absolute probability differences for the given pair of probability vectors
##'
##' @title sum of absolute probability differences
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return an object consisting of height and ylab
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##' 
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##' 
##' diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2, stackHeight = sumOfAbsProbabilityDifferences)
sumOfAbsProbabilityDifferences = function(p1,p2) {
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p1==p2)) {
        obj$height=0;   
    } else {
        obj$height=sum( abs(p1 - p2));
    }
    obj$ylab="Probability";
    return(obj);
}

##' the sum of absolute information content differences for the given pair of probability vectors
##'
##' @title sum of absolute information content differences
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return an object consisting of height and ylab
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##' 
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##' 
##' diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2, stackHeight = sumOfAbsICDifferences)
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

##' the change of information content for the given probability vectors
##'
##' @title the change of information content
##' @param p1 probability vector representing the first symbol distribution
##' @param p2 probability vector representing the second symbol distribution
##' @return an object consisting of height and ylab
##' @export
##' @author Martin Nettling
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##' 
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##' 
##' diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2, stackHeight = lossOfAbsICDifferences)
lossOfAbsICDifferences = function(p1,p2) {
    H1 = informationContent(p1)$height
    H2 = informationContent(p2)$height
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p1==p2)) {
        obj$height=0;   
    } else {
        obj$height=sum( abs(H1*p1 - H2*p2)) / (H1/2 + H2/2) * 100;
    }
    obj$ylab="Change of Information Content [%]";
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
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_names = c("HepG2","MCF7","HUVEC","ProgFib")
##' motifs = list()
##' for (name in motif_names) {
##'   fileName = paste(motif_folder,"/",name,".pwm",sep="")
##'   file = system.file(fileName, package = "DiffLogo")
##'   motifs[[name]] = getPwmFromPwmFile(file)
##' }
##' 
##' pwm1 = motifs[[motif_names[[1]]]]
##' pwm2 = motifs[[motif_names[[2]]]]
##' 
##' diffLogoFromPwm(pwm1 = pwm1, pwm2 = pwm2, stackHeight = shannonDivergence)
shannonDivergence = function(p1,p2) {
    obj=list();
    # if p is identical to q, set height to 0
    if(all(p1==p2)) {
        obj$height=0;   
    } else {
        m = (p1 + p2) / 2
        obj$height=0.5*sum( p1 * (log2(p1) - log2(m)),na.rm=TRUE) + 0.5*sum( p2 * (log2(p2) - log2(m)),na.rm=TRUE);
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
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_name = "HepG2"
##' fileName = paste(motif_folder,"/",motif_name,".pwm",sep="")
##' file = system.file(fileName, package = "DiffLogo")
##' motif = getPwmFromPwmFile(file)
##' seqLogo(pwm = motif, stackHeight = sumProbabilities)
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
##' @examples
##' motif_folder= "extdata/pwm"
##' motif_name = "HepG2"
##' fileName = paste(motif_folder,"/",motif_name,".pwm",sep="")
##' file = system.file(fileName, package = "DiffLogo")
##' motif = getPwmFromPwmFile(file)
##' seqLogo(pwm = motif, stackHeight = informationContent)
informationContent = function (p) {
    obj=list();
    ic = log2(length(p)) + sum(p * log2(p),na.rm=TRUE);
    obj$height=ic;
    obj$ylab="Information Content [bits]";
    return(obj);
}
