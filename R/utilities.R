
##' Creates a matrix-representation of a PWM from a set of sequences
##'
##' @title Create PWM from alignment
##' @param alignment a vector or list of sequences each with equal length
##' @param alphabet of type Alphabet
##' @param pseudoCount the number of pseudo-observations for each character in the alphabet
##' @return PWM as matrix
##' @export
##' @author Hendrik Treutler
##' @examples
##' motif_folder= "extdata/alignments"
##' motif_name = "calamodulin_1"
##' fileName = paste(motif_folder,"/",motif_name,".txt",sep="")
##' file = system.file(fileName, package = "DiffLogo")
##' motif = getPwmFromAlignment(readLines(file), ASN, 1)
##' seqLogo(pwm = motif, alphabet=ASN)
getPwmFromAlignment = function(alignment, alphabet, pseudoCount) {
  alphabetSize = alphabet$size
  alignmentLength = nchar(alignment[[1]])
  numberOfSequences = length(alignment)
  
  pwm = matrix(nrow = alphabetSize, ncol = alignmentLength)
  colnames(pwm) = 1:alignmentLength
  rownames(pwm) = alphabet$chars
  

  for(posIdx in 1:alignmentLength){
    for(charIdx in 1:alphabetSize){
      pwm[charIdx, posIdx] = (length(which(substr(alignment, posIdx, posIdx) == alphabet$chars[[charIdx]])) + pseudoCount) / (numberOfSequences + pseudoCount * alphabetSize)
    }
    pwm[, posIdx] = pwm[, posIdx] / sum(pwm[, posIdx]);
  }
  
  return(pwm);
}
