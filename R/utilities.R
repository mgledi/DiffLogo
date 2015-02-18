
##' Creates a matrix-representation of a PWM from a set of sequences
##'
##' @title Create PWM from alignment
##' @param alignment a vector or list of sequences each with equal length
##' @param alphabet of type Alphabet
##' @param pseudoCount the number of pseudo-observations for each character in the alphabet
##' @return PWM as matrix
##' @export
##' @author Hendrik Treutler
getPwmFromAlignment = function(alignment, alphabet, pseudoCount) {
  alphabetSize = alphabet$size
  alignmentLength = nchar(alignment[[1]])
  numberOfSequences = length(alignment)
  
  pwm = matrix(nrow = alphabetSize, ncol = alignmentLength)
  colnames(pwm) = 1:alignmentLength
  rownames(pwm) = alphabet$chars
  
  for(charIdx in 1:alphabetSize){
    for(posIdx in 1:alignmentLength){
      pwm[charIdx, posIdx] = (length(which(substr(alignment, posIdx, posIdx) == alphabet$chars[[charIdx]])) + pseudoCount) / (numberOfSequences + pseudoCount * alphabetSize)
    }
  }
  
  return(pwm);
}