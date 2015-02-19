
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




source("R/alphabet.R");
source("R/preconditions.R");
source("R/stackHeights.R");
source("R/baseDistr.R");
source("R/utilies.R");
source("R/seqLogo.R");
source("R/diffSeqLogo.R");

## import PWMs
motif_folder = "inst/pwm"
motif_names = c("HepG2","MCF7","HUVEC","ProgFib","NHEK","K562","HeLa-S3","H1-hESC","GM12878")
motifs = list()

for (name in motif_names) {
  fileName = paste(motif_folder,"/",name,".txt",sep="")
  motifs[[name]] = as.matrix(read.delim(fileName,header=F))
}


source("R/diffSeqLogo.R");
png("leafOrder2.png",width=1600,height=1000);
diffLogoTable(motifs,ratio=16/10);
dev.off();
