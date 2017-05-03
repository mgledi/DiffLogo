
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
getPwmFromAlignment = function(alignment, alphabet=NULL, pseudoCount=0) {
    alignment <- gsub("\\s+", "", alignment)
    alignment <- alignment[!grepl(pattern = "^$", x = alignment)]
    
    if(is.null(alphabet))
      alphabet <- getAlphabetFromSequences(alignment)
    
    alphabetSize = alphabet$size
    alignmentLength = nchar(alignment[[1]])
    numberOfSequences = length(alignment)
    
    if(length(unique(unlist(lapply(X = alignment, FUN = nchar)))) > 1)
      stop("Alignment comprises sequences with different lengths.")
    
    pwm = matrix(nrow = alphabetSize, ncol = alignmentLength)
    colnames(pwm) = 1:alignmentLength
    rownames(pwm) = alphabet$chars

    for(posIdx in 1:alignmentLength){
        for(charIdx in 1:alphabetSize){
            pwm[charIdx, posIdx] = (length(which(substr(alignment, posIdx, posIdx) == alphabet$chars[[charIdx]])) + pseudoCount) / (numberOfSequences + pseudoCount * alphabetSize)
        }
        if(sum(pwm[, posIdx]) == 0) {
            # nothing was observed and no pseudoCount was set
            pwm[, posIdx] = rep(1/alphabetSize,alphabetSize);
        } else {
            # normalize pwm
            pwm[, posIdx] = pwm[, posIdx] / sum(pwm[, posIdx]);            
        }
    }
    return(pwm);
}

##' @export
getAlphabetFromSequences <- function(sequences){
    characters <- unique(unlist(strsplit(sequences, split = "")))
    alphabet <- getAlphabetFromCharacters(characters)
    return(alphabet)
}
##' @export
getAlphabetFromCharacters <- function(characters){
    chars <- paste(sort(characters), collapse = "")
  
    dnaRegEx <- paste("^\\-?", paste(sort(DNA$chars), "?", sep = "", collapse = ""), "$", sep = "")
    rnaRegEx <- paste("^\\-?", paste(sort(RNA$chars), "?", sep = "", collapse = ""), "$", sep = "")
    asnRegEx <- paste("^\\-?", paste(sort(ASN$chars), "?", sep = "", collapse = ""), "$", sep = "")
  
    if(grepl(pattern = dnaRegEx, x = chars)){
        return(DNA)
    } else if(grepl(pattern = rnaRegEx, x = chars)){
        return(DNA)
    } else if(grepl(pattern = asnRegEx, x = chars)){
        return(ASN)
    } else{
        warning(paste("Unrecognized alphabet (not DNA, RNA, or ASN):", chars))
        return(FULL_ALPHABET)
    }
}
##' @export
getPwmFromFile <- function(filename){
  extension <- tolower(file_ext(filename))
  
  pwm <- NULL
  error <- NULL
  tryCatch(
    {
      if(extension == "fa" || extension == "fasta")
        pwm <- getPwmFromFastaFile(filename)
      if(extension == "txt" || extension == "text" || extension == "al" || extension == "alignment")
        pwm <- getPwmFromAlignmentFile(filename)
      if(extension == "pwm")
        pwm <- getPwmFromPwmFile(filename)
      if(extension == "pfm")
        pwm <- getPwmFromPfmOrJasparFile(filename)
      if(extension == "motif")
        pwm <- getPwmFromHomerFile(filename)
    }, 
    error = function(e) {
      error <- e
    }
  )
  
  if(!is.null(error)) {
    stop(paste("Could not parse file", basename(filename), ". Error:", error))
  }
  if(is.null(pwm)) {
    stop(paste("The file extension", extension, "of file", basename(filename),"is not supported."))
  }
  
  return(pwm)
}

##' @export
getSequencesFromFastaFile = function(filename) {
    connection = file(filename ,open="r");
    lines = as.vector(read.delim(connection)[,1]);
    close(connection);
    lines = lines[grep("^[^>]",lines)]
    lines = lines[sapply(lines,nchar) > 0]; # remove empty lines
    lines = toupper(lines);
    return(lines);
}


##' @export
getPwmFromFastaFile = function(filename,alphabet=NULL) {
    lines = getSequencesFromFastaFile(filename);
    pwm <- getPwmFromAlignment(alignment = lines, alphabet = alphabet)
    return(pwm);
}

##' @export
getSequencesFromAlignmentFile = function(filename) {
    connection = file(filename ,open="r");
    lines = as.vector(read.delim(connection)[,1]);
    close(connection);
    lines = lines[sapply(lines,nchar) > 0]; # remove empty lines
    lines = toupper(lines);
    return(lines);
}

##' @export
getPwmFromAlignmentFile = function(filename,alphabet=NULL) {
    lines = getSequencesFromAlignmentFile(filename);
    pwm <- getPwmFromAlignment(alignment = lines, alphabet = alphabet)
    return(pwm);
}
##' @export
getPwmFromPwmFile = function(filename) {
    lines = readLines(filename);
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=F));
    close(tc);
    pwm = normalizePWM(pwm);
    return(pwm);
}
##' @export
getPwmFromPfmOrJasparFile = function(filename) {
    lines = readLines(filename);
    lines = lines[grep("^[^>]",lines)]
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    # remove ACGT
    lines = gsub("[ACGT]", "", lines);    
    # remove opening gap
    lines = gsub("\\s+\\[", "", lines);
    lines = gsub("\\s+\\]", "", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=F));
    close(tc);
    pwm = normalizePWM(pwm);
    return(pwm);
}
##' @export
getPwmFromHomerFile = function(filename) {
    # First read lines
    lines = readLines(filename);
    lines = lines[grep("^[^>]",lines)]
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=F));
    close(tc);

    # transpose
    pwm = t(pwm); 
    pwm = normalizePWM(pwm);
    return(pwm);
}

##' @export
normalizePWM = function(pwm) {
    for( i in 1:ncol(pwm)) {
        pwm[,i] = pwm[,i] / sum(pwm[,i]);
    }
    return(pwm);
}

