
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

getAlphabetFromSequences <- function(sequences){
  characters <- unique(unlist(strsplit(sequences, split = "")))
  alphabet <- getAlphabetFromCharacters(characters)
  return(alphabet)
}
getAlphabetFromCharacters <- function(characters){
  chars <- paste(sort(characters), collapse = "")
  if(grepl(pattern = "^A?C?G?T?$", x = chars)){
    return(DNA)
  } else if(grepl(pattern = "^A?C?G?U?$", x = chars)){
    return(DNA)
  } else
    return(ASN)
}

getPwmFromFile <- function(filename){
  extension <- tolower(file_ext(filename))
  #fileContent <- readLines(con = filename)
  
  pwm <- NULL
  if(extension == "fa" || extension == "fasta")
    pwm <- getPwmFromFastaFile(filename)
  if(extension == "txt" || extension == "text" || extension == "al" || extension == "alignment")
    pwm <- getPwmFromAlignmentFile(filename)
  if(extension == "pwm")
    pwm <- getPwmFromPwmFile(filename)
  if(extension == "pfm")
    pwm <- getPwmFromPfmOrJasperFile(filename)
  if(extension == "motif")
    pwm <- getPwmFromHomerFile(filename)
  
  if(is.null(pwm))
    stop(paste("The file extension", extension, "is not supported."))
  
  return(pwm)
}

getPwmFromFastaFile = function(filename) {
    connection = file(filename ,open="r");
    lines = as.vector(read.delim(connection)[,1]);
    close(connection);
    lines = lines[grep("^[^>]",lines)]
    lines = lines[sapply(lines,nchar) > 0]; # remove empty lines
    lines = toupper(lines);
    
    pwm <- getPwmFromAlignment(alignment = lines, alphabet = getAlphabetFromSequences(lines), pseudoCount = 0)
    return(pwm);
}

getPwmFromAlignmentFile = function(filename) {
    connection = file(filename ,open="r");
    lines = as.vector(read.delim(connection)[,1]);
    close(connection);
    lines = lines[sapply(lines,nchar) > 0]; # remove empty lines
    lines = toupper(lines);
    
    pwm <- getPwmFromAlignment(alignment = lines, alphabet = getAlphabetFromSequences(lines), pseudoCount = 0)
    return(pwm);
}

getPwmFromPwmFile = function(filename) {
    lines = readLines(filename);
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=F));
    close(tc);
    pwm = pwm / apply(pwm,2,sum); # normalize pwm
    return(pwm);
}

getPwmFromPfmOrJasperFile = function(filename) {
    lines = readLines(filename);
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=F));
    close(tc);
    pwm = pwm / apply(pwm,2,sum); # normalize pwm
    return(pwm);
}

getPwmFromHomerFile = function(filename) {
    # First read lines
    lines = readLines(filename);
    lines = lines[grep("^[^>]",lines)]
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep = " ",header=F));
    close(tc);

    # transpose
    pwm = t(pwm); 
    pwm = pwm / apply(pwm,2,sum); # normalize pwm
    return(pwm);
}
