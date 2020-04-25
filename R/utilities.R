
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

##' returns the alphabet which fits to the given sequences
##' 
##' @title returns the alphabet which fits to the given sequences
##' @param sequences a character vector of sequences
##' @return an alphabet of type Alphabet
##' 
##' @export
##' @examples 
##' alphabet = getAlphabetFromSequences("AACCGGTT")
getAlphabetFromSequences <- function(sequences){
    characters <- unique(unlist(strsplit(sequences, split = "")))
    alphabet <- getAlphabetFromCharacters(characters)
    return(alphabet)
}
##' returns the alphabet which fits to the given characters
##' 
##' @title returns the alphabet which fits to the given characters
##' @param characters a character vector of characters
##' @return an alphabet of type Alphabet
##' 
##' @export
##' @examples 
##' alphabet = getAlphabetFromSequences(c("A","A","C","C","G","G","T","T"))
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
##' Generates a pwm from a file of different formats. 
##' Supported formats are FASTA files (.fa, .fasta), alignment files (.txt, .text, .al, .alignment), PWM files (.pwm), JASPAR / Position Frequency Matrix files (.pfm), and homer files (.motif).
##' 
##' @title generates a pwm from a file of different formats
##' @param filename the file
##' @return a pwm
##' 
##' @export
##' @examples 
##' fileName = "extdata/pwm/H1-hESC.pwm"
##' file = system.file(fileName, package = "DiffLogo")
##' pwm = getPwmFromFile(file)
getPwmFromFile <- function(filename){
  extension <- tolower(tools::file_ext(filename))
  
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

##' extracts the sequences from a FASTA file
##' 
##' @title extracts the sequences from a FASTA file
##' @param filename the FASTA file
##' @return a vector of sequences
##' @importFrom utils read.delim
##' @export
##' @examples 
##' fileName = "extdata/alignments/F-box_bacteria.seq.fa"
##' file = system.file(fileName, package = "DiffLogo")
##' sequences = getSequencesFromFastaFile(file)
getSequencesFromFastaFile = function(filename) {
    connection = file(filename ,open="r");
    lines = as.vector(read.delim(connection)[,1]);
    close(connection);
    lines = lines[grep("^[^>]",lines)]
    lines = lines[sapply(lines,nchar) > 0]; # remove empty lines
    lines = toupper(lines);
    return(lines);
}


##' generates a pwm from a FASTA file
##' 
##' @title generates a pwm from a FASTA file
##' @param filename the FASTA file
##' @param alphabet the desired alphabet of type Alphabet
##' @return a pwm
##' 
##' @export
##' @examples 
##' fileName = "extdata/alignments/F-box_bacteria.seq.fa"
##' file = system.file(fileName, package = "DiffLogo")
##' pwm = getPwmFromFastaFile(file)
getPwmFromFastaFile = function(filename,alphabet=NULL) {
    lines = getSequencesFromFastaFile(filename);
    pwm <- getPwmFromAlignment(alignment = lines, alphabet = alphabet)
    return(pwm);
}

##' extracts the sequences from an alignment file
##' 
##' @title extracts the sequences from an alignment file
##' @param filename the alignment file
##' @return a vector of sequences
##' @importFrom utils read.delim
##' @export
##' @examples 
##' fileName = "extdata/alignments/calamodulin_1.txt"
##' file = system.file(fileName, package = "DiffLogo")
##' sequences = getSequencesFromAlignmentFile(file)
getSequencesFromAlignmentFile = function(filename) {
    connection = file(filename ,open="r");
    lines = as.vector(read.delim(connection)[,1]);
    close(connection);
    lines = lines[sapply(lines,nchar) > 0]; # remove empty lines
    lines = toupper(lines);
    return(lines);
}

##' generates a pwm from an alignment file
##' 
##' @title generates a pwm from an alignment file
##' @param filename the alignment file
##' @param alphabet the desired alphabet of type Alphabet
##' @return a pwm
##' 
##' @export
##' @examples 
##' fileName = "extdata/alignments/calamodulin_1.txt"
##' file = system.file(fileName, package = "DiffLogo")
##' pwm = getPwmFromAlignmentFile(file)
getPwmFromAlignmentFile = function(filename,alphabet=NULL) {
    lines = getSequencesFromAlignmentFile(filename);
    pwm <- getPwmFromAlignment(alignment = lines, alphabet = alphabet)
    return(pwm);
}

##' generates a pwm from a pwm file
##' 
##' @title generates a pwm from a pwm file
##' @param filename the pwm file
##' @return a pwm
##' @importFrom utils read.delim
##' @export
##' @examples 
##' fileName = "extdata/pwm/H1-hESC.pwm"
##' file = system.file(fileName, package = "DiffLogo")
##' pwm = getPwmFromPwmFile(file)
getPwmFromPwmFile = function(filename) {
    lines = readLines(filename);
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=FALSE));
    close(tc);
    pwm = normalizePWM(pwm);
    return(pwm);
}

##' generates a pwm from a jaspar file
##' 
##' @title generates a pwm from a jaspar file
##' @param filename the jaspar file
##' @return a pwm
##' @importFrom utils read.delim
##' @export
##' @examples 
##' fileName = "extdata/pfm/ctcf_jaspar.pfm"
##' file = system.file(fileName, package = "DiffLogo")
##' pwm = getPwmFromPfmOrJasparFile(file)
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
    pwm = as.matrix(read.delim(tc, sep=" ", header=FALSE));
    close(tc);
    pwm = normalizePWM(pwm);
    return(pwm);
}

##' generates a pwm from a homer file
##' 
##' @title generates a pwm from a homer file
##' @param filename the homer file
##' @return a pwm
##' @importFrom utils read.delim
##' @export
##' @examples 
##' fileName = "extdata/homer/CTCF_Zf_CD4.motif"
##' file = system.file(fileName, package = "DiffLogo")
##' pwm = getPwmFromHomerFile(file)
getPwmFromHomerFile = function(filename) {
    # First read lines
    lines = readLines(filename);
    lines = lines[grep("^[^>]",lines)]
    # replace all whitespaces by one " "
    lines = gsub("\\s+", " ", lines);
    tc = textConnection(lines);
    pwm = as.matrix(read.delim(tc, sep=" ", header=FALSE));
    close(tc);

    # transpose
    pwm = t(pwm); 
    pwm = normalizePWM(pwm);
    return(pwm);
}

##' normalizes the given pwm to column-sums of 1.0
##' 
##' @title normalizes the given pwm
##' @param pwm a pwm
##' @return a normalized pwm
##' 
##' @export
##' @examples 
##' pwm = matrix(1:40, nrow = 4, dimnames = list(c("A","C","G","T"), 1:10))
##' pwm = normalizePWM(pwm)
normalizePWM = function(pwm) {
    for( i in 1:ncol(pwm)) {
        pwm[,i] = pwm[,i] / sum(pwm[,i]);
    }
    return(pwm);
}

