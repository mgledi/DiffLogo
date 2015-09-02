library(cba)
source("R/alphabet.R");
source("R/preconditions.R");
source("R/stackHeights.R");
source("R/baseDistrs.R");
source("R/utilities.R");
source("R/seqLogo.R");
source("R/diffSeqLogo.R");

# ASN BSP
PWMs = list()
lines = readLines(file("inst/extdata/alignments/F-box_bacteria.seq.fa",open="r"));
PWMs[['bacteria']]=getPwmFromAlignment(lines[grep("^[^>]",lines)],alphabet=ASN,pseudoCount=0);
PWMs[['bacteria']][,48]=1/20
lines = readLines(file("inst/extdata/alignments/F-box_fungi.seq.fa",open="r"));
PWMs[['fungi']]=getPwmFromAlignment(lines[grep("^[^>]",lines)],alphabet=ASN,pseudoCount=0);
lines = readLines(file("inst/extdata/alignments/F-box_metazoa.seq.fa",open="r"));
PWMs[['metazoa']]=getPwmFromAlignment(lines[grep("^[^>]",lines)],alphabet=ASN,pseudoCount=0);
lines = readLines(file("inst/extdata/alignments/F-box_viridiplantae.seq.fa",open="r"));
PWMs[['viridiplantae']]=getPwmFromAlignment(lines[grep("^[^>]",lines)],alphabet=ASN,pseudoCount=0);

pdf("fboxMotifs.pdf",width=6, height=4)
seqLogo(PWMs[['bacteria']],alphabet=ASN)
seqLogo(PWMs[['fungi']],alphabet=ASN)
seqLogo(PWMs[['metazoa']],alphabet=ASN)
seqLogo(PWMs[['viridiplantae']],alphabet=ASN)
dev.off();

pdf("fbox_4_long.pdf",width=6*16/10,height=6,compress=T); 
#png("fbox_4_long.png",width=800*16/10, height=800);
diffLogoTable(PWMs,alphabet=ASN,ymin=0,ymax=0); 
dev.off();

png("fbox_4_long_absDiff.png",width=800*16/10, height=800);
diffLogoTable(PWMs,alphabet=ASN,stackHeight=sumOfAbsICDifferences); 
dev.off();
