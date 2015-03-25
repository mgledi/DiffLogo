library(cba)
source("R/alphabet.R");
source("R/preconditions.R");
source("R/stackHeights.R");
source("R/baseDistrs.R");
source("R/utilities.R");
source("R/seqLogo.R");
source("R/diffSeqLogo.R");

## import PWMs
motif_folder = "inst/pwm"
motif_names = c("HepG2","MCF7","HUVEC","ProgFib","NHEK","K562","HeLa-S3","H1-hESC","GM12878")
motif_optimal_order=c(3,1,7,9,6,5,4,2,8);
motifs = list()

for (name in motif_names) {
  fileName = paste(motif_folder,"/",name,".txt",sep="")
  motifs[[name]] = as.matrix(read.delim(fileName,header=F))
}


png("CTCF_9.png",width=1600,height=1000);
	diffLogoTable(motifs,ratio=16/10);
dev.off();

png("CTCF_6.png",width=1600,height=1000);
	diffLogoTable(motifs[motif_optimal_order[c(1,2,3,7,8,9)]],ratio=16/10);
dev.off();

png("CTCF_5.png",width=1600,height=1000);
	diffLogoTable(motifs[motif_optimal_order[c(1,2,3,8,9)]],ratio=16/10);
dev.off();

png("CTCF_4.png",width=1600,height=1000);
	diffLogoTable(motifs[motif_optimal_order[c(1,2,8,9)]],ratio=16/10);
dev.off();

pdf("CTCF_5.pdf",width=16,height=10);
	diffLogoTable(motifs[motif_optimal_order[c(1,2,3,8,9)]],ratio=16/10);
dev.off();
pdf("HUVEC_vs_HeLa-S3.pdf",width=6,height=3);
diffLogoFromPwm(motifs[['HUVEC']],motifs[['HeLa-S3']])
dev.off();

pdf("CTCFseqLogos.pdf",width=6, height=4)
seqLogo(motifs[["HepG2"]])
seqLogo(motifs[["MCF7"]])
seqLogo(motifs[["HUVEC"]])
seqLogo(motifs[["ProgFib"]])
seqLogo(motifs[["NHEK"]])
seqLogo(motifs[["K562"]])
seqLogo(motifs[["HeLa-S3"]])
seqLogo(motifs[["H1-hESC"]])
seqLogo(motifs[["GM12878"]])
dev.off();


# ASN BSP
pwm1=getPwmFromAlignment(readLines(file("inst/alignments/calamodulin_1.txt",open="r")),alphabet=ASN,pseudoCount=0);
pwm2=getPwmFromAlignment(readLines(file("inst/alignments/calamodulin_2.txt",open="r")),alphabet=ASN,pseudoCount=0);

pdf("asnLogos.pdf",width=6, height=4)
seqLogo(pwm1,alphabet=ASN)
seqLogo(pwm2,alphabet=ASN)
dev.off();

pdf("asn_bsp3.pdf",width=6,height=4); diffLogoFromPwm(pwm1,pwm2,alphabet=ASN,ymin=0,ymax=0); dev.off();
