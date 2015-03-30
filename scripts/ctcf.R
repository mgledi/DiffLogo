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


png("CTCF_5.png",width=1600,height=1000);
	diffLogoTable(motifs[motif_optimal_order[c(1,2,3,8,9)]],ratio=16/10);
dev.off();

png("CTCF_5_absDiff.png",width=1600,height=1000);
	diffLogoTable(motifs[motif_optimal_order[c(1,2,3,8,9)]],ratio=16/10,stackHeight=sumOfAbsICDifferences);
dev.off();
