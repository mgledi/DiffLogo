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

png("CTCF_all.png",width=1600,height=1000);
	diffLogoTable(motifs,ratio=16/10);
dev.off()

m1 = motifs[["HUVEC"]]
m2 = motifs[["H1-hESC"]]

pdf("configurations.pdf",height=4,width=6);
diffLogoFromPwm(m1,m2,stackHeight= shannonDivergence , baseDistribution=normalizedDifferenceOfProbabilities )
title("shannonDivergence | normalizedDifferenceOfProbabilities");
diffLogoFromPwm(m1,m2,stackHeight= sumOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities )
title("sumOfAbsICDifferences | normalizedDifferenceOfProbabilities");
diffLogoFromPwm(m1,m2,stackHeight= lossOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities )
title("lossOfAbsICDifferences | normalizedDifferenceOfProbabilities");
diffLogoFromPwm(m1,m2,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=normalizedDifferenceOfProbabilities )
title("sumOfAbsProbabilityDifferences | normalizedDifferenceOfProbabilities");


diffLogoFromPwm(m1,m2,stackHeight= shannonDivergence , baseDistribution=differenceOfICs )
title("shannonDivergence | differenceOfICs");
diffLogoFromPwm(m1,m2,stackHeight= sumOfAbsICDifferences , baseDistribution=differenceOfICs )
title("sumOfAbsICDifferences | differenceOfICs");
diffLogoFromPwm(m1,m2,stackHeight= lossOfAbsICDifferences , baseDistribution=differenceOfICs )
title("lossOfAbsICDifferences | differenceOfICs");
diffLogoFromPwm(m1,m2,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=differenceOfICs )
title("sumOfAbsProbabilityDifferences | differenceOfICs");

dev.off();
