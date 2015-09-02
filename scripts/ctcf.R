library(cba)
source("R/alphabet.R");
source("R/preconditions.R");
source("R/stackHeights.R");
source("R/baseDistrs.R");
source("R/utilities.R");
source("R/seqLogo.R");
source("R/diffSeqLogo.R");

## import PWMs
motif_folder = "inst/extdata/pwm"
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

png("CTCF_all_suboptimal_leafordering.png",width=1600,height=1000);
	diffLogoTable(rev(motifs),ratio=16/10);
dev.off()


pdf("CTCF_all.pdf",width=8*16/10,height=8,compress=T); 
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


e1a = c(1,0,0,0)
e1b = c(0,0,0,1)
e2a = c(0.5,0,0,0.5)
e2b = c(0,0.5,0.5,0)
e3a = c(0.5,0,0,0.5)
e3b = c(0.5,0.0,0.5,0)
e4a = c(0.8,0,0,0.2)
e4b = c(0.2,0,0,0.8)

pdf("tab1.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= shannonDivergence , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-0.5, ymax=0.5  )
diffLogoFromPwm(e2a,e2b,stackHeight= shannonDivergence , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-0.5, ymax=0.5  )
diffLogoFromPwm(e3a,e3b,stackHeight= shannonDivergence , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-0.5, ymax=0.5  )
diffLogoFromPwm(e4a,e4b,stackHeight= shannonDivergence , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-0.5, ymax=0.5  )
dev.off();

pdf("tab2.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= sumOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-2, ymax=2  )
diffLogoFromPwm(e2a,e2b,stackHeight= sumOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-2, ymax=2  )
diffLogoFromPwm(e3a,e3b,stackHeight= sumOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-2, ymax=2  )
diffLogoFromPwm(e4a,e4b,stackHeight= sumOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-2, ymax=2  )
dev.off();

pdf("tab3.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= lossOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-100, ymax=100 )
diffLogoFromPwm(e2a,e2b,stackHeight= lossOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-100, ymax=100 )
diffLogoFromPwm(e3a,e3b,stackHeight= lossOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-100, ymax=100 )
diffLogoFromPwm(e4a,e4b,stackHeight= lossOfAbsICDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-100, ymax=100 )
dev.off();

pdf("tab4.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-1, ymax=1  )
diffLogoFromPwm(e2a,e2b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-1, ymax=1  )
diffLogoFromPwm(e3a,e3b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-1, ymax=1  )
diffLogoFromPwm(e4a,e4b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=normalizedDifferenceOfProbabilities,ymin=-1, ymax=1  )
dev.off();

pdf("tab5.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= shannonDivergence , baseDistribution=differenceOfICs,ymin=-.5, ymax=.5 )
diffLogoFromPwm(e2a,e2b,stackHeight= shannonDivergence , baseDistribution=differenceOfICs,ymin=-.5, ymax=.5 )
diffLogoFromPwm(e3a,e3b,stackHeight= shannonDivergence , baseDistribution=differenceOfICs,ymin=-.5, ymax=.5 )
diffLogoFromPwm(e4a,e4b,stackHeight= shannonDivergence , baseDistribution=differenceOfICs,ymin=-.5, ymax=.5 )
dev.off();

pdf("tab6.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= sumOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-2, ymax=2  )
diffLogoFromPwm(e2a,e2b,stackHeight= sumOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-2, ymax=2  )
diffLogoFromPwm(e3a,e3b,stackHeight= sumOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-2, ymax=2  )
diffLogoFromPwm(e4a,e4b,stackHeight= sumOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-2, ymax=2  )
dev.off();

pdf("tab7.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= lossOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-100, ymax=100 )
diffLogoFromPwm(e2a,e2b,stackHeight= lossOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-100, ymax=100 )
diffLogoFromPwm(e3a,e3b,stackHeight= lossOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-100, ymax=100 )
diffLogoFromPwm(e4a,e4b,stackHeight= lossOfAbsICDifferences , baseDistribution=differenceOfICs,ymin=-100, ymax=100 )
dev.off();

pdf("tab8.pdf",width=3,height=3,compress=T);
diffLogoFromPwm(e1a,e1b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=differenceOfICs,ymin=-1, ymax=1  )
diffLogoFromPwm(e2a,e2b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=differenceOfICs,ymin=-1, ymax=1  )
diffLogoFromPwm(e3a,e3b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=differenceOfICs,ymin=-1, ymax=1  )
diffLogoFromPwm(e4a,e4b,stackHeight= sumOfAbsProbabilityDifferences , baseDistribution=differenceOfICs,ymin=-1, ymax=1  )
dev.off();
