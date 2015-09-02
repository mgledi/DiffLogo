library(cba)
source("R/alphabet.R");
source("R/preconditions.R");
source("R/stackHeights.R");
source("R/baseDistrs.R");
source("R/utilities.R");
source("R/seqLogo.R");
source("R/diffSeqLogo.R");

get.pwm<-function(file,N){
	x<-read.table(file,stringsAsFactors=F)
	
	
	bs<-x[ order(x[,3],decreasing = T),1 ]
	
	
	bs<-bs[1:N]
	
	len<-nchar(bs[1]);
	
	alph<-c(A=1,C=2,G=3,T=4)
	mat<-matrix(alph[unlist(strsplit(x = bs,split = ""))],ncol=len,byrow=T)
	
	pwm<-apply(mat,2,function(a){tab<-table( factor(x = a, levels = 1:4) ); tab/sum(tab)})
	pwm
}

pwm.mad<-get.pwm("inst/extdata/alignments/Mad.txt",1000)
pwm.max<-get.pwm("inst/extdata/alignments/Max.txt",1000)
pwm.myc<-get.pwm("inst/extdata/alignments/Myc.txt",1000)


pdf("ebox.pdf",width=16/1.5, height=10/1.5,compress=T);
    diffLogoTable(list(Mad=pwm.mad,Max=pwm.max,Myc=pwm.myc))
dev.off();

