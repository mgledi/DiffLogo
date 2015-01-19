# set working folder to this directory
# setwd("/data/m.gleditzsch/workspace_git/comparative-sequence-logo/")
source("./diffSeqLogo.R");

motif_folder = "CTCF_motifs";
motif_names = c("GM12878","HeLa-S3","HepG2","HUVEC","K562","MCF7","NHEK","ProgFib");
marSeqLogo=c(1,1.5,0.1,0.1);
margin = 0.02;


pdf(paste("CTCF_motifs_comparison.pdf",sep=""), width=16, height=16);
#dev.off();
CTCF = list();

for (name in motif_names) {
    file = paste("./",motif_folder,"/",name,".txt",sep="");
    CTCF[[name]] = as.matrix(read.delim(file,header=F));
}

dim = length(CTCF);
#dim = 3;
for ( i in 1:dim) {
    motif_i = motif_names[i];
    for ( k in 1:dim) {
        if( i != k ) {
            motif_k = motif_names[k];
            print(paste("plotting ",motif_i," and ",motif_k));
            par(fig=(c(i-1,i,dim-k,dim-k+1) / dim) * (1-margin) + c(margin,margin,0,0), new=TRUE, mar=marSeqLogo)
            diffSeqLogo(CTCF[[ motif_i ]],CTCF[[ motif_k ]],type=3,showSums=F,sparse=T,ymin=-0.4, ymax=0.4)
        }
    }
}
# add names
par(fig=c(0,1,0,1)* (1-margin) + c(margin,margin,0,0), new=TRUE, mar=c(0,0,0,0))
plot(NA,ylim=c(0,dim),xlim=c(0,dim),,xaxt="n",yaxt="n",xaxs="i",yaxs="i",bty="n") #xaxt="n",yaxt="n",
axis(4, pos=-0.19, at=(1:dim) - 0.5, labels = rev(motif_names[1:dim]), tick = F, cex.axis=1.5)
axis(3, pos=dim-.10, at= (1:dim) - 0.5, labels = motif_names[1:dim], tick = F, cex.axis=1.5)

dev.off();
