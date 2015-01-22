# set working folder to this directory
# setwd("/data/m.gleditzsch/workspace_git/comparative-sequence-logo/")
source("./diffSeqLogo.R");

motif_folder = "CTCF_motifs";
#motif_names = c("GM12878","HeLa-S3","HepG2","HUVEC","K562","MCF7","NHEK","ProgFib");
motif_names = c("HepG2","HUVEC","MCF7","ProgFib");
marSeqLogo=c(1,1.5,0.1,0.1);
margin = 0.03;
ratio=16/10;

ymin=-0.22; ymax=0.22;
#ymin=-0.4; ymax=0.4;

if( !exists("type") ) {
  stop("Please define \"type\" before continuing");
}  



CTCF = list();

pdf(paste("CTCF_motifs.pdf",sep=""), width=6, height=4.5);
for (name in motif_names) {
    file = paste("./",motif_folder,"/",name,".txt",sep="");
    CTCF[[name]] = as.matrix(read.delim(file,header=F));

    seqLogo(CTCF[[name]])
    
}
dev.off();

# TODO: formulate this as function of a list of TFBSs
dim = length(CTCF);

#pdf(paste("CTCF_motifs_comparison_short_",type,".pdf",sep=""), width=dim*2, height=dim*2/ratio);
diffSeqLogoScatterPlot(CTCF,type=type,margin=margin);
#dev.off();
